#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Vector-friendly fonts for PDF/SVG
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"


def infer_species_from_filename(path: str) -> str:
    """
    Extract species name from filename.
    Example:
        g2400.t1_Nkunapalari.matrix  ->  Nkunapalari
    """
    base = os.path.basename(path)
    m = re.match(r".+\.t1_([A-Za-z]+)\.matrix$", base)
    if not m:
        raise ValueError(f"Filename not in expected format: {base}")
    return m.group(1)


def read_one_matrix(path: str, pos_only: bool = True) -> pd.DataFrame:
    """
    Read a single *.matrix file.
    Returns a DataFrame with index = site, column = mean_alt_prop.
    """
    df = pd.read_csv(path, sep="\t", engine="python", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    if df.shape[1] < 3:
        raise ValueError(f"{path}: expect at least 3 columns (CHROM, POS, samples...)")

    meta = df.iloc[:, :2].copy()
    meta.columns = ["CHROM", "POS"]

    vals = df.iloc[:, 2:].apply(pd.to_numeric, errors="coerce")
    mean_vals = vals.mean(axis=1, skipna=True)

    if pos_only:
        site = meta["POS"].astype(str)
    else:
        site = meta["CHROM"].astype(str) + ":" + meta["POS"].astype(str)

    out = pd.DataFrame({"site": site, "mean": mean_vals})
    out = (
        out.dropna(subset=["mean"], how="all")
           .drop_duplicates(subset=["site"])
           .set_index("site")
    )
    return out


def tidy_species_label(sp: str) -> str:
    """
    Convert species names into pretty labels.
    Example:
        Nkunapalari -> N. kunapalari
    """
    if sp.startswith("N"):
        return "N. " + sp[1:].lower()
    return sp


def plot_one_gene(gene_id: str, matrix_files, ax_tree, ax_heat, cmap, method="average",
                  metric="euclidean", yticksize=10):
    """
    Plot the heatmap and hierarchical clustering tree for a single gene
    onto given axes: ax_tree (left) and ax_heat (right).
    """
    mats = []
    order_sites = None
    for f in matrix_files:
        sp = infer_species_from_filename(f)
        df = read_one_matrix(f, pos_only=True)
        df.columns = [sp]
        if order_sites is None:
            order_sites = df.index.tolist()
        mats.append(df)

    M = pd.concat(mats, axis=1)
    M = M.loc[order_sites, :]

    data_for_dist = np.nan_to_num(M.T.values, nan=0.0)
    if data_for_dist.shape[0] > 1:
        Z = linkage(pdist(data_for_dist, metric=metric), method=method)
        row_leaves = leaves_list(Z)
        species_order = [M.columns[i] for i in row_leaves]
    else:
        Z = None
        species_order = list(M.columns)

    # Reorder rows based on clustering
    M = M.loc[:, species_order].T

    # Draw left dendrogram
    ax = ax_tree
    if Z is not None:
        dendrogram(
            Z,
            orientation="left",
            ax=ax,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="black",
            link_color_func=lambda k: "black",
        )
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_facecolor("white")
    else:
        ax.axis("off")

    # Gene ID label on top
    pos = ax.get_position()
    ax.figure.text(pos.x0, pos.y1 + 0.005, gene_id, fontsize=10,
                   ha="left", va="bottom")

    # Draw heatmap
    axh = ax_heat
    im = axh.imshow(
        M.values,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
    )

    axh.set_yticks(np.arange(M.shape[0]))
    ylabels = [tidy_species_label(s) for s in M.index]
    axh.set_yticklabels(ylabels, fontsize=yticksize)
    axh.yaxis.tick_right()
    axh.yaxis.set_label_position("right")
    for lab in axh.get_yticklabels():
        lab.set_fontstyle("italic")
    axh.set_ylabel("")

    # X-axis ticks = genomic positions
    ncols = M.shape[1]
    step = max(1, ncols // 15)
    xticks = np.arange(0, ncols, step)
    axh.set_xticks(xticks)
    axh.set_xticklabels(
        [M.columns[i] for i in xticks],
        rotation=60,
        ha="right",
        fontsize=8,
    )
    axh.set_xlabel("genome position", fontsize=9)

    return im, axh


def main():
    parser = argparse.ArgumentParser(
        description="Draw multi-gene heatmaps (4 genes per page)."
    )
    parser.add_argument("--indir", required=True,
                        help="directory containing *.matrix files")
    parser.add_argument("--outpdf", default="genes_heatmap.pdf",
                        help="output multi-page PDF")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--yticksize", type=int, default=10)
    args = parser.parse_args()

    # List of genes to process
    genes = [
        "g10125", "g11479", "g11757", "g11895",
        "g1430", "g1431", "g1506", "g16961",
        "g18887", "g22785", "g23295", "g23303",
        "g23322", "g23338", "g23417", "g23627",
        "g2400", "g3616", "g5068", "g5622",
        "g5727", "g5917", "g5964", "g5970",
        "g7844", "g864", "g9063",
    ]

    cmap = mpl.cm.get_cmap(args.cmap).copy()
    cmap.set_bad(color="white")

    genes_per_page = 4
    n_pages = (len(genes) + genes_per_page - 1) // genes_per_page

    with PdfPages(args.outpdf) as pdf:
        for page_idx in range(n_pages):
            start = page_idx * genes_per_page
            end = min((page_idx + 1) * genes_per_page, len(genes))
            genes_this_page = genes[start:end]
            n_rows = genes_per_page

            # Slightly widen figure to preserve species labels
            fig = plt.figure(figsize=(9, 11))
            gs = gridspec.GridSpec(
                nrows=n_rows,
                ncols=2,
                width_ratios=[1.0, 5.0],
                height_ratios=[1.0] * n_rows,
                wspace=0.02,
                hspace=0.25,
                figure=fig,
            )

            # Adjust margins and spacing so labels don't overlap
            fig.subplots_adjust(
                left=0.08,
                right=0.97,
                bottom=0.06,
                top=0.98,
                hspace=0.55,
                wspace=0.05,
            )

            last_im = None
            last_axh_used = None

            for i in range(n_rows):
                row_gene_idx = start + i
                ax_tree = fig.add_subplot(gs[i, 0])
                ax_heat = fig.add_subplot(gs[i, 1])

                if row_gene_idx >= len(genes):
                    ax_tree.axis("off")
                    ax_heat.axis("off")
                    continue

                gene = genes[row_gene_idx]
                pattern = os.path.join(args.indir, f"{gene}.t1_*.matrix")
                files = sorted(glob.glob(pattern))
                if not files:
                    raise SystemExit(f"No files for gene {gene} with pattern {pattern}")

                im, axh = plot_one_gene(
                    gene_id=gene,
                    matrix_files=files,
                    ax_tree=ax_tree,
                    ax_heat=ax_heat,
                    cmap=cmap,
                    yticksize=args.yticksize,
                )
                last_im = im
                last_axh_used = axh

            # Draw colorbar on the last subplot of the page
            if last_im is not None and last_axh_used is not None:
                cax = inset_axes(
                    last_axh_used,
                    width="14%", height="22%",
                    loc="lower left",
                    bbox_to_anchor=(0.0, -0.55, 1.0, 1.0),
                    bbox_transform=last_axh_used.transAxes,
                    borderpad=0.0,
                )
                cbar = fig.colorbar(last_im, cax=cax, orientation="horizontal")
                cbar.set_ticks([0.0, 0.5, 1.0])
                cbar.ax.set_xticklabels(["0.00", "0.50", "1.00"])

            # bbox_inches="tight" ensures species labels are not truncated
            pdf.savefig(fig, dpi=args.dpi, bbox_inches="tight")
            plt.close(fig)

    print(f"Saved multi-page PDF to {args.outpdf}")


if __name__ == "__main__":
    main()
