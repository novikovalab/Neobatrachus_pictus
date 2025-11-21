#########################
# Shared SNPs along the genome
#########################

library(GenomicRanges)
library(genomation)   # transform GFF file to GRanges
library(karyoploteR)
library(RColorBrewer)
library(data.table)

## ---- user paths (edit to your own structure) ----
# Base project directory
proj_dir      <- "~/Neobatrachus_project"

# Working directory for shared SNP processing
work_dir      <- file.path(proj_dir, "shared_uniq_processing")

# Genome annotation (GFF)
gff_file      <- file.path(proj_dir, "Dovetail_genome", "augustus_hints2sources.sort.gff")

# Cytoband / chromosome size file
cytoband_file <- file.path(proj_dir, "Dovetail_genome", "biggest_chrs.txt")

# Shared SNP tables
shared_all6_file   <- file.path(proj_dir, "shared_uniq", "shared_all6_genic_regions_clean.txt")
shared_no_pic_file <- file.path(proj_dir, "shared_uniq", "shared_alb_pel_ful_sut_wil_genicregions.txt")
shared_no_pel_file <- file.path(proj_dir, "shared_uniq", "shared_alb_pic_ful_sut_wil_genicregions.txt")
shared_no_wil_file <- file.path(proj_dir, "shared_uniq", "shared_alb_pic_pel_ful_sut_genicregions.txt")
shared_no_sut_file <- file.path(proj_dir, "shared_uniq", "shared_alb_pic_pel_ful_wil_genicregions.txt")
shared_no_ful_file <- file.path(proj_dir, "shared_uniq", "shared_alb_pic_pel_sut_wil_genicregions.txt")
shared_no_alb_file <- file.path(proj_dir, "shared_uniq", "shared_pic_pel_ful_sut_wil_genicregions.txt")
## -------------------------------------------------

setwd(work_dir)

## ----- load genome annotation and SNPs (all six species) -----

gff_data   <- gffToGRanges(gff_file, filter = NULL, zero.based = FALSE, ensembl = FALSE)
all_genes  <- gff_data[gff_data$type == "gene"]

snps_all_six <- read.table(shared_all6_file, stringsAsFactors = FALSE)
colnames(snps_all_six) <- c("chr", "pos", "function", "gene")

test <- makeGRangesFromDataFrame(
  snps_all_six,
  keep.extra.columns   = TRUE,
  ignore.strand        = TRUE,
  seqnames.field       = "chr",
  start.field          = "pos",
  end.field            = "pos",
  starts.in.df.are.0based = FALSE
)

## ----- build custom genome (cytobands) -----

cytoband_N <- read.table(
  cytoband_file,
  colClasses = c("character", "numeric", "numeric"),
  sep = "\t"
)
custom.genome <- toGRanges(data.frame(
  chr   = cytoband_N$V1,
  start = cytoband_N$V2,
  end   = cytoband_N$V3
))

## =====================
## 1) Genome-wide SNP density
## =====================

# 1A. Histogram-style density (plot.type = 4)
png("segregatingSNP_density_histogram.png",
    width  = 15 * 400,
    height = 10 * 400,
    res    = 400)

pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin  <- 20

kp <- plotKaryotype(
  genome           = custom.genome,
  plot.type        = 4,
  ideogram.plotter = NULL,
  labels.plotter   = NULL,
  plot.params      = pp,
  main             = "SNP Density"
)
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt = 45)
kpPlotDensity(kp, test, window.size = 10e6, col = "#ddaacc")
dev.off()

# 1B. Ideogram density (SNPs + genes)
png("segregatingSNP_density_ideogram.png",
    width  = 15 * 400,
    height = 10 * 400,
    res    = 400)

kp <- plotKaryotype(
  genome           = custom.genome,
  plot.type        = 6,
  ideogram.plotter = NULL,
  cex              = 1.2
)
kp <- kpPlotDensity(
  kp, test,
  window.size  = 0.5e6,
  data.panel   = "ideogram",
  col          = "#3388FF",
  border       = "#3388FF",
  r0           = 0.5,
  r1           = 1.5
)
kp <- kpPlotDensity(
  kp, all_genes,
  window.size  = 0.5e6,
  data.panel   = "ideogram",
  col          = "orange",
  border       = "orange",
  r0           = 0.5,
  r1           = 0
)
legend("bottomright",
       col    = c("#3388FF", "orange"),
       legend = c("shared segregating SNPs", "genes"),
       pch    = 19,
       cex    = 1.2)
dev.off()

# PDF version
pdf("segregatingSNP_density_ideogram.pdf", width = 15)
kp <- plotKaryotype(
  genome           = custom.genome,
  plot.type        = 6,
  ideogram.plotter = NULL,
  cex              = 1.2
)
kp <- kpPlotDensity(
  kp, test,
  window.size  = 0.5e6,
  data.panel   = "ideogram",
  col          = "#3388FF",
  border       = "#3388FF",
  r0           = 0.5,
  r1           = 1.5
)
kp <- kpPlotDensity(
  kp, all_genes,
  window.size  = 0.5e6,
  data.panel   = "ideogram",
  col          = "orange",
  border       = "orange",
  r0           = 0.5,
  r1           = 0
)
legend("bottomright",
       col    = c("#3388FF", "orange"),
       legend = c("shared segregating SNPs", "genes"),
       pch    = 19,
       cex    = 1.2)
dev.off()

## =====================
## 2) 5-way combinations (missing one species each time)
## =====================

colors <- brewer.pal(6, "Dark2")

# SNPs shared by all but pictus
snps_no_pic <- read.table(shared_no_pic_file, stringsAsFactors = FALSE)
colnames(snps_no_pic) <- c("chr", "pos", "function", "gene")
GRnopic <- makeGRangesFromDataFrame(
  snps_no_pic,
  keep.extra.columns = TRUE,
  ignore.strand      = TRUE,
  start.field        = "pos",
  end.field          = "pos"
)

# SNPs shared by all but pelobatoides
snps_no_pel <- read.table(shared_no_pel_file, stringsAsFactors = FALSE)
colnames(snps_no_pel) <- c("chr", "pos", "function", "gene")
GRnopel <- makeGRangesFromDataFrame(
  snps_no_pel,
  keep.extra.columns = TRUE,
  ignore.strand      = TRUE,
  start.field        = "pos",
  end.field          = "pos"
)

# SNPs shared by all but wilsmorei
snps_no_wil <- read.table(shared_no_wil_file, stringsAsFactors = FALSE)
colnames(snps_no_wil) <- c("chr", "pos", "function", "gene")
GRnowil <- makeGRangesFromDataFrame(
  snps_no_wil,
  keep.extra.columns = TRUE,
  ignore.strand      = TRUE,
  start.field        = "pos",
  end.field          = "pos"
)

# SNPs shared by all but sutor
snps_no_sut <- read.table(shared_no_sut_file, stringsAsFactors = FALSE)
colnames(snps_no_sut) <- c("chr", "pos", "function", "gene")
GRnosut <- makeGRangesFromDataFrame(
  snps_no_sut,
  keep.extra.columns = TRUE,
  ignore.strand      = TRUE,
  start.field        = "pos",
  end.field          = "pos"
)

# SNPs shared by all but fulvus
snps_no_ful <- read.table(shared_no_ful_file, stringsAsFactors = FALSE)
colnames(snps_no_ful) <- c("chr", "pos", "function", "gene")
GRnoful <- makeGRangesFromDataFrame(
  snps_no_ful,
  keep.extra.columns = TRUE,
  ignore.strand      = TRUE,
  start.field        = "pos",
  end.field          = "pos"
)

# SNPs shared by all but albipes
snps_no_alb <- read.table(shared_no_alb_file, stringsAsFactors = FALSE)
colnames(snps_no_alb) <- c("chr", "pos", "function", "gene")
GRnoalb <- makeGRangesFromDataFrame(
  snps_no_alb,
  keep.extra.columns = TRUE,
  ignore.strand      = TRUE,
  start.field        = "pos",
  end.field          = "pos"
)

# Quick ideogram overlay (not saved as file here)
kp <- plotKaryotype(
  genome           = custom.genome,
  plot.type        = 6,
  ideogram.plotter = NULL,
  cex              = 1.2
)
kp <- kpPlotDensity(kp, GRnopic, window.size = 0.5e6,
                    data.panel = "ideogram", col = colors[1], border = colors[1],
                    r0 = 0.5, r1 = 1.5)
kp <- kpPlotDensity(kp, GRnopel, window.size = 0.5e6,
                    data.panel = "ideogram", col = colors[2], border = colors[2],
                    r0 = 0.5, r1 = 1.5)
kp <- kpPlotDensity(kp, GRnowil, window.size = 0.5e6,
                    data.panel = "ideogram", col = colors[3], border = colors[3],
                    r0 = 0.5, r1 = 1.5)
kp <- kpPlotDensity(kp, GRnosut, window.size = 0.5e6,
                    data.panel = "ideogram", col = colors[4], border = colors[4],
                    r0 = 0.5, r1 = 1.5)
kp <- kpPlotDensity(kp, GRnoful, window.size = 0.5e6,
                    data.panel = "ideogram", col = colors[5], border = colors[5],
                    r0 = 0.5, r1 = 1.5)
kp <- kpPlotDensity(kp, GRnoalb, window.size = 0.5e6,
                    data.panel = "ideogram", col = colors[6], border = colors[6],
                    r0 = 0.5, r1 = 1.5)

kp <- kpPlotDensity(kp, test, window.size = 0.5e6,
                    data.panel = "ideogram", col = "darkgrey", border = "darkgrey",
                    r0 = 0.5, r1 = 0)

## Main 5-way density panel
png("segregatingSNP_density_ideogram_all5way.png",
    width  = 15 * 400,
    height = 10 * 400,
    res    = 400)

pp <- getDefaultPlotParams(plot.type = 2)
pp$data1height <- 50

kp <- plotKaryotype(
  genome           = custom.genome,
  plot.type        = 1,
  ideogram.plotter = NULL
)

# Background: all shared SNPs
kpPlotRegions(kp, data = test, col = "#CCCCCC44", border = "#CCCCCC44")

# Stacked density tracks for each 5-way combination
kpPlotDensity(kp, GRnopic, data.panel = 1, col = colors[1],
              window.size = 0.5e5, r0 = 0,      r1 = 0.167)
kpPlotDensity(kp, GRnopel, data.panel = 1, col = colors[2],
              window.size = 0.5e5, r0 = 0.167,  r1 = 0.33)
kpPlotDensity(kp, GRnowil, data.panel = 1, col = colors[3],
              window.size = 0.5e5, r0 = 0.33,   r1 = 0.5)
kpPlotDensity(kp, GRnosut, data.panel = 1, col = colors[4],
              window.size = 0.5e5, r0 = 0.5,    r1 = 0.667)
kpPlotDensity(kp, GRnoful, data.panel = 1, col = colors[5],
              window.size = 0.5e5, r0 = 0.667,  r1 = 0.833)
kpPlotDensity(kp, GRnoalb, data.panel = 1, col = colors[6],
              window.size = 0.5e5, r0 = 0.833,  r1 = 1)

leg <- c(
  expression(paste("no  ", italic(N.), " ", italic(pictus))),
  expression(paste("no  ", italic(N.), " ", italic(pelobatoides))),
  expression(paste("no  ", italic(N.), " ", italic(wilsmorei))),
  expression(paste("no  ", italic(N.), " ", italic(sutor))),
  expression(paste("no  ", italic(N.), " ", italic(fulvus))),
  expression(paste("no  ", italic(N.), " ", italic(albipes)))
)

legend("right", col = colors, legend = leg, pch = 19, cex = 1.2)
dev.off()

## =====================
## 3) Zoom in on a density-rich region (example)
## =====================

zoom.region <- toGRanges(data.frame("Scaffold_17238", 160e6, 230e6))

kp <- plotKaryotype(
  plot.type        = 2,
  genome           = custom.genome,
  chromosomes      = "Scaffold_17238",
  ideogram.plotter = NULL,
  zoom             = zoom.region
)
kp <- kpPlotDensity(kp, data = test, window.size = 0.5e5)
kpAddBaseNumbers(kp)
kpPlotDensity(kp, data = all_genes, data.panel = 2, window.size = 0.5e5)
kpAxis(kp,
       ymax = kp$latest.plot$computed.values$max.density,
       cex  = 0.8)
kpAbline(kp,
         h    = 0.3 * kp$latest.plot$computed.values$max.density,
         lty  = 2,
         ymax = kp$latest.plot$computed.values$max.density)

outlier_windows_index <- which(
  kp$latest.plot$computed.values$density >
    0.3 * kp$latest.plot$computed.values$max.density
)
outlier_windows <- kp$latest.plot$computed.values$windows[outlier_windows_index]
kpPlotRegions(kp, data = outlier_windows, col = NA, data.panel = 2)

kp <- plotKaryotype(
  plot.type        = 2,
  genome           = custom.genome,
  chromosomes      = "Scaffold_17238",
  ideogram.plotter = NULL,
  zoom             = zoom.region
)
kp <- kpPlotDensity(kp, test,
                    window.size = 0.5e6,
                    col         = "#3388FF",
                    border      = "#3388FF",
                    data.panel  = 1)
kp <- kpPlotDensity(kp, all_genes,
                    window.size = 0.5e6,
                    data.panel  = "ideogram",
                    col         = "orange",
                    border      = "orange",
                    r0          = 0.5,
                    r1          = 0)
legend("bottomright",
       col    = c("#3388FF", "orange"),
       legend = c("shared segregating SNPs", "genes"),
       pch    = 19,
       cex    = 1.2)

# Extract SNPs in outlier windows for downstream GO analysis
outlier_snps <- subsetByOverlaps(test, outlier_windows)
# 'outlier_snps' can now be used for GO enrichment of genes containing these SNPs.
