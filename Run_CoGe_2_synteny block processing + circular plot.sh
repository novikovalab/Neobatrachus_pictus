## Step 1: Build synteny blocks from CoGe starts/ends

library(dplyr)

# Read start and end coordinates extracted in Bash
mystarts <- read.delim("x.starts", header = FALSE)
myends   <- read.delim("x.ends",   header = FALSE)

# Extra line that should be appended to 'myends'
# (equivalent to manually adding the last line from the CoGe file)
extraline <- c(
  1.4452, 0.0018,
  "a52485_NW_016683643.1",
  "NW_016683643.1||64988||68412||cds-NP_001072450.1||1||CDS||2471794141||80653||87.67",
  64988, 68412,
  "b59218_Scaffold_3532",
  "Scaffold_3532||4231944||4234578||g1407||-1||CDS||3198286918||6938||87.67",
  4234578, 4231944,
  1.000000e-250, 216
)

# Append the extra line so starts and ends have matching lengths
myends <- rbind(myends, extraline)

# Keep only relevant columns: chr and start/end in each species
mystarts <- mystarts[, c(3, 5, 7, 9)]
myends   <- myends[,   c(3, 5, 7, 9)]

# Combine starts and ends
mystart_ends        <- cbind(mystarts, myends)
mystart_ends_blocks <- mystart_ends[, c(1, 2, 6, 3, 4, 8)]
colnames(mystart_ends_blocks) <- c(
  "xenopus_chr", "xenopus_start", "xenopus_end",
  "neobatrachus_chr", "neobatrachus_start", "neobatrachus_end"
)

# Clean species-specific prefixes from chromosome names
mystart_ends_blocks$xenopus_chr      <- gsub("a52485_", "", mystart_ends_blocks$xenopus_chr)
mystart_ends_blocks$neobatrachus_chr <- gsub("b59218_", "", mystart_ends_blocks$neobatrachus_chr)

# Write unfiltered synteny blocks (all scaffolds)
write.table(
  mystart_ends_blocks,
  file      = "synteny.unflipped.txt",
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE,
  sep       = "\t"
)

## Optional: filter for largest Xenopus scaffolds (done manually beforehand)
## The filtered file is assumed to be saved as:
##   ~/CoGe_analysis/synteny.unflipped.biggestScX.filt.txt

data <- read.delim(
  "~/CoGe_analysis/synteny.unflipped.biggestScX.filt.txt",
  colClasses = c("character", "numeric", "numeric",
                 "character", "numeric", "numeric"),
  header     = TRUE
)

## Step 2: Prepare cytobands / chromosome size data

# Cytoband-style file for Neobatrachus:
# columns: chr, start (0-based), end
cytoband_N <- read.table(
  "~/Dovetail_genome/biggest_chrs.txt",
  colClasses = c("character", "numeric", "numeric"),
  sep        = "\t"
)

# Cytoband file for Xenopus:
# columns: chr, start, something, end (we only use selected columns)
cytoband_X <- read.delim(
  "~/CoGe_analysis/cytobands_X.txt",
  colClasses = c("character", "numeric", "numeric", "numeric"),
  header     = FALSE
)

# Combine cytobands of both species.
# Set all start coordinates to 0 to standardize.
library(data.table)
cytoband <- rbindlist(list(cytoband_N, cytoband_X[, c(1, 2, 4)]))
cytoband$V2 <- rep(0, nrow(cytoband))

## Step 3: Circular synteny plot with circlize

library(circlize)
library(RColorBrewer)

png(
  "synteny_plot_fullcolor.png",
  width  = 15 * 400,
  height = 10 * 400,
  res    = 400
)

# Clear any previous circos plot configuration
circos.clear()

# Basic circos plot parameters
circos.par(
  "start.degree" = 177,
  "gap.degree"   = rep(4, 23),
  canvas.ylim    = c(-1, 1),
  canvas.xlim    = c(-1, 1),
  clock.wise     = TRUE,
  track.margin   = c(0, 0),
  track.height   = 0.1
)

# Ensure chromosomes are plotted in the order of appearance in cytoband
cytoband[[1]] <- factor(cytoband[[1]], levels = cytoband$V1)

# Initialize genomic layout without default axes/labels
circos.genomicInitialize(
  cytoband,
  plotType              = NULL,
  tickLabelsStartFromZero = TRUE,
  axis.labels.cex       = 0.3 * par("cex"),
  labels.cex            = 0.9 * par("cex"),
  track.height          = 0.1
)

# Outer track: chromosome segments and labels
circos.track(
  ylim       = c(0, 0.02),
  bg.col     = c(rep("#7fbc41", 13), rep("#276419", 10)),  # colors for each species
  bg.border  = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    chr  <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    # Add chromosome label
    circos.text(
      mean(xlim), mean(ylim) + 0.06, chr,
      cex       = 1,
      facing    = "inside",
      niceFacing = TRUE
    )
    # Add axis (Mb scale)
    circos.axis(
      major.tick    = TRUE,
      minor.ticks   = 0,
      labels        = c("", "10MB", "20MB", "30MB"),
      labels.facing = "inside",
      labels.cex    = 0.8 * par("cex")
    )
  }
)

# Colors for synteny links
mycolors <- colorRampPalette(brewer.pal(8, "Set1"), alpha = TRUE)(13)

# Get list of Neobatrachus scaffolds / chromosomes
Pictus_scaffolds <- unique(data$neobatrachus_chr)

# Draw synteny links, one color per Neobatrachus scaffold
for (i in seq_len(13)) {
  data_scaffold <- dplyr::filter(data, neobatrachus_chr == Pictus_scaffolds[i])
  circos.genomicLink(
    data_scaffold[, 1:3],
    data_scaffold[, 4:6],
    col    = mycolors[i],
    border = mycolors[i]
  )
}

# Add legend for species
legend(
  x      = 0.8,
  y      = 1,
  legend = c("N. pictus", "X. tropicalis"),
  col    = c("#7fbc41", "#276419"),
  lty    = 1,
  lwd    = 7
)

dev.off()

## done
