##################
# PCA admixture tree barplot
##################

### Step 1: PCA ###
# Set working directory
setwd("<PATH_TO_EV_ANALYSIS>")

# Input data
# Only first scaffold
data_sc1 <- read.table("<PATH_TO_PARSED_VCF>", header = TRUE, sep = ",")

# 4fds sites
data_all <- read.table("<PATH_TO_4FDS_FILE>", header = TRUE, sep = ",")

# Remove rows with NA
data <- na.omit(data_all)

# Remove problematic sample Nwilsmorei_15834 (identified as erroneous earlier)
data <- data[, -90]

# Transpose genotype matrix: rows = samples, columns = SNPs
transformed <- t(data[, 3:dim(data)[2]])

# Prepare color table for PCA plot
library(RColorBrewer)
colors <- brewer.pal(9, "Paired")

species <- c(
  "Nalbipes", "Naquilonius", "Nfulvus", "Nkunapalari",
  "Npelobatoides", "Npictus", "Nsudellae", "Nsutor", "Nwilsmorei"
)

mycol <- rownames(transformed)
PCA_legend <- c(19, 17, 19, 17, 19, 19, 17, 19, 19)

# Assign colors for each species
for (i in 1:9) {
  S <- species[i]
  mycol[grep(S, rownames(transformed))] <- colors[i]
}

colortable <- cbind(rownames(transformed), mycol)

# Set point shapes (pch): triangles for tetraploids
mydot <- rep(19, times = dim(transformed)[1])
for (i in 1:9) {
  S <- species[i]
  if (S == "Nsudellae" | S == "Nkunapalari" | S == "Naquilonius") {
    mydot[grep(S, rownames(transformed))] <- 17
  }
}

################################
# Run PCA
################################
library(factoextra)

data.pca <- prcomp(transformed)
my <- data.pca$x                      # PCA coordinates of samples
data.ind <- get_pca_ind(data.pca)
eig.val <- get_eigenvalue(data.pca)

# PDF with explained variance + PCA plot
pdf('PCA_4fdsites_allsamples_labels.pdf', width = 14)
par(mfrow = c(1, 2))

barplot(eig.val$variance.percent[1:10],
        names.arg = paste('PC', 1:10, sep = ""),
        ylab = 'Percentage of explained variance')

plot(my[, 1], my[, 2], col = adjustcolor(mycol, alpha.f = 0.5),
     xlab = 'PCA1', ylab = 'PCA2', pch = mydot, cex = 1.5)

legend('topright', legend = species, text.col = colors,
       bty = "n", pch = PCA_legend)

dev.off()

################################
# PCA with sample labels
################################
library(basicPlotteR)

png("PCA_4fdsites_allsamples_labels.png",
    width = 9 * 300,
    height = 6 * 300,
    res = 300)

plot(my[, 1], my[, 2], col = mycol,
     xlab = 'PCA1', ylab = 'PCA2',
     pch = mydot, cex = 1.5)

addTextLabels(my[, 1], my[, 2], rownames(my),
              cex.label = 0.7, col.label = mycol,
              lty = 3, col.line = 'dark grey')

legend('topright', legend = species, text.col = colors,
       bty = "n", pch = PCA_legend, cex = 1)

dev.off()



#############################
### Step 2: Tree assembly ###
#############################
library(ape)

# Distance matrix using Euclidean distance
mydist <- dist(transformed)

# Visualize distance matrix
library(factoextra)
fviz_dist(mydist)

pdf('mydist.pdf', width = 14)
fviz_dist(mydist)
dev.off()

# Build and visualize hierarchical and NJ trees
pdf('tree.pdf', width = 14)
par(mfrow = c(1, 2))

plot(hclust(mydist))                        # hierarchical clustering tree

njtree <- nj(mydist)                        # neighbor-joining tree
plot(root(njtree, outgroup = "Npelobatoides_135619"))
dev.off()

# Only NJ tree as PNG
png("NJtree.png",
    width = 21 * 400,
    height = 13 * 400,
    res = 400)

njtree <- nj(mydist)
plot(root(njtree, outgroup = "Npelobatoides_135619"))

dev.off()



################################
### Step 3: Admixture analysis
################################
setwd("<PATH_TO_ADMIXTURE_FOLDER>")
outputfolder = "<PATH_TO_ADMIX_OUTPUT>"

# Choose optimal K using CV error
Kval <- read.table("CV_k.txt")

plot(Kval$V1, Kval$V5, xlab = 'K', ylab = 'Cross-validation error')
# The smallest CV error determines optimal K

K = 8   # Insert chosen K

tbl <- read.table(paste(outputfolder,
                        "Neobatrachus.admixture20maf2.", K, ".Q",
                        sep = ''))

tbl <- tbl[-88, ]  # Remove erroneous sample



################################
# Plot admixture barplots
################################
pdf(file = paste("admixture20maf2", K, "originalorder.pdf", sep='.'), width = 16)

par(mar = c(10, 4, 4, 2))
colors <- brewer.pal(K, "RdYlBu")

barplot(t(as.matrix(tbl)), col = colors,
        ylab = "Ancestry", border = NA,
        cex.names = 0.7, las = 2, main = K)

dev.off()



##############################################
### Step 4: Combine admixture with NJ tree ###
##############################################

# Prepare NJ tree (rooted)
mydist <- dist(transformed)
hierarchclust <- hclust(mydist)
njtree <- nj(mydist)
njtree_root <- root(njtree, outgroup = "Npelobatoides_135619")

plot(njtree_root, align.tip.label = TRUE)

# Extract sample order from tree
library(dendextend)

tree_order <- order.hclust(hierarchclust)

is_tip <- njtree_root$edge[, 2] <= length(njtree_root$tip.label)
ordered_tips <- njtree_root$edge[is_tip, 2]

reorder_tbl <- tbl[ordered_tips, ]


##############################################
# Combine NJ tree + Admixture barplot
##############################################
png(paste("hierarch+admix", K, "njtree.png", sep = '.'),
    width = 21 * 200,
    height = 13 * 200,
    res = 200)

par(mfrow = c(1, 2), mar = c(2.5, 0, 1, 4))

plot(root(njtree, outgroup = "Npelobatoides_135619"), align.tip.label = TRUE)

colors <- brewer.pal(K, "RdYlBu")

barplot(t(as.matrix(reorder_tbl)), col = colors,
        border = NA,
        names.arg = vector(mode = "character", length = 87),
        las = 2, horiz = TRUE)

dev.off()
