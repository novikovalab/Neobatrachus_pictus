########################################
## 1. Load required packages
########################################
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(patchwork)  # or use cowplot

########################################
## 2. Define and load the tree
########################################

# Example Newick string
tree_text <- "((((((Nsutor,Nwilsmorei),Nalbipes),Nfulvus),Npictus),Npelobatoides),Ldumerilii);"
tree <- read.tree(text = tree_text)

# Generate the left tree (roundrect layout, ignore branch lengths)
p_left <- ggtree(tree, layout = "roundrect", branch.length = "none") +
  geom_tiplab(align = TRUE, offset = 0.5, linesize = NA)
p_left

# Generate the right tree
p_right <- ggtree(tree, layout = "roundrect", branch.length = "none") +
  geom_tiplab(align = TRUE, offset = 0.5, linesize = NA)

# Extract data as data frames (avoids “Invalid edge matrix” warnings)
d_left <- as.data.frame(p_left$data)
d_right <- as.data.frame(p_right$data)

# Ensure tip labels are character type
d_left$label <- as.character(d_left$label)
d_right$label <- as.character(d_right$label)

# Mirror the right tree: reverse x-axis and shift to the right (offset = 5)
d_right$x <- max(d_right$x) - d_right$x + max(d_left$x) + 5

########################################
## 3. Read CSV data and compute link frequencies
########################################

# Modify the CSV path according to your environment
df_raw <- read.csv("Tetra_assignments.csv", stringsAsFactors = FALSE)

# Function to clean and standardize taxon names from the CSV
clean_name <- function(x) {
  x <- tolower(x)
  x <- gsub("^n\\.?\\s*", "", x)  # remove leading n. or n
  x <- gsub("\\.", "", x)         # remove dots
  x <- trimws(x)
  dict <- c(
    "sutor"        = "Nsutor",
    "wilsmorei"    = "Nwilsmorei",
    "albipes"      = "Nalbipes",
    "fulvus"       = "Nfulvus",
    "pictus"       = "Npictus",
    "pelobatoides" = "Npelobatoides",
    "dumerilii"    = "Ldumerilii"
  )
  ifelse(x %in% names(dict), dict[x], x)
}

# Apply cleaning to partner columns
df <- df_raw %>%
  mutate(
    nsu1_partner = clean_name(nsu1_partner),
    nsu2_partner = clean_name(nsu2_partner)
  )

# Count occurrences and compute relative frequencies
df_freq <- df %>%
  group_by(nsu1_partner, nsu2_partner) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(freq = count / sum(count))

# Extract tip coordinates for left and right trees
tips_left  <- d_left  %>% filter(isTip) %>% select(label, x, y)
tips_right <- d_right %>% filter(isTip) %>% select(label, x, y)

# Merge link frequencies with tree coordinates
df_freq <- df_freq %>%
  left_join(tips_left,  by = c("nsu1_partner" = "label")) %>%
  rename(x_left = x, y_left = y) %>%
  left_join(tips_right, by = c("nsu2_partner" = "label")) %>%
  rename(x_right = x, y_right = y)

# Print merged table for checking
print(df_freq)

########################################
## 4. Define species color mapping
########################################
species_colors <- c(
  "Nsutor"        = "#66c2a5",
  "Nwilsmorei"    = "#ffff99",
  "Nalbipes"      = "#fc8d62",
  "Nfulvus"       = "#8da0cb",
  "Npictus"       = "#e78ac3",
  "Npelobatoides" = "grey70",
  "Ldumerilii"    = "grey40"
)

########################################
## 5. Build combined plot
########################################

d_left_tip  <- d_left  %>% filter(isTip)
d_right_tip <- d_right %>% filter(isTip)

combined_plot <- ggplot() +
  # --- Base layer: connecting lines between trees ---
  geom_segment(
    data = df_freq,
    aes(x = x_left, y = y_left, xend = x_right, yend = y_right, size = freq),
    color = "grey40", alpha = 1, show.legend = TRUE
  ) +
  # --- Left tree ---
  geom_tree(data = d_left, layout = "roundrect") +
  geom_label(
    data = d_left_tip,
    aes(x = x, y = y, label = label, fill = label),
    label.size = 0,
    fontface = "italic",
    color = "black",
    alpha = 1,
    position = position_nudge(x = -0.5)
  ) +
  # --- Right tree ---
  geom_tree(data = d_right, layout = "roundrect") +
  geom_label(
    data = d_right_tip,
    aes(x = x, y = y, label = label, fill = label),
    label.size = 0,
    fontface = "italic",
    color = "black",
    alpha = 1,
    position = position_nudge(x = 0.5)
  ) +
  # --- Apply species colors ---
  scale_fill_manual(values = species_colors, guide = "none") +
  # --- Size legend for link frequency ---
  scale_size_continuous(
    name = "Relative occurrence",
    limits = c(0,1),
    breaks = c(0.1, 0.2, 0.3, 0.5, 1.0),
    range = c(0.5, 3),
    guide = guide_legend(override.aes = list(linetype = 1, color = "grey40"))
  ) +
  ggtitle(expression(italic("N. kunapalari"))) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "horizontal",
    plot.title = element_text(size = 18, hjust = 0.5)
  )

print(combined_plot)

########################################
## 6. Save output figures
########################################
ggsave("N_sudellae_FINAL.png", combined_plot, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("N_sudellae_FINAL.svg", combined_plot, width = 10, height = 6, dpi = 300, bg = "white")
