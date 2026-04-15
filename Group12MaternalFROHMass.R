# Alexandra Kelly
# Final Maternal FROH vs Pup Mass Gain
# April 10, 2026

setwd("/Users/alexandrakelly/Desktop/EEB331")
getwd()

library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)

# Recalculate FROH

# Constant
genome_size <- 2319465413

# Read in ROH output from bcftools
lines <- readLines("Marmots_maf3_geno90_mind20_ROH.txt")
lines <- lines[!grepl("^#", lines)]
rg_lines <- lines[grepl("^RG", lines)]

# Extract ROH segments
roh_segments <- read.table(
  text = rg_lines,
  header = FALSE,
  stringsAsFactors = FALSE
)

# Assign column names
colnames(roh_segments) <- c(
  "type",
  "sample",
  "chr",
  "start",
  "end",
  "length_bp",
  "n_markers",
  "quality"
)

# Summarize ROH by individual
roh_summary <- roh_segments %>%
  mutate(length_bp = as.numeric(end) - as.numeric(start)) %>%
  group_by(sample) %>%
  summarise(
    NSEG = n(),
    total_bp = sum(length_bp, na.rm = TRUE),
    mean_bp = mean(length_bp, na.rm = TRUE),
    mean_length = total_bp / NSEG,
    .groups = "drop"
  )

# Calculate FROH
roh_summary <- roh_summary %>%
  mutate(FROH = total_bp / genome_size)

# Extract numeric uid from sample names
roh_summary <- roh_summary %>%
  mutate(uid = str_extract(sample, "[0-9]+_[0-9]+"))

# Keep only uid and FROH
roh_clean <- roh_summary %>%
  select(uid, FROH)

# Read in metadata

ped_data <- read.csv("marmot_Parentage_Informed_Pedigree.csv")
LH <- read.csv("RADseq_marmot_metadata.csv")
meta_pedigree <- read_excel("RADseq_pedigree_metadata.xlsx")
mass_metadata <- read_excel("marmot_mass_metadata.xlsx")

# Prepare mass data
meta_mass <- mass_metadata %>%
  select(uid, year, massjun, massaug, nb_mass) %>%
  mutate(mass_gain = massaug - massjun)

# Pair pups and mothers
pup_mom <- ped_data %>%
  select(id, dam) %>%
  filter(!is.na(dam) & dam != "") %>%
  rename(
    pup_id = id,
    mom_id = dam
  ) %>%
  mutate(
    mom_uid = str_extract(mom_id, "[0-9]+_[0-9]+"),
    pup_uid = str_extract(pup_id, "[0-9]+_[0-9]+")
  ) %>%
  select(pup_uid, mom_uid)

# Merge maternal FROH onto each pup
combined_froh_id <- merge(
  pup_mom,
  roh_clean,
  by.x = "mom_uid",
  by.y = "uid",
  all.x = TRUE
)

# Rename maternal FROH column
colnames(combined_froh_id)[colnames(combined_froh_id) == "FROH"] <- "mom_FROH"

# Reorder columns
combined_froh_id <- combined_froh_id[, c("pup_uid", "mom_uid", "mom_FROH")]

# Remove rows missing pup or mom uid
combined_froh_id <- combined_froh_id %>%
  filter(!is.na(pup_uid), !is.na(mom_uid))

# Merge with mass data and apply stricter filter
analysis_mass_froh <- combined_froh_id %>%
  inner_join(meta_mass, by = c("pup_uid" = "uid")) %>%
  filter(
    !is.na(mom_FROH),
    !is.na(massjun),
    !is.na(massaug),
    !is.na(mass_gain),
    !is.na(nb_mass),
    nb_mass >= 2
  )

# Check sample size
nrow(analysis_mass_froh)
# 48 individuals

# Distribution checks
shapiro.test(analysis_mass_froh$mom_FROH)
# W = 0.83714, p-value = 1.008e-05
shapiro.test(analysis_mass_froh$mass_gain)
# W = 0.97789, p-value = 0.4943
shapiro.test(analysis_mass_froh$massjun)
# W = 0.96296, p-value = 0.1329
shapiro.test(analysis_mass_froh$massaug)
# W = 0.97025, p-value = 0.259

# Histograms
par(mfrow = c(1, 4))
hist(
  analysis_mass_froh$mom_FROH,
  main = "Distribution of Mom FROH",
  xlab = "Mom FROH",
  col = "lightblue"
)
hist(
  analysis_mass_froh$mass_gain,
  main = "Distribution of Mass Gain",
  xlab = "Mass Gain (g)",
  col = "lightgreen"
)
hist(
  analysis_mass_froh$massjun,
  main = "Distribution of June Mass",
  xlab = "June Mass (g)",
  col = "tomato"
)
hist(
  analysis_mass_froh$massaug,
  main = "Distribution of August Mass",
  xlab = "August Mass (g)",
  col = "khaki"
)

# QQ plots
par(mfrow = c(1, 4))
qqnorm(analysis_mass_froh$mom_FROH, main = "QQ Plot: Mom FROH")
qqline(analysis_mass_froh$mom_FROH, col = "red")

qqnorm(analysis_mass_froh$mass_gain, main = "QQ Plot: Mass Gain")
qqline(analysis_mass_froh$mass_gain, col = "red")

qqnorm(analysis_mass_froh$massjun, main = "QQ Plot: June Mass")
qqline(analysis_mass_froh$massjun, col = "red")

qqnorm(analysis_mass_froh$massaug, main = "QQ Plot: August Mass")
qqline(analysis_mass_froh$massaug, col = "red")

# Reset plotting layout
par(mfrow = c(1, 1))

# Correlation tests
# Use Spearman consistently for all three traits
# exact = FALSE avoids tie warnings

sp_gain <- cor.test(
  analysis_mass_froh$mom_FROH,
  analysis_mass_froh$mass_gain,
  method = "spearman",
  exact = FALSE
)

sp_aug <- cor.test(
  analysis_mass_froh$mom_FROH,
  analysis_mass_froh$massaug,
  method = "spearman",
  exact = FALSE
)

sp_jun <- cor.test(
  analysis_mass_froh$mom_FROH,
  analysis_mass_froh$massjun,
  method = "spearman",
  exact = FALSE
)

sp_gain
# S = 23702, p-value = 0.04838
# rho: -0.2864757 
sp_aug
# S = 22614, p-value = 0.12
# rho: -0.227441 
sp_jun
# S = 15675, p-value = 0.3114
# rho: 0.1492174 

# Final figures for submission

# Figure 1: Maternal FROH vs Pup Mass Gain
ggplot(analysis_mass_froh, aes(x = mom_FROH, y = mass_gain)) +
  geom_point(alpha = 0.6, color = "#2c7bb6", size = 2.5) +
  geom_smooth(method = "loess", se = TRUE, color = "#d7191c", linewidth = 0.8) +
  annotate(
    "text",
    x = max(analysis_mass_froh$mom_FROH, na.rm = TRUE),
    y = max(analysis_mass_froh$mass_gain, na.rm = TRUE),
    label = paste0(
      "Spearman rho = ", round(sp_gain$estimate, 3),
      "\np = ", signif(sp_gain$p.value, 3)
    ),
    hjust = 1,
    vjust = 1,
    size = 4,
    fontface = "italic"
  ) +
  labs(
    title = "Maternal FROH vs. Pup Seasonal Mass Gain",
    subtitle = "Filtered to pups with complete June/August mass data and nb_mass >= 2",
    x = expression(Maternal ~ F[ROH]),
    y = "Pup Mass Gain (g)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Figure 2: Maternal FROH vs Pup August Mass
ggplot(analysis_mass_froh, aes(x = mom_FROH, y = massaug)) +
  geom_point(alpha = 0.6, color = "#2c7bb6", size = 2.5) +
  geom_smooth(method = "loess", se = TRUE, color = "#d7191c", linewidth = 0.8) +
  annotate(
    "text",
    x = max(analysis_mass_froh$mom_FROH, na.rm = TRUE),
    y = max(analysis_mass_froh$massaug, na.rm = TRUE),
    label = paste0(
      "Spearman rho = ", round(sp_aug$estimate, 3),
      "\np = ", signif(sp_aug$p.value, 3)
    ),
    hjust = 1,
    vjust = 1,
    size = 4,
    fontface = "italic"
  ) +
  labs(
    title = "Maternal FROH vs. Pup August Mass",
    subtitle = "Filtered to pups with complete June/August mass data and nb_mass >= 2",
    x = expression(Maternal ~ F[ROH]),
    y = "Pup August Mass (g)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Figure 3: Maternal FROH vs Pup June Mass
ggplot(analysis_mass_froh, aes(x = mom_FROH, y = massjun)) +
  geom_point(alpha = 0.6, color = "#2c7bb6", size = 2.5) +
  geom_smooth(method = "loess", se = TRUE, color = "#d7191c", linewidth = 0.8) +
  annotate(
    "text",
    x = max(analysis_mass_froh$mom_FROH, na.rm = TRUE),
    y = max(analysis_mass_froh$massjun, na.rm = TRUE),
    label = paste0(
      "Spearman rho = ", round(sp_jun$estimate, 3),
      "\np = ", signif(sp_jun$p.value, 3)
    ),
    hjust = 1,
    vjust = 1,
    size = 4,
    fontface = "italic"
  ) +
  labs(
    title = "Maternal FROH vs. Pup June Mass",
    subtitle = "Filtered to pups with complete June/August mass data and nb_mass >= 2",
    x = expression(Maternal ~ F[ROH]),
    y = "Pup June Mass (g)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
