# Group 12: Alexandra Kelly and Jamil Fayad
# Maternal FROH vs. pup mass

# Setup

setwd("/Users/alexandrakelly/Desktop/EEB331")

library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)

# Total genome size in base pairs, used to calculate FROH
genome_size <- 2319465413

# Read and summarize ROH data

# Read bcftools ROH output
lines <- readLines("Marmots_maf3_geno90_mind20_ROH.txt")

# Remove comment lines and keep only RG entries (ROH segments)
rg_lines <- lines[!grepl("^#", lines) & grepl("^RG", lines)]

# Convert ROH segment lines into a data frame
roh_segments <- read.table(
  text = rg_lines,
  header = FALSE,
  stringsAsFactors = FALSE
)

# Assign column names based on file structure
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

# Summarize ROH per individual and calculate FROH
roh_clean <- roh_segments %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    length_bp = end - start,
    uid = str_extract(sample, "[0-9]+_[0-9]+")
  ) %>%
  group_by(sample, uid) %>%
  summarise(
    NSEG = n(),
    total_bp = sum(length_bp, na.rm = TRUE),
    mean_bp = mean(length_bp, na.rm = TRUE),
    mean_length = total_bp / NSEG,
    FROH = total_bp / genome_size,
    .groups = "drop"
  ) %>%
  select(uid, FROH)

# Read metadata

# Parentage-informed pedigree
ped_data <- read.csv("marmot_Parentage_Informed_Pedigree.csv")

# Mass metadata
mass_metadata <- read_excel("marmot_mass_metadata.xlsx")

# Keep only variables needed for analysis
meta_mass <- mass_metadata %>%
  select(uid, year, massjun, massaug, nb_mass) %>%
  mutate(mass_gain = massaug - massjun)

# Match pups to mothers

# Extract pup and mother IDs from pedigree file
pup_mom <- ped_data %>%
  select(id, dam) %>%
  filter(!is.na(dam), dam != "") %>%
  transmute(
    pup_uid = str_extract(id, "[0-9]+_[0-9]+"),
    mom_uid = str_extract(dam, "[0-9]+_[0-9]+")
  ) %>%
  filter(!is.na(pup_uid), !is.na(mom_uid))

# Attach each mother's FROH value to her pup
combined_froh_id <- pup_mom %>%
  left_join(roh_clean, by = c("mom_uid" = "uid")) %>%
  rename(mom_FROH = FROH)

# Build final analysis dataset

analysis_mass_froh <- combined_froh_id %>%
  inner_join(meta_mass, by = c("pup_uid" = "uid")) %>%
  filter(
    !is.na(mom_FROH),
    !is.na(massjun),
    !is.na(massaug),
    !is.na(mass_gain),
    !is.na(nb_mass)
  )

# Check final sample size
nrow(analysis_mass_froh)

# Distribution checks

# Shapiro-Wilk tests
shapiro.test(analysis_mass_froh$mom_FROH)
shapiro.test(analysis_mass_froh$mass_gain)
shapiro.test(analysis_mass_froh$massjun)
shapiro.test(analysis_mass_froh$massaug)

# Histograms
par(mfrow = c(1, 4))
hist(analysis_mass_froh$mom_FROH,
     main = "Distribution of Mom FROH",
     xlab = "Mom FROH",
     col = "lightblue")

hist(analysis_mass_froh$mass_gain,
     main = "Distribution of Mass Gain",
     xlab = "Mass Gain (g)",
     col = "lightgreen")

hist(analysis_mass_froh$massjun,
     main = "Distribution of June Mass",
     xlab = "June Mass (g)",
     col = "tomato")

hist(analysis_mass_froh$massaug,
     main = "Distribution of August Mass",
     xlab = "August Mass (g)",
     col = "khaki")

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

# Reset plotting window
par(mfrow = c(1, 1))

# Correlation tests

# Spearman is used consistently here because that became the final approach
# exact = FALSE avoids warnings when ties are present

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

# Print results
sp_gain
sp_aug
sp_jun

# Final figures

# Figure 1: Maternal FROH vs pup seasonal mass gain
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
    subtitle = "Filtered to pups with complete June/August mass data",
    x = expression(Maternal ~ F[ROH]),
    y = "Pup Mass Gain (g)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Figure 2: Maternal FROH vs pup August mass
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
    subtitle = "Filtered to pups with complete June/August mass data",
    x = expression(Maternal ~ F[ROH]),
    y = "Pup August Mass (g)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Figure 3: Maternal FROH vs pup June mass
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
    subtitle = "Filtered to pups with complete June/August mass data",
    x = expression(Maternal ~ F[ROH]),
    y = "Pup June Mass (g)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))