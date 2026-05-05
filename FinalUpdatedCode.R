# Combined code of Alex Kelly and Jamil Fayad
# Maternal FROH vs. pup mass traits
# Our main question is: Does maternal inbreeding having a negative effect on mass of the pups
# June mass is the emergence mass
# August mass is end of season mass 
# Seasonal mass gained is August mass - June mass

# --------------
# Set up and loading the necessary packages 
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
# For regressions:
library(lmerTest) 
library(lme4)
library(patchwork) # for creating multipanels
# for AIC comparisons: 
library(MuMIn)  
library(kableExtra)
library(knitr)
getwd() # making sure that we are in the correct directory 

# Total genome size in base pairs, used to calculate FROH
genome_size <- 2319465413

# --------------
# Importing the ROH data from ADROIT (computer cluster) and cleaning it up

# Read bcftools ROH output (our filtered marmot library)
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

# --------------
# Importing our metadata 

# Parentage-informed pedigree
ped_data <- read.csv("marmot_Parentage_Informed_Pedigree.csv")

# Mass metadata
mass_metadata <- read_excel("marmot_mass_metadata.xlsx")

# Keep only variables needed for analysis
meta_mass <- mass_metadata %>%
  select(uid, year, massjun, massaug, nb_mass) %>%
  mutate(mass_gain = massaug - massjun)

# --------------
# Matching the pups to their mothers --> using SNP Pedigree that was constructed in ADROIT 

# Extract pup and mother IDs from pedigree file
pup_mom <- ped_data %>%
  select(id, dam) %>%
  filter(!is.na(dam), dam != "") %>%
  transmute(
    pup_uid = str_extract(id, "[0-9]+_[0-9]+"), # filter the IDs to have the correct format in order to compare to each other
    mom_uid = str_extract(dam, "[0-9]+_[0-9]+")
  ) %>%
  filter(!is.na(pup_uid), !is.na(mom_uid))%>%
  filter(mom_uid != pup_uid)

# Attach each mother's FROH value to her pup
combined_froh_id <- pup_mom %>%
  left_join(roh_clean, by = c("mom_uid" = "uid")) %>%
  rename(mom_FROH = FROH)

# --------------
# Constructing the analysis data set in order to compare maternal FROH to pup masses 
analysis_mass_froh <- combined_froh_id %>%
  inner_join(meta_mass, by = c("pup_uid" = "uid")) %>%
  filter(
    !is.na(mom_FROH),
    !is.na(massjun),
    !is.na(massaug),
    !is.na(mass_gain),
    !is.na(nb_mass),
    !is.na(year)
  )

# --------------
# Running the regressions, linear mixed models, and the AIC comparisons

# AIC analysis tells us if the linear mixed models accounting for year as random effect is a better model than the regular linear model between mass variables and maternal FROH
# If the delta AICc is less than 2 then there is no difference between the linear mixed model and our regular basic linear model

### Basic linear model of mass gain with just FROH as fixed effect
base_mod <- lm(mass_gain ~ mom_FROH, data = analysis_mass_froh)
summary(base_mod) # p-value: 0.149

# Linear mixed effects model of mass gain with FROH as fixed effect and year born as random effect
mod_yrborn <- lmer(mass_gain ~ mom_FROH + (1 | year), data = analysis_mass_froh, REML = FALSE)

# Conduct AIC analyses on mass gain 
aic_dataframe <- as.data.frame(model.sel(base_mod, mod_yrborn))

# Keep only useful columns
aic_dataframe <- aic_dataframe[, c("df", "logLik", "AICc", "delta", "weight")]

aic_dataframe$AICc <- round(aic_dataframe$AICc, 1)
aic_dataframe$delta <- round(aic_dataframe$delta, 3)

# Rename columns and models
colnames(aic_dataframe) <- c("K", "logLik", "AICc", "Delta AICc", "Weight")

aic_dataframe$Model <- c("Base model + year born (random)",
                 "Base model: Mass gain ~ FROH") 

# Order by best AICc value 
aic_dataframe <- aic_dataframe[order(aic_dataframe$AICc), ]

# Table formatting
kable(aic_dataframe[, c("Model", "K", "AICc", "Delta AICc")]) %>% 
  kable_styling(position = "left", full_width = FALSE)
# 0.797 = delta AIC when comparing the normal linear model of FROH against mass gain with the mixed model accounting for year 

# Best model summary
summary(mod_yrborn) # p = 0.261

### Basic linear regression comparing mass in august with mom FROH
base_linear_aug <- lm(massaug ~ mom_FROH, data = analysis_mass_froh)
summary(base_linear_aug) # p-value: 0.129

# Linear mixed model to compare August mass with FROH with year as random effect
aug_ml_year <- lmer(massaug ~ mom_FROH + (1 | year), data = analysis_mass_froh, REML = FALSE)

# Conduct AIC analyses for August mass
aic_dataframe_aug <- as.data.frame(model.sel(base_linear_aug, aug_ml_year))

# Keep only useful columns
aic_dataframe_aug <- aic_dataframe_aug[, c("df", "logLik", "AICc", "delta", "weight")]

aic_dataframe_aug$AICc <- round(aic_dataframe_aug$AICc, 1)
aic_dataframe_aug$delta <- round(aic_dataframe_aug$delta, 3)

# Rename columns and models
colnames(aic_dataframe_aug) <- c("K", "logLik", "AICc", "Delta AICc", "Weight")

aic_dataframe_aug$Model <- c("Base model: August mass ~ FROH",
                             "Base model + year born (random)")

# Order by best AICc value 
aic_dataframe_aug <- aic_dataframe_aug[order(aic_dataframe_aug$AICc), ]

# Table formatting
kable(aic_dataframe_aug[, c("Model", "K", "AICc", "Delta AICc")]) %>% 
  kable_styling(position = "left", full_width = FALSE)
# 1.844 = delta AIC when comparing the normal linear model of FROH against August mass with the mixed model accounting for year 

summary(aug_ml_year) # p = 0.172

### Basic linear regression comparing mass in June against with mom FROH
base_linear_jun <- lm(massjun ~ mom_FROH, data = analysis_mass_froh)
summary(base_linear_jun) # p-value: 0.887

# Linear mixed model to compare June mass with FROH with year as random effect
jun_ml_year <- lmer(massjun ~ mom_FROH + (1 | year), data = analysis_mass_froh, REML = FALSE)

# Conduct AIC analyses for June mass
aic_dataframe_jun <- as.data.frame(model.sel(base_linear_jun, jun_ml_year))

# Keep only useful columns
aic_dataframe_jun <- aic_dataframe_jun[, c("df", "logLik", "AICc", "delta", "weight")]

aic_dataframe_jun$AICc <- round(aic_dataframe_jun$AICc, 1)
aic_dataframe_jun$delta <- round(aic_dataframe_jun$delta, 3)

# Rename columns and models
colnames(aic_dataframe_jun) <- c("K", "logLik", "AICc", "Delta AICc", "Weight")

aic_dataframe_jun$Model <- c("Base model + year born (random)",
                         "Base model: June mass ~ FROH") 

# Order by best AICc value 
aic_dataframe_jun <- aic_dataframe_jun[order(aic_dataframe_jun$AICc), ]

# Table formatting
kable(aic_dataframe_jun[, c("Model", "K", "AICc", "Delta AICc")]) %>% 
  kable_styling(position = "left", full_width = FALSE)
# 3.505 = delta AIC when comparing the normal linear model of FROH against June mass with the mixed model accounting for year

summary(jun_ml_year) # p = 0.993

# --------------
# Running distributions checks to see if we have parametric or non-parametric
# Shapiro-Wilk tests
shapiro.test(analysis_mass_froh$mom_FROH)   # p = 7.37e-06 (non-normal)
shapiro.test(analysis_mass_froh$mass_gain)  # p = 0.480 (normal)
shapiro.test(analysis_mass_froh$massjun)    # p = 0.044 (non-normal)
shapiro.test(analysis_mass_froh$massaug)    # p = 0.220 (normal)

# Histograms
par(mfrow = c(1, 4))
hist(analysis_mass_froh$mom_FROH, main = "Mom FROH", xlab = "Mom FROH", col = "lightblue")
hist(analysis_mass_froh$mass_gain, main = "Mass Gain", xlab = "Mass Gain (g)", col = "lightgreen")
hist(analysis_mass_froh$massjun, main = "June Mass", xlab = "June Mass (g)", col = "tomato")
hist(analysis_mass_froh$massaug, main = "August Mass", xlab = "August Mass (g)", col = "khaki")

# QQ plots
par(mfrow = c(1, 4))
for (var in c("mom_FROH", "mass_gain", "massjun", "massaug")) {
  qqnorm(analysis_mass_froh[[var]], main = paste("QQ:", var))
  qqline(analysis_mass_froh[[var]], col = "red")
}
par(mfrow = c(1, 1))

# --------------
# Running correlation tests
# Spearman Rank Correlation tests
sp_gain <- cor.test(analysis_mass_froh$mom_FROH, analysis_mass_froh$mass_gain, method = "spearman", exact = FALSE)
sp_aug  <- cor.test(analysis_mass_froh$mom_FROH, analysis_mass_froh$massaug, method = "spearman", exact = FALSE)
sp_jun  <- cor.test(analysis_mass_froh$mom_FROH, analysis_mass_froh$massjun, method = "spearman", exact = FALSE)

sp_gain # rho = -0.328, p = 0.017
sp_aug # rho = -0.272, p = 0.049
sp_jun # rho = 0.096, p = 0.496

# Pearson (check for mass_gain and massaug, which were normally distributed)
cor.test(analysis_mass_froh$mom_FROH, analysis_mass_froh$mass_gain, method = "pearson")  # p = 0.149
cor.test(analysis_mass_froh$mom_FROH, analysis_mass_froh$massaug, method = "pearson")    # p = 0.129

# --------------
# Making multipanel figures
shared_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank()
  )

# Panel A: Mass Gain
massgain_panel <- ggplot(analysis_mass_froh, aes(x = mom_FROH, y = mass_gain)) +
  geom_point(alpha = 0.6, color = "#2c7bb6", size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#d7191c", linewidth = 0.8) +
  annotate("text",
           x = max(analysis_mass_froh$mom_FROH, na.rm = TRUE),
           y = max(analysis_mass_froh$mass_gain, na.rm = TRUE),
           label = paste0("Spearman rho = ", round(sp_gain$estimate, 3),
                          "\nModel: linear model(Mass gain ~ maternal FROH)",
                          "\nModel p = 0.149"),
           hjust = 1, vjust = 1, size = 3.2, fontface = "italic"
  ) +
  labs(title = "A) Seasonal mass gain",
       subtitle = "Seasonal mass gain = August mass - June mass",
       x = expression(Maternal ~ F[ROH]),
       y = "Mass Gain (g)") +
  shared_theme

# Panel B: August Mass
augustmass_panel <- ggplot(analysis_mass_froh, aes(x = mom_FROH, y = massaug)) +
  geom_point(alpha = 0.6, color = "#2c7bb6", size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#d7191c", linewidth = 0.8) +
  annotate("text",
           x = max(analysis_mass_froh$mom_FROH, na.rm = TRUE),
           y = max(analysis_mass_froh$massaug, na.rm = TRUE),
           label = paste0("Spearman rho = ", round(sp_aug$estimate, 3),
                          "\nModel: linear model(August mass ~ maternal FROH)",
                          "\nModel p = 0.129"),
           hjust = 1, vjust = 1, size = 3.2, fontface = "italic"
  ) +
  labs(title = "B) August mass",
       x = expression(Maternal ~ F[ROH]),
       y = "August Mass (g)") +
  shared_theme

# Combine into multi-panel
combined_fig <- massgain_panel + augustmass_panel +
  plot_annotation(
    title = "Effect of Maternal FROH on Pup Mass Traits",
    subtitle = "Yellow-bellied marmots (Marmota flaviventris), n = 53 pups",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray40")
    )
  ) +
  plot_layout(ncol = 2)

combined_fig

# June mass against maternal FROH was not included in the figure because it was highly insignificant. 
# The correlation tests along with the regressions showcase that there is no real correlation between maternal FROH and June mass 
# even when accounting for year born. The mass gain and August mass when compared to maternal FROH show some promise because 
# of the linear trend, and the spearman correlation results. Overall, based on our data the pups emergence mass (June mass) 
# is not correlated with maternal FROH and therefore is excluded from our figures and results. 

