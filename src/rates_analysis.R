# load dependencies
library(tidyverse)
library(lme4)


# read in dataset
dat <- read_csv("data/Output/Out.csv",
    col_types = list(Is_Synonymous = col_factor(), Px = col_factor())
)

# filter out fluctuations
dat <- filter(dat, Pr_fixation > 0)

# evo rate by synonymity
evo_by_syn_m <- lmer(Evo_rate ~ Is_Synonymous + (1 | Px), dat)

evo_by_syn_plot <- ggplot(dat, aes(x = as.factor(Is_Synonymous), y = Evo_rate)) +
    geom_boxplot() +
    xlab("cyl")
