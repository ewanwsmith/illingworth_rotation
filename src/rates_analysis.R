# load dependencies
library(readr)
library(tidyverse)
library(lme4)
library(ggplot2)
library(ggpubr)
library(viridis)
library(scales)

# read in variants dataset
dat <- read_csv("data/Output/Out.csv",
    col_types = list(Is_Synonymous = readr::col_factor(), Px = readr::col_factor(), Protein = readr::col_factor())
)

# add Px level for Px with no mutations
dat$Px <- factor(dat$Px, levels = c("data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                                    "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                                    "data/CAMP004884", "data/CAMP007136", "data/Kemp"))
dat$Px <- str_remove(dat$Px, "data/")

#convert Px name to Px number for cleaner plotting
pxn <- read_csv("data/reference/Px_Names_to_n.csv")
dat <- left_join(dat, pxn, join_by(Px == Px_Name))
dat$Px_n <- as.factor(dat$Px_n)

#create patient colour palette for consistency between plots
PxPallete <- c("2" = viridis(8, alpha = 0.8, direction = -1)[1],
               "3" = viridis(8, alpha = 0.8, direction = -1)[2],
               "4" = viridis(8, alpha = 0.8, direction = -1)[3],
               "5" = viridis(8, alpha = 0.8, direction = -1)[4],
               "6" = viridis(8, alpha = 0.8, direction = -1)[5],
               "7" = viridis(8, alpha = 0.8, direction = -1)[6],
               "8" = viridis(8, alpha = 0.8, direction = -1)[7],
               "9" = viridis(8, alpha = 0.8, direction = -1)[8])


# filter out fluctuations
dat <- tibble(dat)
dat <- filter(dat, Pr_fixation > 0)


# mutations by Px
Px_count <- dat %>% group_by(Px_n, Is_Synonymous, .drop = FALSE) %>% count()
Px_count <- drop_na(Px_count)

mutations_by_px_plot = ggplot(Px_count, aes(fill=Is_Synonymous, y = n, x = as.factor(Px_n))) + 
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Patient", y = "Variant Count") +
  theme_linedraw() +
  scale_y_continuous(breaks = pretty(Px_count$n, n = 10)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.8, name = "", labels = c("Synonymous", "Non-synonymous")) +
  theme(legend.position = "bottom")
plot(mutations_by_px_plot)


# evo rate by synonymity
evo_by_syn_m <- lmer(Evo_rate ~ Is_Synonymous + (1 | Px), dat)

evo_by_syn_plot <- ggplot(dat, aes(x = as.factor(Is_Synonymous), y = Evo_rate, fill = Is_Synonymous)) +
  geom_violin(width=0.4, alpha = 0.8) +
  geom_boxplot(width=0.1, color="black", alpha = 0) +
  geom_point() +
  scale_x_discrete(labels = c('Non-Synonymous','Synonymous')) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "Mean rate of evolution") +
  ylim(0, 0.45) +
  stat_compare_means(aes(label = ..p.signif..), label.x = 1.5, label.y = 0.35) +
  theme_linedraw() +
  theme(legend.position = "none")
plot(evo_by_syn_plot)

# variants within genome

# read in positions dataset
pos <- read_csv("data/reference/protein_positions.csv",
                col_types = list(Protein = readr::col_factor())
)

dat$Protein = factor(dat$Protein, levels = levels(pos$Protein))

# join pos to dat for ordering by genome position
pos <- rename(pos, Position_range = Position)
dat <- left_join(dat, pos, join_by(Protein))

# where are the variants ? 
# position 
position_plot <- dat %>% 
  ggplot(aes(y = as.factor(Px_n), x = Position, color = Px_n)) +
  geom_point() +
  scale_y_discrete(drop = FALSE) +
  scale_color_manual(values = PxPallete) +
  labs(y = "Patient") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ggtitle("All variants")
plot(position_plot)

s_position_plot <- dat %>% filter(dat$Is_Synonymous == "Yes") %>%
  ggplot(aes(y = as.factor(Px_n), x = Position, color = Px_n)) +
  geom_point() +
  scale_y_discrete(drop = FALSE) +
  scale_color_manual(values = PxPallete) +
  labs(y = "Patient") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ggtitle("Synonymous variants")
plot(s_position_plot)

ns_position_plot <- dat %>% filter(dat$Is_Synonymous == "No") %>%
  ggplot(aes(y = as.factor(Px_n), x = Position, color = Px_n)) +
  geom_point() +
  scale_y_discrete(drop = FALSE) +
  scale_color_manual(values = PxPallete) +
  labs(y = "Patient") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ggtitle("Non-synonymous variants")
plot(ns_position_plot)

# protein
# all variants
Protein_by_px_count <- dat %>% group_by(Px_n, Protein, .drop = FALSE) %>% count()
Protein_by_px_count <- left_join(Protein_by_px_count, pos, join_by(Protein))


protein_by_px_plot <- Protein_by_px_count %>%
  arrange(Position_range) %>%
  mutate(Protein = factor(Protein, levels = Protein)) %>% 
  ggplot(aes(fill = as.factor(Px_n), y = n, x = forcats::fct_rev(Protein))) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = PxPallete, name = "Patient") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(x = "Protein", y = "Variant count") +
  theme(legend.position = "right", legend.direction = "vertical") +
  coord_flip() +
  ggtitle("All variants")
plot(protein_by_px_plot)

# synonymous variants
s <- filter(dat, Is_Synonymous == "Yes")
Protein_by_px_count_s <- s %>% group_by(Px_n, Protein, .drop = FALSE) %>% count()
Protein_by_px_count_s <- left_join(Protein_by_px_count_s, pos, join_by(Protein))

protein_by_px_plot_s <- Protein_by_px_count_s %>%
  mutate(Protein = factor(Protein, levels = Protein)) %>%
  mutate(Px_n = factor(Px_n, levels = 1:9)) %>% 
  ggplot(aes(fill = as.factor(Px_n), y = n, x = forcats::fct_rev(Protein))) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = PxPallete, name = "Patient") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0:10)) +
  labs(x = "Protein", y = "Variant count") +
  theme(legend.position = "right", legend.direction = "vertical") +
  coord_flip() +
  ggtitle("Synonymous variants")
plot(protein_by_px_plot_s)

# non-synonymous variants
ns <- filter(dat, Is_Synonymous == "No")
Protein_by_px_count_ns <- ns %>% group_by(Px_n, Protein, .drop = FALSE) %>% count()
Protein_by_px_count_ns <- left_join(Protein_by_px_count_ns, pos, join_by(Protein))

protein_by_px_plot_ns <- Protein_by_px_count_ns %>%
  mutate(Protein = factor(Protein, levels = Protein)) %>%
  mutate(Px_n = factor(Px_n, levels = 1:9)) %>% 
  ggplot(aes(fill = as.factor(Px_n), y = n, x = forcats::fct_rev(Protein))) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = PxPallete, name = "Patient") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0:10)) +
  labs(x = "Protein", y = "Variant count") +
  theme(legend.position = "right", legend.direction = "vertical") +
  coord_flip() +
  ggtitle("Non-synonymous variants")
plot(protein_by_px_plot_ns)


                  