# =============================================================================
# Figure 2 — Publications, pathogens, species × database, world map
# Panels A–D assembled with patchwork
# Output: Figure2.pdf / Figure2.png
# =============================================================================

library(countrycode)
library(dplyr)
library(forcats)
library(ggnewscale)
library(ggplot2)
library(maps)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(scales)
library(stringr)
library(tidyr)

path_xlsx <- "Supp_Appendix_1_final.xlsx"

df_review <- read_excel(path_xlsx)
names(df_review) <- make.names(trimws(names(df_review)), unique = TRUE)

df_review$Country[df_review$Country == "Scotland"]        <- "UK"
df_review$Country[df_review$Country == "Nothern Ireland"] <- "UK"

# ── Panel A: publications by period and species ───────────────────────────────

data_pub <- df_review %>%
  mutate(
    year       = as.numeric(Year),
    year_range = cut(year,
                     breaks = c(2005, 2010, 2015, 2020, 2025),
                     right  = FALSE,
                     labels = c("2005–2009", "2010–2014",
                                "2015–2019", "2020–2024")),
    Species_plot = case_when(
      is.na(Species) | str_trim(Species) == "" ~ NA_character_,
      str_detect(Species, "\\*")               ~ "Multiple species",
      Species == "Ducks"                       ~ "Poultry",
      TRUE                                     ~ str_trim(Species)
    )
  ) %>%
  filter(!is.na(year_range), !is.na(Species_plot)) %>%
  count(year_range, Species_plot, name = "n")

species_order <- data_pub %>%
  group_by(Species_plot) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(total) %>%
  pull(Species_plot)

data_pub <- data_pub %>%
  mutate(Species_plot = factor(Species_plot, levels = species_order))

n_sp  <- length(levels(data_pub$Species_plot))
pal_sp <- scales::alpha(brewer.pal(max(3, min(n_sp, 8)), "Set3"), 0.8)

plot_pub <- ggplot(data_pub, aes(x = year_range, y = n, fill = Species_plot)) +
  geom_col(width = 0.65, color = "gray30", linewidth = 0.3) +
  scale_fill_manual(values = pal_sp, name = "Animals") +
  labs(x = "Publication year", y = "Number of studies") +
  theme_bw() +
  theme(
    axis.text          = element_text(color = "gray25"),
    legend.text        = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    axis.ticks         = element_blank()
  )

# ── Panel B: pathogens lollipop ───────────────────────────────────────────────

pathogen_count <- df_review %>%
  select(Pathogen) %>%
  mutate(Pathogen = if_else(Pathogen %in% c("Non-specified", "Non-specific", "Unspecified", NA_character_),
                            "Non pathogen-specific", Pathogen)) %>%
  filter(!is.na(Pathogen), Pathogen != "") %>%
  separate_rows(Pathogen, sep = "\\*") %>%
  mutate(Pathogen = str_trim(Pathogen)) %>%
  filter(Pathogen != "") %>%
  group_by(Pathogen) %>%
  summarise(n_studies = n(), .groups = "drop") %>%
  arrange(desc(n_studies))

plot_pat <- ggplot(pathogen_count,
                   aes(x = reorder(Pathogen, n_studies), y = n_studies)) +
  geom_segment(aes(xend = Pathogen, yend = 0, color = n_studies),
               linewidth = 1.2, linetype = "solid", alpha = 0.7) +
  geom_point(aes(color = n_studies), size = 5, alpha = 0.8) +
  coord_flip() +
  labs(x = "Pathogen", y = "Number of studies", color = "Number\nof studies") +
  scale_color_gradient(low = "#3498db", high = "#e74c3c") +
  theme_minimal() +
  theme(
    axis.text.y        = element_text(size = 10),
    axis.text.x        = element_text(size = 11),
    axis.title         = element_text(size = 12, face = "bold"),
    legend.position    = "right",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

# ── Panel C: sunburst (species × database × region) ──────────────────────────

df_donut <- df_review %>%
  mutate(Country = case_when(
    Country %in% c("Scotland", "Nothern Ireland") ~ "UK",
    TRUE ~ Country
  )) %>%
  mutate(Database_category = case_when(
    Database_category == "National system"   ~ "Official tracking database",
    Database_category == "State-managed"     ~ "Official tracking database",
    Database_category == "Production system" ~ "Company-held database",
    TRUE ~ Database_category
  )) %>%
  filter(Database_category != "", !is.na(Database_category)) %>%
  separate_rows(Species, sep = "\\*") %>%
  mutate(Species = str_trim(Species)) %>%
  mutate(
    Country_std   = case_when(Country == "UK"  ~ "United Kingdom",
                              Country == "USA" ~ "United States",
                              TRUE ~ Country),
    Continent_raw = countrycode(Country_std, "country.name", "continent"),
    Continent     = case_when(
      Continent_raw != "Americas"                               ~ Continent_raw,
      Country_std %in% c("Canada", "United States", "Mexico")  ~ "North America",
      TRUE ~ "South America"
    )
  )

top3     <- df_donut %>% count(Species, sort = TRUE) %>% slice_head(n = 3) %>% pull(Species)
df_donut <- df_donut %>%
  mutate(Species = if_else(Species %in% top3, Species, "Other animals"))

keep_cont <- df_donut %>%
  filter(!is.na(Continent)) %>%
  count(Continent, sort = TRUE) %>%
  slice_head(n = 3) %>%
  pull(Continent)

df_donut <- df_donut %>%
  mutate(Continent2 = case_when(
    Continent %in% c("North America", "South America") ~ Continent,
    Continent %in% keep_cont                           ~ Continent,
    TRUE ~ "Other regions"
  ))

make_ring  <- function(dat) dat %>% mutate(ymax = cumsum(value), ymin = lag(ymax, default = 0))
pal_ramp   <- function(n, lo, hi) colorRampPalette(c(lo, hi))(n)

r1    <- df_donut %>% count(Species, name = "value") %>% arrange(Species) %>% make_ring()
offs1 <- r1 %>% select(Species, off = ymin)

r2 <- df_donut %>%
  count(Species, Database_category, name = "value") %>%
  arrange(Species, Database_category) %>%
  group_by(Species) %>%
  mutate(ymax_local = cumsum(value), ymin_local = lag(ymax_local, default = 0)) %>%
  ungroup() %>%
  left_join(offs1, by = "Species") %>%
  mutate(ymin = ymin_local + off, ymax = ymax_local + off)
offs2 <- r2 %>% select(Species, Database_category, off2 = ymin)

r3 <- df_donut %>%
  count(Species, Database_category, Continent = Continent2, name = "value") %>%
  arrange(Species, Database_category, Continent) %>%
  group_by(Species, Database_category) %>%
  mutate(ymax_local = cumsum(value), ymin_local = lag(ymax_local, default = 0)) %>%
  ungroup() %>%
  left_join(offs2, by = c("Species", "Database_category")) %>%
  mutate(ymin = ymin_local + off2, ymax = ymax_local + off2)

sp_levels <- c(sort(top3), "Other animals"); sp_levels <- sp_levels[sp_levels %in% r1$Species]
db_levels <- c("Official tracking database", "Company-held database", "Mixed")
db_levels <- db_levels[db_levels %in% r2$Database_category]
ct_levels <- c("Europe","North America","South America","Asia","Africa","Oceania","Other regions")
ct_levels <- ct_levels[ct_levels %in% r3$Continent]

r1 <- r1 %>% mutate(Species_f = factor(Species, levels = sp_levels))
r2 <- r2 %>% mutate(DB_f      = factor(Database_category, levels = db_levels))
r3 <- r3 %>% mutate(CT_f      = factor(Continent, levels = ct_levels))

sp_cols <- setNames(pal_ramp(length(sp_levels), "mistyrose",  "firebrick"), sp_levels)
db_cols <- setNames(pal_ramp(length(db_levels), "#FAEBD7",    "burlywood4"), db_levels)
ct_cols <- setNames(pal_ramp(length(ct_levels), "lightblue",  "navy"),      ct_levels)

plot_donut <- ggplot() +
  geom_rect(data = r3,
            aes(xmin = 2.0, xmax = 3.0, ymin = ymin, ymax = ymax, fill = CT_f),
            color = "white", linewidth = 0.15) +
  scale_fill_manual(values = ct_cols, name = "Region") +
  new_scale_fill() +
  geom_rect(data = r2,
            aes(xmin = 1.2, xmax = 1.95, ymin = ymin, ymax = ymax, fill = DB_f),
            color = "white", linewidth = 0.15) +
  scale_fill_manual(values = db_cols, name = "Database") +
  new_scale_fill() +
  geom_rect(data = r1,
            aes(xmin = 0.5, xmax = 1.15, ymin = ymin, ymax = ymax, fill = Species_f),
            color = "white", linewidth = 0.15) +
  scale_fill_manual(values = sp_cols, name = "Animals") +
  coord_polar(theta = "y") +
  xlim(0, 3.05) +
  theme_void() +
  theme(legend.position = "right")

# ── Panel D: world map ────────────────────────────────────────────────────────

tab_nb   <- table(df_review$Country)
df_nb    <- data.frame(region = names(tab_nb), nb = as.numeric(tab_nb))
mp       <- map_data("world")
mp_data  <- left_join(mp, df_nb, by = "region")

plot_mapp <- ggplot(mp_data, aes(x = long, y = lat, group = group, fill = nb)) +
  geom_polygon() +
  scale_fill_viridis_c(
    option  = "viridis",
    breaks  = c(1, 5, 10, 15, 20, 25),
    limits  = c(0, max(mp_data$nb, na.rm = TRUE)),
    na.value = "grey80"
  ) +
  labs(fill = "Number of\nstudies") +
  theme_void() +
  theme(
    legend.key   = element_rect(colour = "black", fill = NA),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 12)
  )

# ── Assembly ──────────────────────────────────────────────────────────────────

fig2 <- (plot_pub | plot_pat) / (plot_donut | plot_mapp) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave("Figure2.pdf", fig2, width = 14, height = 10, device = "pdf")
ggsave("Figure2.png", fig2, width = 14, height = 10, dpi = 300)
