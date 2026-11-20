
library(tidyverse)
library(arrow)
library(PHEindicatormethods)
library(gtsummary)
library(prais)
library(segmented)
library(scales)
library(patchwork)
library(officer)
library(flextable)
library(readxl)
library(ggbump)   # install.packages("ggbump")
# Prep -----
source("helper_fun.R")
study_years <- 2005:2024
sepsis_label <- "A41"
regions_num <- c("1", "2", "3", "4", "5")
regions_all <- c("0", regions_num)

age_levels <- c(
  "0", "1_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
  "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69",
  "70_74", "75_79", "80"
)
region_levels <- c("North", "Northeast", "Southeast", "South", "Central-west")
region_ord <- c("Brazil", region_levels)
race_levels <- c("White", "Black", "Mixed", "Yellow", "Indigenous", "Missing")
region_lbl <- c(
  "0" = "Brazil", "1" = "North", "2" = "Northeast",
  "3" = "Southeast", "4" = "South", "5" = "Central-west"
)
region_colors <- c(
  Brazil = "#111111", 
  North = "#5E2D30FF",
  Northeast = "#008E90FF", 
  Southeast = "#1C77A3FF",
  South = "#C5A387FF",
  `Central-west` = "#E9D097FF"
)
comorb_colors <- c(
  "dm" = "#E69F00FF",
  "ckd" = "#56B4E9FF",
  "has" = "#009E73FF"
)
X_BY <- 3
PT_SIZE <- 1.3
LN_SIZE <- 1.2
RIBBON_A <- 0.12
# ---- 1. Load data -----------------------------------------------------------
dados_comorb <- read_parquet("dados_comorb0524_151126.parquet")
br_pop_years <- read_parquet("br_pop_years.parquet")
brpop_2010 <- read_excel("br_pop_2010.xlsx")


# ---- 2. Transformations ----------------------------------------------------
dados_comorb <- dados_comorb |>
  mutate(
    age_group = case_when(
      age < 1 ~ "0",     age < 5 ~ "1_4",   age < 10 ~ "5_9",
      age < 15 ~ "10_14", age < 20 ~ "15_19", age < 25 ~ "20_24",
      age < 30 ~ "25_29", age < 35 ~ "30_34", age < 40 ~ "35_39",
      age < 45 ~ "40_44", age < 50 ~ "45_49", age < 55 ~ "50_54",
      age < 60 ~ "55_59", age < 65 ~ "60_64", age < 70 ~ "65_69",
      age < 75 ~ "70_74", age < 80 ~ "75_79", age >= 80 ~ "80"
    ),
    age_bin = if_else(age < 60, "<60", "\u226560"), # 2 faixas -> Tabela 1
    sex = if_else(SEXO == 2, "female", "male")
  )
dados_comorb <- dados_comorb |> 
  mutate(
    
    immediate_cause_clean = if_else(immediate_cause_grouped %in% c("E10","E11","E14"),"DM",
                                    immediate_cause_grouped)
  )
# Presentation copy for the gtsummary tables (labels + ordering).
dados_tab <- dados_comorb |>
  mutate(
    age_group = factor(age_group, levels = age_levels),
    sex = factor(sex, levels = c("female", "male")),
    region_lbl = factor(
      recode(as.character(region),
        `1` = "North",     `2` = "Northeast",
        `3` = "Southeast", `4` = "South",
        `5` = "Central-west"
      ),
      levels = region_levels
    ),
    race_ethnicity = forcats::fct_na_value_to_level(
      as.factor(race_ethnicity), "Missing"
    ),
    race_ethnicity = factor(race_ethnicity, levels = race_levels),
    number_commorbidities = factor(number_commorbidities, levels = c(1, 2, 3))
  )


dados_comorb <- dados_comorb |>
  mutate(
    all = 1
  )
# ---- 3. Age-standardized mortality rates (DSR) -----------------------------
# Reference population: IBGE 2010 census
std_2010 <- brpop_2010 |>
  pivot_longer(cols = c(male, female), names_to = "sex", values_to = "std_pop") |>
  dplyr::select(age_group, sex, std_pop)

# One DSR run for a given death-row filter; returns Brazil + 5 regions.
study_years <- 2005:2024

# população: 1 linha por estrato, Brasil = soma das regiões
pop_reg <- br_pop_years |>
  rename(region = regiao, year_death = ANO) |>
  mutate(region = as.character(region)) |>
  filter(region %in% regions_num, year_death %in% study_years) |>
  group_by(region, age_group, sex, year_death, pop) |>
  summarise(pop = sum(pop), .groups = "drop") # just sorted/organised dt
pop_grid <- bind_rows(
  pop_reg,
  pop_reg |> group_by(age_group, sex, year_death) |>
    summarise(pop = sum(pop), .groups = "drop") |> mutate(region = "0")
)

# Tables ####
# Table 1 #####
table1_final <- tbl_merge(
  tbls = list(
    tbl1_base(dados_tab |> 
                mutate(group="Underlying cause")),
    tbl1_base(dados_tab |> filter(has == 1) |> 
                mutate(
                  group = if_else(substr(CAUSABAS, 1, 3) %in% c("I10"),
                                  "Underlying cause", "Elsewhere on certificate"
                  )
                )),
    tbl1_base(dados_tab |> filter(dm == 1) |> 
                mutate(group = if_else(substr(CAUSABAS, 1, 3) %in% c("E10","E11","E14"),
                                       "Underlying cause", "Elsewhere on certificate"
                ))
    ),
    tbl1_base(dados_tab |> filter(ckd == 1) |> 
                mutate(group = if_else(substr(CAUSABAS, 1, 3) %in% c("N18"),
                                       "Underlying cause", "Elsewhere on certificate"
                ))
    )
  ),
  tab_spanner = c(
    "**Overall**", "**Hypertension (I10)**",
    "**Diabetes mellitus (E10-E14)**",
    "**Chronic kidney disease (N18)**"
  )
)


# Table 2 #####
table2_final <- tbl_merge(
  tbls = list(
    tbl2_base(dados_tab),
    tbl2_base(dados_tab |> filter(has == 1)),
    tbl2_base(dados_tab |> filter(dm == 1)),
    tbl2_base(dados_tab |> filter(ckd == 1))
  ),
  tab_spanner = c(
    "**Overall**", "**Hypertension (I10)**",
    "**Diabetes mellitus (E10-E14)**",
    "**Chronic kidney disease (N18)**"
  )
)

## Male####

table2_final_male <- tbl_merge(
  tbls = list(
    tbl2_base(dados_tab |> filter(sex=="male")),
    tbl2_base(dados_tab|> filter(sex=="male") |> filter(has == 1)),
    tbl2_base(dados_tab|> filter(sex=="male") |> filter(dm == 1)),
    tbl2_base(dados_tab|> filter(sex=="male") |> filter(ckd == 1))
  ),
  tab_spanner = c(
    "**Overall**", "**Hypertension (I10)**",
    "**Diabetes mellitus (E10-E14)**",
    "**Chronic kidney disease (N18)**"
  )
)
## Female ####
table2_final_female <- tbl_merge(
  tbls = list(
    tbl2_base(dados_tab |> filter(sex=="female")),
    tbl2_base(dados_tab|> filter(sex=="female") |> filter(has == 1)),
    tbl2_base(dados_tab|> filter(sex=="female") |> filter(dm == 1)),
    tbl2_base(dados_tab|> filter(sex=="female") |> filter(ckd == 1))
  ),
  tab_spanner = c(
    "**Overall**", "**Hypertension (I10)**",
    "**Diabetes mellitus (E10-E14)**",
    "**Chronic kidney disease (N18)**"
  )
)
dsr_all <- bind_rows( 
  dsr_comorb_sepsis("dm") |> mutate(comorbidity = "Diabetes mellitus"),
  dsr_comorb_sepsis("has") |> mutate(comorbidity = "Hypertension"),
  dsr_comorb_sepsis("ckd") |> mutate(comorbidity = "Chronic kidney disease"),
  dsr_comorb_sepsis("all") |> mutate(comorbidity = "Any")
)

dsr_all |> 
  mutate(
    rate = paste0(round(value,2)," (",round(lowercl,2), " to ", round(uppercl,2),")")
  )|> clipr::write_clip()

# before 2020
apc_agestd_bef <- dsr_all |> filter(year_death<2020) |> 
  ungroup() |>
  group_by(comorbidity, region) |>
  group_modify(~ apc_from_rate(.x)) |>
  ungroup()
apc_agestd_bef <- apc_agestd_bef |> 
  mutate(
    apc = paste0(round(APC,2)," (",round(APC_lower,2), " to ", round(APC_upper,2),")")
  )
apc_agestd_bef |> clipr::write_clip()

# all period

apc_agestd <- dsr_all |>
  ungroup() |>
  group_by(comorbidity, region) |>
  group_modify(~ apc_from_rate(.x)) |>
  ungroup()
apc_agestd <- apc_agestd |> 
  mutate(
    apc = paste0(round(APC,2)," (",round(APC_lower,2), " to ", round(APC_upper,2),")")
  )
apc_agestd |> clipr::write_clip()




# ---- Figure 1: #####

prop_ill <- bind_rows(
  underlying_prop_all(dados_comorb |> filter(has == 1), "Ill defined") |>
    mutate(cause = "has"),
  underlying_prop_all(dados_comorb |> filter(dm == 1), "Ill defined") |>
    mutate(cause = "dm"),
  underlying_prop_all(dados_comorb |> filter(ckd == 1), "Ill defined") |>
    mutate(cause = "ckd")
) |>
  ggplot(aes(year_death, prop, color = cause)) +
  geom_line(
    linewidth = LN_SIZE,
    alpha = 0.7
  ) +
  geom_point(
    size = PT_SIZE,
    alpha = 0.7
  ) +
  geom_point(
    size = PT_SIZE, stroke = .2, shape = 1, color = "black"
  ) +
  scale_color_manual(
    values = comorb_colors,
    labels = c(
      "dm" = "Diabetes\nMellitus",
      "ckd" = "Chronic Kidney\nDisease",
      "has" = "Hypertension"
    )
  ) +
  scale_x_year() +
  scale_y_continuous(
    limits = c(0, .4),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "",
    color = "Commorbidity",
    y = "Proportion of ill-defined as\nimmediate cause of death"
  ) +
  guides(color = guide_legend(override.aes = list(
    linetype = 0,
    size = 2
  ))) +
  theme_paper()
prop_ill

prop_sepsis <- bind_rows(
  underlying_prop_all(dados_comorb |> filter(has == 1), "A41") |>
    mutate(cause = "has"),
  underlying_prop_all(dados_comorb |> filter(dm == 1), "A41") |>
    mutate(cause = "dm"),
  underlying_prop_all(dados_comorb |> filter(ckd == 1), "A41") |>
    mutate(cause = "ckd")
) |>
  ggplot(aes(year_death, prop, color = cause)) +
  geom_line(
    linewidth = LN_SIZE,
    alpha = 0.7
  ) +
  geom_point(
    size = PT_SIZE,
    alpha = 0.7
  ) +
  geom_point(
    size = PT_SIZE, stroke = .2, shape = 1, color = "black"
  ) +
  scale_color_manual(
    values = comorb_colors,
    labels = c(
      "dm" = "Diabetes\nMellitus",
      "ckd" = "Chronic Kidney\nDisease",
      "has" = "Hypertension"
    )
  ) +
  scale_x_year() +
  scale_y_continuous(
    limits = c(0, .4),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "Year",
    color = "Commorbidity",
    y = "Proportion of sepsis as \nimmediate cause of death"
  ) +
  guides(color = guide_legend(override.aes = list(
    linetype = 0,
    size = 2
  ))) +
  theme_paper()



dsr_plot <- dsr_all |>
  ungroup() |>
  transmute(
    year = year_death,
    region = factor(region_lbl[as.character(region)], levels = region_ord),
    rate = value, low = lowercl, high = uppercl,
    comorbidity = comorbidity
  )
dsr_plot <- dsr_plot |>
  mutate(label_var = if_else(year == 2024, region, NA_character_))
brazil_dsr <- dsr_plot |> filter(
  region == "Brazil",
  comorbidity != "Any"
)
p_all <- ggplot(brazil_dsr, aes(year, rate, color = comorbidity, group = comorbidity)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, color = NA, fill = "gray70") +
  geom_line(
    linewidth = LN_SIZE,
    alpha = 0.7
  ) +
  geom_point(
    size = PT_SIZE,
    alpha = 0.7
  ) +
  geom_point(
    size = PT_SIZE, stroke = .2, shape = 1, color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Diabetes mellitus" = "#E69F00FF",
      "Chronic kidney disease" = "#56B4E9FF",
      "Hypertension" = "#009E73FF"
    ),
    labels = c(
      "dm" = "Diabetes\nMellitus",
      "ckd" = "Chronic Kidney\nDisease",
      "has" = "Hypertension"
    )
  ) +
  scale_x_year()+
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    limits = c(0, 25)
  ) +
  labs(
    x = NULL,
    color = "",
    y = "Age-standardised mortality\nrate per 100,000 (sepsis)"
  ) +
  guides(color = guide_legend(override.aes = list(
    linetype = 0,
    size = 2
  ))) +
  theme_paper()


stack_panels(list(prop_ill, prop_sepsis,p_all))

# sup figure 1 region ####
pill_has <- prop_ill_region(dados_comorb |> filter(has == 1))
pill_dm <- prop_ill_region(dados_comorb |> filter(dm == 1))
pill_ckd <- prop_ill_region(dados_comorb |> filter(ckd == 1))

pill_ckd/pill_dm/(pill_has & labs(x="Year"))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))



# Figure 2 ####



p_has <- plot_by_region("Hypertension")
p_dm <- plot_by_region("Diabetes mellitus")
p_ckd <- plot_by_region("Chronic kidney disease")

p_ckd/p_dm/(p_has & labs(x="Year"))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))


# Table 3 - full-period APC (wide by comorbidity) #####
apc_wide <- apc_agestd |>
  mutate(
    Region = factor(region_lbl[as.character(region)], levels = region_ord),
    cell = sprintf("%.2f (%.2f to %.2f)", APC, APC_lower, APC_upper)
  ) |>
  dplyr::select(comorbidity, Region, cell) |>
  pivot_wider(names_from = comorbidity, values_from = cell) |>
  arrange(Region)
apc_ft <- flextable(apc_wide) |>
  fontsize(size = 9, part = "all") |>
  add_header_row(
    values = c("", "Age-standardized APC (%, 95% CI)"),
    colwidths = c(1, ncol(apc_wide) - 1)
  ) |>
  set_table_properties(width = 1, layout = "autofit")

# sup table 1
apc_wide_bef <- apc_agestd_bef |>
  mutate(
    Region = factor(region_lbl[as.character(region)], levels = region_ord),
    cell = sprintf("%.2f (%.2f to %.2f)", APC, APC_lower, APC_upper)
  ) |>
  dplyr::select(comorbidity, Region, cell) |>
  pivot_wider(names_from = comorbidity, values_from = cell) |>
  arrange(Region)
apc_bef_ft <- flextable(apc_wide_bef) |>
  fontsize(size = 9, part = "all") |>
  add_header_row(
    values = c("", "Age-standardized APC (%, 95% CI)"),
    colwidths = c(1, ncol(apc_wide) - 1)
  ) |>
  set_table_properties(width = 1, layout = "autofit")



dsr_wide <- dsr_all |>
  ungroup() |>
  transmute(
    Comorbidity = comorbidity, Year = as.integer(year_death),
    region = factor(region_lbl[as.character(region)], levels = region_ord),
    cell = sprintf("%.1f (%.1f-%.1f)", value, lowercl, uppercl)
  ) |>
  pivot_wider(names_from = region, values_from = cell) |>
  arrange(Comorbidity, Year)
dsr_ft <- flextable(dsr_wide) |>
  fontsize(size = 8, part = "all") |>
  set_table_properties(width = 1, layout = "autofit") |>
  merge_v(j = "Comorbidity")









rank_df <- dados_comorb |>
  filter(!immediate_cause_grouped %in% c("Other", "Ill defined")) |>
  count(year_death, cause = immediate_cause_grouped) |>
  mutate(rank = rank(-n, ties.method = "first"), .by = year_death) |>
  filter(rank <= 5) |>
  mutate(periodo_num = as.numeric(factor(year_death))) |>
  arrange(cause, periodo_num) |>
  # new segment whenever a cause leaves and re-enters the top 5
  mutate(seg = cumsum(periodo_num - lag(periodo_num, default = first(periodo_num) - 1) > 1),
         grp = paste(cause, seg),
         .by = cause)

# geom_bump needs >= 2 points per line
bump_df <- rank_df |> filter(n() > 1, .by = grp)

# label each segment at its last point
lab_df <- rank_df |> slice_max(periodo_num, by = grp)

causes <- sort(unique(rank_df$cause))

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "grey75",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")

pal <- setNames(okabe_ito[seq_along(causes)], causes)

ggplot(rank_df, aes(x = periodo_num, y = rank, colour = cause)) +
  geom_bump(data = bump_df, aes(group = grp), linewidth = 1.2, smooth = 8) +
  geom_point(size = 4) +
  geom_text(data = lab_df, aes(label = cause),
            hjust = 0, nudge_x = 0.15, size = 3) +
  scale_y_reverse(breaks = 1:5) +
  scale_x_continuous(breaks = unique(rank_df$periodo_num),
                     labels = levels(factor(rank_df$year_death)),
                     expand = expansion(add = c(0.3, 1.2))) +
  labs(x = NULL, y = "Rank") +
  scale_colour_manual(values = pal)+
  theme_minimal() +
  theme(legend.position = "none", panel.grid.minor = element_blank())


