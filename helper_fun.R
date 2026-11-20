# Functions

dsr_run <- function(df_deaths) {
  d <- df_deaths |>
    mutate(region = as.character(region)) |>
    filter(region %in% regions_num)
  deaths <- bind_rows(
    d |> count(age_group, region, sex, year_death, name = "deaths"),
    d |> mutate(region = "0") |>
      count(age_group, region, sex, year_death, name = "deaths")
  )
  pop_grid |>
    left_join(deaths, by = c("age_group", "region", "sex", "year_death")) |>
    mutate(deaths = tidyr::replace_na(deaths, 0L)) |>
    left_join(std_2010, by = c("age_group", "sex")) |>
    group_by(year_death, region) |>
    PHEindicatormethods::calculate_dsr(x = deaths, n = pop, stdpop = std_pop)
}
# Wrapper: sepsis (A41) as immediate cause among deaths of one comorbidity.
dsr_comorb_sepsis <- function(flag) {
  dsr_run(dados_comorb |>
            filter(.data[[flag]] == 1, immediate_cause_grouped == sepsis_label))
}

tbl1_base <- function(df) {
  vars <- c("sex", "age", "age_bin", "race_ethnicity", "region_lbl","number_commorbidities","group")
  labs <- list(
    sex ~ "Sex", age ~ "Age (years)",
    age_bin ~ "Age group", # <- era age_group
    race_ethnicity ~ "Race/ethnicity",
    region_lbl ~ "Geographic region",
    number_commorbidities ~ "Number of comorbidities")
  df |>
    dplyr::select(all_of(vars)) |>
    tbl_summary(
      label = labs,
      statistic = list(
        all_continuous() ~ "{median} ({p25}-{p75})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_categorical() ~ c(0, 1)
    )
}
ord_by_freq <- function(x) {
  factor(x, levels = names(sort(table(x), decreasing = TRUE)))
}

tbl2_base <- function(df) {
  df |>
    transmute(cause = ord_by_freq(immediate_cause_clean)) |>
    tbl_summary(
      label = list(cause ~ "Immediate cause"),
      statistic = list(all_categorical() ~ "{n} ({p}%)"),
      digits = all_categorical() ~ c(0, 1),
      missing = "no"
    )
}


apc_from_rate <- function(df) {
  df <- df |>
    filter(is.finite(value), value > 0) |>
    arrange(year_death)
  pw <- prais_winsten(log(value) ~ year_death, data = df, index = "year_death")
  co <- summary(pw)$coefficients
  beta <- co["year_death", "Estimate"]
  se <- co["year_death", "Std. Error"]
  pval <- co["year_death", "Pr(>|t|)"]
  tibble(
    APC = (exp(beta) - 1) * 100,
    APC_lower = (exp(beta - 1.96 * se) - 1) * 100,
    APC_upper = (exp(beta + 1.96 * se) - 1) * 100,
    p_value = ifelse(pval < 0.001, "<0.001",
                     formatC(pval, format = "f", digits = 3)
    )
  )
}

apc_segmented_one <- function(df) {
  df <- df |>
    filter(is.finite(value), value > 0) |>
    arrange(year_death)
  if (nrow(df) < 6) {
    return(tibble())
  }
  
  fmt_p <- function(p) {
    ifelse(p < 0.001, "<0.001",
           formatC(p, format = "f", digits = 3)
    )
  }
  
  # Step 1: detect breakpoint using OLS (segmented requires lm)
  fit_ols <- lm(log(value) ~ year_death, data = df)
  
  out <- tryCatch(
    withCallingHandlers(
      segmented::segmented(fit_ols, seg.Z = ~year_death, npsi = 1),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
  
  # Step 2a: no breakpoint — single Prais-Winsten on full series
  if (is.null(out) || is.null(out$psi)) {
    pw <- prais::prais_winsten(log(value) ~ year_death,
                               data = df,
                               index = "year_death"
    )
    co <- summary(pw)$coefficients
    b <- co["year_death", "Estimate"]
    s <- co["year_death", "Std. Error"]
    p <- co["year_death", "Pr(>|t|)"]
    return(tibble(
      period     = sprintf("%d-%d", min(df$year_death), max(df$year_death)),
      APC        = (exp(b) - 1) * 100,
      APC_lower  = (exp(b - 1.96 * s) - 1) * 100,
      APC_upper  = (exp(b + 1.96 * s) - 1) * 100,
      p_value    = fmt_p(p),
      breakpoint = NA_real_
    ))
  }
  
  # Step 2b: breakpoint found — split data and run Prais-Winsten on each segment
  psi <- out$psi[, "Est."]
  
  df1 <- df |> filter(year_death <= round(psi))
  df2 <- df |> filter(year_death >= round(psi))
  
  pw_seg <- function(d) {
    # need at least 4 points for Prais-Winsten to be meaningful
    if (nrow(d) < 4) {
      return(list(b = NA, s = NA, p = NA))
    }
    pw <- tryCatch(
      prais::prais_winsten(log(value) ~ year_death,
                           data = d,
                           index = "year_death"
      ),
      error = function(e) NULL
    )
    if (is.null(pw)) {
      return(list(b = NA, s = NA, p = NA))
    }
    co <- summary(pw)$coefficients
    list(
      b = co["year_death", "Estimate"],
      s = co["year_death", "Std. Error"],
      p = co["year_death", "Pr(>|t|)"]
    )
  }
  
  s1 <- pw_seg(df1)
  s2 <- pw_seg(df2)
  
  tibble(
    period = c(
      sprintf("%d-%.0f", min(df$year_death), psi),
      sprintf("%.0f-%d", psi, max(df$year_death))
    ),
    APC = c((exp(s1$b) - 1) * 100, (exp(s2$b) - 1) * 100),
    APC_lower = c(
      (exp(s1$b - 1.96 * s1$s) - 1) * 100,
      (exp(s2$b - 1.96 * s2$s) - 1) * 100
    ),
    APC_upper = c(
      (exp(s1$b + 1.96 * s1$s) - 1) * 100,
      (exp(s2$b + 1.96 * s2$s) - 1) * 100
    ),
    p_value = fmt_p(c(s1$p, s2$p)),
    breakpoint = psi
  )
}

# Plots ####
plot_comorb <- function(comorb_name) {
  d <- dsr_plot |> filter(comorbidity == comorb_name)
  ggplot(d, aes(year, rate, color = region, fill = region)) +
    geom_ribbon(aes(ymin = low, ymax = high), color = NA, alpha = RIBBON_A) +
    geom_line(aes(linewidth = (region == "Brazil"))) +
    geom_point(size = PT_SIZE, alpha = 0.9) +
    scale_linewidth_manual(values = c(
      `TRUE` = LN_SIZE + 0.1,
      `FALSE` = LN_SIZE - 0.4
    ), guide = "none") +
    scale_color_manual(values = region_colors, drop = FALSE) +
    scale_fill_manual(values = region_colors, drop = FALSE) +
    scale_x_year() +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    labs(
      title = comorb_name, x = NULL,
      y = "Age-standardized sepsis (A41) rate\nper 100,000 population"
    ) +
    theme_paper()
}

# age std rate per region 
plot_by_region <- function(name){
  dsr_plot |> filter(comorbidity == name) |> 
    ggplot(
      aes(year, rate, color = region, group = region)
    ) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, color = NA, fill = "gray70") +
    geom_line(
      linewidth = 0.5,
      alpha = 0.8
    ) +
    geom_point(
      size = PT_SIZE,
      alpha = 0.8
    ) +
    geom_point(
      size = PT_SIZE, stroke = .2, shape = 1, color = "black"
    ) +
    scale_color_manual(values = region_colors)+
    scale_x_continuous(
      breaks = c(2005, 2009, 2013, 2017, 2021, 2024),
      limits = c(2005, 2028)
    ) +
    scale_y_continuous(
      labels = label_number(accuracy = 1),
      limits = c(0, 40)
    ) +
    labs(
      x = NULL,
      color = "Region",
      y = "Age-standardized mortality rate\nper 100,000"
    ) +
    guides(color = guide_legend(override.aes = list(
      linetype = 0,
      size = 2
    ))) +
    ggrepel::geom_label_repel(
      aes(
        color = region,
        label = label_var
      ),
      size              = 2.5,        # ← smaller text
      label.padding     = unit(0.1, "lines"),  # ← tighter box
      label.size        = 0.2,        # ← thinner box border
      direction         = "y",
      xlim              = c(2025.1, NA),
      hjust             = 1,
      segment.size      = 0.5,
      segment.alpha     = .5,
      segment.linetype  = "dotted",
      min.segment.length = 0,
      box.padding       = .6,         # ← more separation between labels
      segment.curvature = -0.1,
      segment.ncp       = 3,
      segment.angle     = 20,
      max.overlaps      = Inf         # ← never drop labels
    )+
    theme_paper() +
    theme(legend.position = "none")
}


plot_underlying <- function(df, title_text, subtitle_text, ymax = 0.4) {
  # df$group <- factor(df$group,
  #                    levels = c("Underlying cause", "Elsewhere on certificate"))
  ggplot(df, aes(year_death, prop, color = group)) +
    geom_line(linewidth = LN_SIZE) +
    geom_point(size = PT_SIZE) +
    # scale_color_manual(values = c(`Underlying cause` = "#1F77B4",
    #                               `Elsewhere on certificate` = "#FF7F0E")) +
    scale_x_year() +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, ymax)) +
    labs(
      title = title_text, subtitle = subtitle_text, x = NULL,
      y = "Proportion of deaths"
    ) +
    theme_paper()
}


prop_ill_region <- function(data) {
  underlying_prop_region(data,cause =  "Ill defined") |>
    mutate(
      region = fct_recode(region,
                          "North"="1",
                          "Northeast"="2",
                          "Central-west"="5",
                          "Southeast"="3",
                          "South"="4")
    ) |> 
    mutate(label_var = if_else(year_death == 2024, region, NA_character_)) |> 
    ggplot(
      aes(year_death, prop, color = region, group = region)
    ) +
    geom_line(
      linewidth = 0.5,
      alpha = 0.8
    ) +
    geom_point(
      size = PT_SIZE,
      alpha = 0.8
    ) +
    geom_point(
      size = PT_SIZE, stroke = .2, shape = 1, color = "black"
    ) +
    scale_color_manual(values = region_colors)+
    scale_x_continuous(
      breaks = c(2005, 2009, 2013, 2017, 2021, 2024),
      limits = c(2005, 2028)
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, .5)
    ) +
    labs(
      x = NULL,
      color = "Region",
      y = "Proportion of ill-defined as\nimmediate cause of death"
    ) +
    guides(color = guide_legend(override.aes = list(
      linetype = 0,
      size = 2
    ))) +
    ggrepel::geom_label_repel(
      aes(
        color = region,
        label = label_var
      ),
      size              = 2.5,        # ← smaller text
      label.padding     = unit(0.1, "lines"),  # ← tighter box
      label.size        = 0.2,        # ← thinner box border
      direction         = "y",
      xlim              = c(2025.1, NA),
      hjust             = 1,
      segment.size      = 0.5,
      segment.alpha     = .5,
      segment.linetype  = "dotted",
      min.segment.length = 0,
      box.padding       = .6,         # ← more separation between labels
      segment.curvature = -0.1,
      segment.ncp       = 3,
      segment.angle     = 20,
      max.overlaps      = Inf         # ← never drop labels
    )+
    theme_paper() +
    theme(legend.position = "none")
  
}

# Stacks any list of plots, legend only on the last, A/B/C tags.
stack_panels <- function(plots) {
  n <- length(plots)
  plots <- lapply(seq_along(plots), function(i) {
    plots[[i]] + theme(legend.position = if (i == n) "bottom" else "none")
  })
  patchwork::wrap_plots(plots, ncol = 1) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold"))
}


underlying_prop_cause <- function(df_filter, icd_codes, cause = "A41") {
  bind_rows(df_filter |>
              mutate(group = if_else(substr(CAUSABAS, 1, 3) %in% icd_codes,
                                     "Underlying cause", "Elsewhere on certificate"
              )) |>
              count(group, year_death, immediate_cause_grouped) |>
              group_by(group, year_death) |>
              mutate(prop = n / sum(n)) |>
              ungroup() |>
              filter(immediate_cause_grouped == cause))
}
underlying_prop_all <- function(df_filter, cause = "A41") {
  bind_rows(df_filter |>
              mutate(group = "all") |>
              count(group, year_death, immediate_cause_grouped) |>
              group_by(group, year_death) |>
              mutate(prop = n / sum(n)) |>
              ungroup() |>
              filter(immediate_cause_grouped == cause))
}

underlying_prop_region <- function(df_filter, cause = "A41") {
  bind_rows(df_filter |>
              count(year_death, region, immediate_cause_grouped) |>
              group_by(year_death, region) |>
              mutate(prop = n / sum(n)) |>
              ungroup() |>
              filter(immediate_cause_grouped == cause))
}


scale_x_year <- function()
  scale_x_continuous(breaks = c(2005, 2008, 2011, 2014, 2017, 2020, 2024), limits = c(2005, 2024))
#Theme
theme_paper <- function() {
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "Roboto"),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.major.y = element_line(color = "grey92"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey30", size = 10),
      legend.position = "bottom",
      legend.key.width = unit(18, "pt"),
      legend.text = element_text(size = 10),
      axis.title.y = element_text(size = 9, color = "grey20"),
      axis.text = element_text(color = "grey20")
    )
}

