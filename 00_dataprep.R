# Load required packages
library(tidytable)
library(stringr)
# usr fun#### faster calculation for age
age_fun <- function(from, to){
  from_lt <- as.POSIXlt(from)
  to_lt <- as.POSIXlt(to)
  age = to_lt$year - from_lt$year
  ifelse(to_lt$mon < from_lt$mon |
           (to_lt$mon == from_lt$mon & to_lt$mday < from_lt$mday),
         age -1, age)
}

# Define comorbidity codes and labels
comorbids <- c("E10", "E11", "E14", "I10", "N18")
dados_raw  <- arrow::read_parquet("sim0524.parquet")

# Convert age ####
dados_raw <- dados_raw |> 
  mutate(
    DTNASC = lubridate::dmy(DTNASC),
    DTOBITO = lubridate::dmy(DTOBITO),
    IDADE = as.numeric(IDADE),
    age_years_form = case_when(
      IDADE<400 ~ 0,
      IDADE>=400 & IDADE<=599 ~ IDADE-400
    ),
    age_calc_nasc = age_fun(DTNASC, DTOBITO),
    age = case_when(
      is.na(age_calc_nasc) ~ age_years_form,
      age_years_form!=age_calc_nasc ~ age_years_form,
      TRUE ~ age_calc_nasc
    )) 

# icd <- readxl::read_excel("icd_clean.xlsx")
# icd <- icd |> select(
#   ICD10_3_Code,
#   ICD10_Code,
#   Chapter_No,
#   Group_Code,
#   Group_Desc
# )
# icd <- icd |>  filter(!duplicated(ICD10_Code))
# icd <- icd |>  mutate(ICD10_Code = stringr::str_remove(ICD10_Code, "\\."))
dados  <- dados_raw  |> 
  mutate(across(starts_with("LINHA"), ~ stringr::str_remove(., "^\\*")))|> 
  mutate(
    LINHAA = if_else(
      stringr::str_detect(LINHAA,"\\*"),
      strex::str_before_first(LINHAA,"\\*"),
      LINHAA
    ),
    LINHAB = if_else(
      stringr::str_detect(LINHAB,"\\*"),
      strex::str_before_first(LINHAB,"\\*"),
      LINHAB
    ),
    LINHAC = if_else(
      stringr::str_detect(LINHAC,"\\*"),
      strex::str_before_first(LINHAC,"\\*"),
      LINHAC
    ),
    LINHAD = if_else(
      stringr::str_detect(LINHAD,"\\*"),
      strex::str_before_first(LINHAD,"\\*"),
      LINHAD
    ),# keep only the first ICD
    LINHAA = na_if(LINHAA,""),
    LINHAB = na_if(LINHAB,""),
    LINHAC = na_if(LINHAC,""),
    LINHAD = na_if(LINHAD,""),
    # LINHAB = case_when(
    #   LINHAB == LINHAA ~ NA, 
    #   TRUE ~ LINHAB),
    # LINHAC = case_when(
    #   LINHAC == LINHAB ~ NA,
    #   TRUE ~ LINHAC),
    # LINHAD = case_when(
    #   LINHAD == LINHAC ~ NA,
    #   TRUE ~ LINHAD), #fix
    # LINHAA = coalesce(LINHAA,LINHAB,LINHAC, LINHAD),
    # LINHAB = if_else(LINHAB == LINHAA, NA, LINHAB),
    # LINHAB = coalesce(LINHAB,LINHAC, LINHAD),
    # LINHAC = if_else(LINHAC == LINHAB, NA, LINHAC),
    # LINHAC = coalesce(LINHAC, LINHAD),
    # LINHAD = if_else(LINHAD == LINHAC, NA, LINHAD), #fix    
    # LINHAA = coalesce(LINHAA,LINHAB,LINHAC, LINHAD),
    # LINHAB = if_else(LINHAB == LINHAA, NA, LINHAB),
    # LINHAB = coalesce(LINHAB,LINHAC, LINHAD),
    # LINHAC = if_else(LINHAC == LINHAB, NA, LINHAC),
    # LINHAC = coalesce(LINHAC, LINHAD),
    # LINHAD = if_else(LINHAD == LINHAC, NA, LINHAD), #fix
    linha_a_ill = case_when(
      stringr::str_detect(LINHAA,"R95") ~ "R95",
      stringr::str_detect(LINHAA,"^R") ~ "Ill defined",
      TRUE ~ LINHAA
    ),
    linha_a_ill = case_when(
      linha_a_ill =="Ill defined" ~ "Ill defined",
      TRUE ~ stringr::str_sub(linha_a_ill,1,3)
    )
  )

# Indicador de cada comorbidade

dados  <- dados  |> 
  mutate(
    e10 = case_when(
      str_detect(LINHAB,"E10") |
        str_detect(LINHAC,"E10")|
        str_detect(LINHAD,"E10") |
        str_detect(CAUSABAS,"E10")|
        str_detect(LINHAII,"E10")~ 1,
      TRUE ~ 0
    ),
    e11 = case_when(
      str_detect(LINHAB,"E11") |
        str_detect(LINHAC,"E11")|
        str_detect(LINHAD,"E11") |
        str_detect(CAUSABAS,"E11")|
        str_detect(LINHAII,"E11")~ 1,
      TRUE ~ 0
    ),
    e14 = case_when(
      str_detect(LINHAB,"E14") | 
        str_detect(LINHAC,"E14")|
        str_detect(LINHAD,"E14") |
        str_detect(CAUSABAS,"E14")|
        str_detect(LINHAII,"E14")~ 1,
      TRUE ~ 0
    ),
    has = case_when(
      str_detect(LINHAB,"I10") |
        str_detect(LINHAC,"I10")|
        str_detect(LINHAD,"I10") |
        str_detect(CAUSABAS,"I10")|
        str_detect(LINHAII,"I10")~ 1,
      TRUE ~ 0
    ),
    ckd = case_when(
      str_detect(LINHAB,"N18") |
        str_detect(LINHAC,"N18")|
        str_detect(LINHAD,"N18") |
        str_detect(CAUSABAS,"N18")|
        str_detect(LINHAII,"N18")~ 1,
      TRUE ~ 0
    ),
    dm = case_when(
      rowSums(across(c(e10:e14)))>=1 ~ 1,
      TRUE ~ 0
    ),
    number_commorbidities = rowSums(across(c(has,ckd,dm))),
    e10_linea = case_when(
      str_detect(LINHAA,"E10")~ 1,
      TRUE ~ 0
    ),
    e11_linea = case_when(
      str_detect(LINHAA,"E11")~ 1,
      TRUE ~ 0
    ),
    e14_linea = case_when(
      str_detect(LINHAA,"E14") ~ 1,
      TRUE ~ 0
    ),
    has_linea = case_when(
      str_detect(LINHAA,"I10") ~ 1,
      TRUE ~ 0
    ),
    ckd_linea = case_when(
      str_detect(LINHAA,"N18") ~ 1,
      TRUE ~ 0
    ),
    dm_linea = case_when(
      rowSums(across(c(e10_linea:e14_linea)))>=1 ~ 1,
      TRUE ~ 0
    ),
    number_commorbidities_linea = rowSums(across(c(has_linea,ckd_linea,dm_linea))),
  )



# line_a_comp <- dados |>
#   mutate(
#     year_death = lubridate::year(as.Date(DTOBITO)),
#     cat = case_when(
#     substr(LINHAA, 1, 3) == "A41"            ~ "Sepsis",
#     substr(LINHAA, 1, 3) %in% c("J96") ~ "Resp failure",
#     substr(LINHAA, 1, 1) == "R"   & LINHAA!="R95"           ~ "Ill-defined",
#     .default = "Other specified"
#   )) |>
#   count(year_death, cat) |>
#   mutate(prop = n / sum(n), .by = c(year_death))
# 
# ggplot2::ggplot(line_a_comp, ggplot2::aes(year_death, prop, fill = cat)) +
#   ggplot2::geom_area() 


# Filter possible errors #####

dados <- dados |> 
  tidylog::filter(number_commorbidities>=1) |> 
  tidylog::filter(!is.na(linha_a_ill)) |> 
  tidylog::filter(age<110, age>=0) |> 
  tidylog::filter(SEXO!=0,
                  SEXO!=9,
                  !is.na(SEXO))

# filter: removed 19,295,035 rows (75%), 6,382,022 rows remaining
# filter: removed 252,601 rows (4%), 6,129,421 rows remaining
# filter: removed 2,163 rows (<1%), 6,127,258 rows remaining
# filter: removed 428 rows (<1%), 6,126,830 rows remaining

dados <- dados  |> 
  mutate(
    sex = factor(SEXO, levels = c(1, 2), labels = c("Male", "Female")),
    year_death = lubridate::year(as.Date(DTOBITO)),
    region = str_sub(CODMUNRES,1,1),
    age_group = if_else(age>=60,"≥60","<60"),
    race_ethnicity = case_when(
      RACACOR == "1" ~ "White",
      RACACOR == "2" ~ "Black",
      RACACOR == "3" ~ "Yellow",
      RACACOR == "4" ~ "Mixed",
      RACACOR == "5" ~ "Indigenous",
      TRUE ~ NA_character_
    )
  )


dados_comorb <- dados |>
  tidylog::filter(number_commorbidities >= 1) |>
  mutate(
    top_codes = list(levels(forcats::fct_lump(linha_a_ill, n = 11))),
    immediate_cause_grouped = forcats::fct_other(
      linha_a_ill,
      keep = union(setdiff(top_codes[[1]], "Other"), comorbids)
    )
  ) |>
  dplyr::select(-top_codes)

arrow::write_parquet(dados_comorb,"dados_comorb0524_151126.parquet")
