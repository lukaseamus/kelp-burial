# 1. Data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

# UK
coulter <- here("Granulometry", "Coulter.csv") %>%
  read_csv() %T>%
  print()

# Australia
master <- here("Granulometry", "Mastersizer.csv") %>%
  read_csv() %T>%
  print()

# UK and Australian particle size data are from different instruments,
# so are output differently.

# 1.2 Wrangle bins ####

coulter %<>%
  mutate(
    Lower = Bin,
    Upper = lead(Bin), # These data are leaded
    Mid = (Lower + Upper)/2
  ) %>%
  drop_na(Percentage) %T>%
  print()


master %<>% # Use only technical averages
  filter(Sample %>% str_detect("Average")) %>%
  pivot_longer(
    cols = -c(Sample, Datetime, Mode, 
              starts_with("Mean"), starts_with("SD"),
              Skewness, Kurtosis),
    names_to = "Bin",
    values_to = "Percentage"
  ) %>%
  mutate(
    Bin = Bin %>% as.numeric(),
    Lower = lag(Bin), # These data are lagged
    Upper = Bin,
    Mid = (Lower + Upper)/2
  ) %>%
  drop_na(Percentage) %T>%
  print()

# 1.3 Combine data ####
granulo <- bind_rows(
  UK = coulter %>%
    mutate(Sample = Sample %>% as.character()) %>%
    select(Sample, Percentage, Lower, Upper, Mid),
  Australia = master %>%
    select(Sample, Percentage, Lower, Upper, Mid),
  .id = "Location"
) %>%
  mutate(
    Sediment = case_when(
      Location == "UK" ~ "Silt",
      Location == "Australia" &
        Sample %>% str_detect("1_1") ~ "Coarse sand",
      Location == "Australia" &
        Sample %>% str_detect("1_2") ~ "Fine sand"
    )
  ) %T>%
  print()

# 1.4 Simulate particles ####
# There are 6 decimal places in some percentages,
# but simulating enough particles is not possible,
# so rounding is necessary.
granulo %<>% 
  mutate(Particles = round(Percentage / 100 * 1e6)) %>%
  rowwise() %>%
  mutate(Size = list( runif(Particles, Lower, Upper) )) %>%
  ungroup() %>%
  unnest(Size) %T>%
  print()

granulo %>%
  ggplot() +
    geom_density(aes(Size), colour = NA, fill = "black") +
    facet_wrap(~ Sample) +
    theme_minimal()
# Very right-skewed

granulo %>%
  ggplot() +
    geom_density(aes(Size), colour = NA, fill = "black") +
    scale_x_log10() +
    facet_wrap(~ Sample) +
    theme_minimal()
# Approximately lognormal


# 1.5 Calculate summaries by sediment type ####
size_summary <- granulo %>%
  mutate(log_Size = log10(Size),
         Phi = -log(Size/1e3, base = 2)) %>%
  group_by(Sediment) %>%
  summarise(
    across(
      c(Size, log_Size, Phi),
      list(mean = mean, sd = sd, median = median)
    ),
    n = n_distinct(Sample)
  ) %T>%
  print()
# Log mean and median are close, so the lognormal
# distribution is good.

require(glue)
size_summary %<>%
  mutate(
    across(
      where(is.numeric),
      ~ .x %>% signif(2)
    ),
    log_Size = glue("{log_Size_mean} ± {log_Size_sd} ({n})"),
    Phi = glue("{Phi_mean} ± {Phi_sd} ({n})")
  ) %>%
  select(Sediment, log_Size, Phi) %T>%
  print()

size_summary %>%
  write_csv(here("Tables", "Table_1a.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = size_summary) %>%
  print(target = here("Tables", "Table_1a.docx"))
