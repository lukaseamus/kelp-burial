# 1. Data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

# Core meta-data
cores <- here("Carbon", "Cores.csv") %>%
  read_csv() %T>%
  print()

# Sediment dry mass
mass <- here("Carbon", "Mass.csv") %>%
  read_csv() %T>%
  print()

# Organic mass
acid <- here("Carbon", "Acidification.csv") %>%
  read_csv() %>%
  mutate(
    Residual = (Final - Tube) / Initial, # Proportion of mass that is not CaCO3
    CaCO3 = 1 - Residual # Proportion of mass that is CaCO3
  ) %>%
  select(ID, Residual, CaCO3) %T>%
  print()

# Sediment carbon
sed_carb <- here("Carbon", "d13C_Sediment.csv") %>%
  read_csv() %>%
  filter(Reliable) %T>% # Take only reliable samples
  print()

sed_carb %>%
  group_by(ID) %>%
  count() %>%
  filter(n > 1) %>%
  print(n = 20)
# 18 reliable samples have duplicated IDs 

# Check if they are true replicates that may be averaged
sed_carb %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  mutate(Replicate = row_number()) %>%
  select(ID, Replicate, Mass) %>%
  pivot_wider(names_from = Replicate, values_from = Mass) %>%
  mutate(
    Unique = n_distinct(c_across(1:3), na.rm = TRUE) ==
      sum(!is.na(c_across(1:3)))
  )
# All duplicates are true, unique replicates

# I'm taking the mean of samples mass, carbon content and d13C,
# and picking the first 
sed_carb %<>%
  group_by(ID) %>%
  summarise(
    across(where(is.numeric), ~mean(.x)),
    across(c(where(is.character), where(is.logical)), first),
    n = n() # Save number of reps
  ) %T>%
  print()

sed_carb %>%
  group_by(ID) %>%
  count() %>%
  filter(n > 1)
# No more duplicates

# Kelp carbon
kelp_carb <- here("Carbon", "d13C_Kelp.csv") %>%
  read_csv() %T>%
  print()

# 1.2 Combine sediment data ####
sed <- mass %>%
  full_join(acid) %>%
  full_join(
    sed_carb %>% 
      rename(Sample = Mass) %>%
      mutate(
        Acidified = if_else(
          ID %>% str_detect("inorganic"),
          FALSE, TRUE
        ),
        ID = ID %>% 
          str_remove("_inorganic")
      )
  ) %>%
  mutate(ID_core = ID %>% str_extract("^[A-Za-z]+(?:_[A-Za-z]+)?_\\d+_\\d+")) %>%
  full_join(cores %>% rename(ID_core = ID)) %T>%
  print()

# 1.3 Extract and calculate variables ####
# Experimental units and duration are contained in ID and sampling dates.
# The variable C is %C of the acidified sample, but this is not %Corg
# of sediment dry mass. For this I need to account for CaCO3 mass lost
# during acidification. To do so, I convert C from % to proportion and 
# multiply by the porportion of dry mass left after acidification.
# CaCO3 can be converted to Cinorg as CaCO3 / 100.0869 * 12.0107.

sed %<>%
  mutate(
    Experiment = if_else(
      ID %>% str_detect("L"),
      "United Kingdom", "Australia"
    ),
    Tank = if_else(
      ID %>% str_count("_") == 3,
      ID %>% str_split_i("_", 2),
      ID %>% str_split_i("_", 3)
    ) %>% as.numeric(),
    Core = if_else(
      ID %>% str_count("_") == 3,
      ID %>% str_split_i("_", 3),
      ID %>% str_split_i("_", 4)
    ) %>% as.numeric(),
    Depth = if_else(
      ID %>% str_count("_") == 3,
      ID %>% str_split_i("_", 4),
      ID %>% str_split_i("_", 5)
    ) %>% as.numeric() - 0.5,
    Corg = C / 100 * Residual,
    Cinorg = CaCO3 / 100.0869 * 12.0107,
    Added = Added %>% dmy(),
    Sampled = Sampled %>% dmy(),
    Day = Added %--% Sampled / ddays(),
    Day = if_else(Species == "Control", 0, Day)
  ) %>%
  arrange(Experiment, Species, Tank, Core, Depth) %T>%
  print()

# 1.4 Calculate source material summaries ####
# 1.4.1 Baseline sediment ####
baseline_summary <- sed %>%
  filter(Species == "Control") %>%
  select(Sediment, Cinorg, Corg, d13C) %>%
  mutate(Cinorg = Cinorg * 100, # Convert to percentage
         Corg = Corg * 100) %>%
  group_by(Sediment) %>%
  summarise(
    across(
      everything(),
      list(
        mean = ~ .x %>% mean(na.rm = T),
        sd = ~ .x %>% sd(na.rm = T),
        median = ~ .x %>% median(na.rm = T),
        n = ~ sum(!is.na(.x))
      )
    )
  ) %T>%
  print()
# Mean and median are practically identical in all cases.
# Summarise for Table 1:
require(glue)
baseline_summary %<>%
  mutate(
    # Round
    across(
      where(is.numeric),
      ~ .x %>% signif(2)
    ),
    # Glue
    Cinorg = glue("{Cinorg_mean} ± {Cinorg_sd} ({Cinorg_n})"),
    Corg = glue("{Corg_mean} ± {Corg_sd} ({Corg_n})"),
    d13C = glue("{d13C_mean} ± {d13C_sd} ({d13C_n})")
  ) %>%
  select(Sediment, Cinorg, Corg, d13C) %T>%
  print()

baseline_summary %>%
  write_csv(here("Tables", "Table_1b.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = baseline_summary) %>%
  print(target = here("Tables", "Table_1b.docx"))

# 1.4.2 Enriched kelp ####
kelp_summary <- kelp_carb %>%
  filter(Date == "26.09.25" & Species != "Saccharina latissima") %>%
  select(Species, C, d13C) %>%
  group_by(Species) %>%
  summarise(
    across(
      everything(),
      list(mean = mean, sd = sd, median = median)
    ),
    n = n()
  ) %T>%
  print()
# Mean and median are practically identical in all cases.
# Summarise for Table 2:
kelp_summary %<>%
  mutate(
    # Round
    across(
      where(is.numeric),
      ~ case_when(
          .x < 100 ~ .x %>% signif(2),
          .x < 1e3 ~ .x %>% signif(3),
          T ~ .x %>% signif(4)
        )
    ),
    # Glue
    C = glue("{C_mean} ± {C_sd}"),
    d13C = glue("{d13C_mean} ± {d13C_sd}")
  ) %>%
  select(Species, C, d13C) %T>%
  print()

kelp_summary %>%
  write_csv(here("Tables", "Table_2.csv"))

read_docx() %>%
  body_add_table(value = kelp_summary) %>%
  print(target = here("Tables", "Table_2.docx"))

# 1.5 Filter data for model ####
# # Filter only cores that have been processed
# sed %<>%
#   filter(
#     (Experiment == "United Kingdom" | Depth < 10.5) &
#       !is.na(Organic)
#   ) %T>%
#   print()
  
sed %<>%
  filter(
    (Experiment == "United Kingdom" | Depth < 10.5) &
      !is.na(d13C) & Acidified
  ) %T>%
  print()

# 1.6 Visualise data ####
sed %>%
  ggplot() +
    geom_point(aes(Corg * 100, Depth, colour = Sediment), alpha = 0.2, shape = 16) +
      scale_y_reverse() +
      facet_grid(~ Species) +
      theme_minimal()

sed %>% 
  filter(Species != "Control") %>%
  ggplot() +
    geom_point(aes(Day, Corg * 100, colour = Sediment), alpha = 0.2, shape = 16) +
      facet_grid(Depth ~ Species, scales = "free_x") +
      theme_minimal()

sed %>%
  ggplot() +
    geom_point(aes(d13C, Depth, colour = Sediment), alpha = 0.2, shape = 16) +
      scale_y_reverse() +
      facet_grid(~ Species) +
      theme_minimal()

sed %>% 
  filter(Species != "Control") %>%
  ggplot() +
    geom_point(aes(Day, d13C, colour = Sediment), alpha = 0.2, shape = 16) +
      facet_grid(Depth ~ Species, scales = "free_x") +
      theme_minimal()

sed %>%
  ggplot() +
    geom_point(aes(CaCO3 * 100, Depth, colour = Sediment), alpha = 0.2, shape = 16) +
      scale_y_reverse() +
      facet_grid(~ Species) +
      theme_minimal()

sed %>% 
  filter(Species != "Control") %>%
  ggplot() +
    geom_point(aes(Day, CaCO3 * 100, colour = Sediment), alpha = 0.2, shape = 16) +
      facet_grid(Depth ~ Species, scales = "free_x") +
      theme_minimal()


