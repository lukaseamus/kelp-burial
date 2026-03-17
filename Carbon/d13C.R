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

# 1.6 Assign control to treatments at t0 ####
sed %<>%
  mutate(
    # Assign controls to species
    Species = case_when(
      Species == "Control" & Experiment == "Australia" ~
        "Ecklonia radiata",
      ID_core %in% c("L_1_1", "L_2_1", "L_2_3") ~ "Laminaria hyperborea",
      ID_core %in% c("L_1_2", "L_2_2", "L_2_7") ~ "Laminaria ochroleuca",
      TRUE ~ Species
    ),
    # Round kelp load
    Load = case_when(
      Kelp == 0 ~ 0,
      Kelp < 2 ~ 1,
      Kelp < 8 ~ 7.5,
      Kelp < 11 ~ 10,
      Kelp > 11 ~ 15
    ),
    # Assign controls to loads
    Load = case_when(
      Load == 0 & Experiment == "Australia" ~ 10,
      ID_core %in% c("L_1_1", "L_1_2") ~ 1,
      ID_core %in% c("L_2_1", "L_2_2") ~ 7.5,
      ID_core %in% c("L_2_3", "L_2_7") ~ 15,
      TRUE ~ Load
    )
  ) %T>%
  print()

# 1.6 Calculate fraction of kelp ####
sed %<>%
  full_join( # Sediment source material
    sed %>% 
      filter(Day == 0) %>%
      summarise(d13C_sed = mean(d13C),
                .by = c(Sediment, Depth))
  ) %>%
  full_join( # Kelp source material
    kelp_carb %>%
      filter(Date == "26.09.25" & Species != "Saccharina latissima") %>%
      summarise(Corg_kelp = mean(C/100),
                d13C_kelp = mean(d13C),
                .by = Species)
  ) %>%
  full_join( # Empty bag mass
    mass %>%
      filter(ID %>% str_detect("blank") & !is.na(Mass)) %>%
      summarise(Blank = mean(Mass), .by = ID) %>%
      mutate(
        Experiment = if_else(
          ID %>% str_detect("L"),
          "United Kingdom", "Australia"
        )
      ) %>% 
      select(-ID)
  ) %>%
  mutate(
    # Calculate kelp fraction of sediment organic carbon (mixing model)
    Fraction = ( d13C - d13C_sed ) / ( d13C_kelp - d13C_sed ),
    # Correct total sediment mass for empty bag mass
    Mass_c = Mass - Blank,
    # Infer kelp fraction at t0
    Fraction = case_when(
      Day == 0 & Depth == 0.5 ~ # All kelp is initially on surface
        Load * Corg_kelp / (Load * Corg_kelp + Mass_c * Corg),
      Day == 0 & Depth > 1 ~ 0, # No kelp below
      TRUE ~ Fraction
    )
  ) %T>%
  print()

sed %$% range(Fraction)
# Some negative fractions, meaning the sediment baseline d13C is
# greater than the sediment sample d13C. Negative values can be 
# replaced with zero.

sed %<>%
  mutate(
    Fraction = if_else( Fraction < 0 , 0 , Fraction )
  ) %T>%
  print()

sed %$% range(Fraction) # Fixed

# 1.7 Calculate remaining kelp carbon ####
sed %<>%
  mutate(
    # Calculate added kelp carbon
    C_kelp_initial = Kelp * Corg_kelp, 
    # Calculate remaining kelp carbon
    C_kelp_final = Mass_c * Corg * Fraction, 
    # Calculate remaining proportion
    Proportion = C_kelp_final / C_kelp_initial,
    # Proportion is 1 at t0 on surface, 0 below
    Proportion = case_when(
      Day == 0 & Depth == 0.5 ~ 1,
      Day == 0 & Depth > 1 ~ 0,
      TRUE ~ Proportion
    )
  ) %T>%
  print()

sed %$% range(Proportion) # All good

# 1.8 Sum across horizons per core ####
sed_sum <- sed %>%
  summarise(
    across(c(Mass, Mass_c, Proportion), sum),
    .by = c(ID_core, Species, Sediment, Load, Day)
  ) %T>%
  print(n = 50)

sed_sum %<>% # Convert categorical variables to factors
  mutate(Species = Species %>% fct(), 
         Sediment = Sediment %>% fct())

# 1.9 Calculate final kelp carbon stock ####
final <- sed %>%
  filter(Day == max(Day),
         .by = c(Species, Sediment, Load)) %>%
  # Cores were 10 cm in diameter, so radius is 0.05 m
  mutate(C_m2 = C_kelp_final / ( pi * 0.05^2 )) %T>%
  print()

final %<>% # Convert categorical variables to factors
  mutate(Species = Species %>% fct(), 
         Sediment = Sediment %>% fct())

# 2. Models ####
# 2.1 Time ####
# 2.1.1 Prior simulation ####
require(extraDistr)
tibble(n = 1:1e3,
       alpha_mu = rnorm( 1e3 , log(0.1) , 1 ),
       alpha_sigma_sp = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       alpha_sigma_sed = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       beta = rnorm( 1e3 , 0 , 0.1 ),
       gamma_mu = rnorm( 1e3 , qlogis(0.1) , 1 ),
       gamma_sigma_sp = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       gamma_sigma_sed = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       delta = rnorm( 1e3 , 0 , 0.1 ),
       nu = rgamma( 1e3 , 30^2 / 20^2 , 30 / 20^2 )) %>%
  mutate(
    k = exp(
      alpha_mu +
        rnorm( n() , 0 , alpha_sigma_sp ) + 
        rnorm( n() , 0 , alpha_sigma_sed ) + 
        beta * 10 # for 10 g kelp
    ),
    r = plogis(
      gamma_mu +
        rnorm( n() , 0 , gamma_sigma_sp ) + 
        rnorm( n() , 0 , gamma_sigma_sed ) + 
        delta * 10 # for 10 g kelp
    )
  ) %>%
  expand_grid(Day = sed_sum %$% 
                seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    mu = (1 - r) * exp( -k * Day ) + r,
    Proportion = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %>%
  pivot_longer(cols = c(mu, Proportion),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.1.2 Stan model ####
require(cmdstanr)
time_model <- here("Carbon", "Stan", "time.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
time_samples <- time_model$sample(
          data = sed_sum %>%
            filter(Day != 0) %>% # remove t0
            mutate(Load = Load - 10) %>% # centre load
            select(Day, Proportion, Species, 
                   Sediment, Load) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

# Save draws
time_samples$draws() %>%
  write_rds(here("Carbon", "RDS", "time_samples.rds"))
time_samples$draws(format = "df") %>%
  write_rds(here("Carbon", "RDS", "time_samples_df.rds"))

# 2.1.3 Model checks ####
# Rhat
time_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000145.

# Chains
require(bayesplot)
time_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

# Pairs
time_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "gamma_mu", "beta", "delta"))
# No correlation

# 2.1.4 Prior-posterior comparison ####
source("functions.R")
time_prior <- prior_samples(
  model = time_model,
  data = sed_sum %>%
    filter(Day != 0) %>%
    mutate(Load = Load - 10) %>%
    select(Day, Proportion, Species, 
           Sediment, Load) %>%
    compose_data()
)

time_prior %>% 
  prior_posterior_draws(
    posterior_samples = time_samples,
    group = sed_sum %>% select(Species, Sediment),
    parameters = c("alpha_mu", "alpha_sigma_sp", "alpha_sp[Species]",
                   "alpha_sigma_sed", "alpha_sed[Sediment]", 
                   "beta", "gamma_mu", "gamma_sigma_sp", "gamma_sp[Species]",
                   "gamma_sigma_sed", "gamma_sed[Sediment]", 
                   "delta", "nu"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Sediment"
  )
# No real difference between species and sediments. No change with load.

# 2.1.5 Parameters ####
time_prior_posterior_global <- time_prior %>% 
  prior_posterior_draws(
    posterior_samples = time_samples,
    parameters = c("alpha_mu", "alpha_sigma_sp", "alpha_sigma_sed", "beta", 
                   "gamma_mu", "gamma_sigma_sp", "gamma_sigma_sed", "delta", 
                   "nu"),
    format = "short"
  ) %>% 
  mutate( # Calculate k and r for new kelp species and sediment types
    k = exp(
      alpha_mu + 
        rnorm( n() , 0 , alpha_sigma_sp ) +
        rnorm( n() , 0 , alpha_sigma_sed ) +
        beta * 0 # for 10 g kelp because load is centred on 10
    ), 
    t0.5 = log(2) / k, # half-life
    r = plogis(
      gamma_mu + 
        rnorm( n() , 0 , gamma_sigma_sp ) +
        rnorm( n() , 0 , gamma_sigma_sed ) +
        delta * 0
    )
  ) %T>%
  print()

time_prior_posterior <- sed_sum %>% # Get treatment combinations
  distinct(Species, Sediment, Load) %>%
  cross_join( # Join global parameters
    time_prior %>% 
      prior_posterior_draws(
        posterior_samples = time_samples,
        parameters = c("alpha_mu", "beta", "gamma_mu", "delta", "nu"),
        format = "short"
      )
  ) %>%
  full_join( # Join species parameters
    time_prior %>% 
      prior_posterior_draws(
        posterior_samples = time_samples,
        group = sed_sum %>% select(Species),
        parameters = c("alpha_sp[Species]", "gamma_sp[Species]"),
        format = "short"
      )
  ) %>%
  full_join( # Join sediment parameters
    time_prior %>% 
      prior_posterior_draws(
        posterior_samples = time_samples,
        group = sed_sum %>% select(Sediment),
        parameters = c("alpha_sed[Sediment]", "gamma_sed[Sediment]"),
        format = "short"
      )
  ) %>%
  mutate( # Calculate k and r for existing treatments
    Load_c = Load - 10, # centred load
    k = exp( alpha_mu + alpha_sp + alpha_sed + beta * Load_c ),
    t0.5 = log(2) / k, # half-life
    r = plogis( gamma_mu + gamma_sp + gamma_sed + delta * Load_c )
  ) %>% # Embed prior in treatments
  filter(Species == "Ecklonia radiata" & Sediment == "Coarse sand" &
           Load == 10 & distribution == "prior" |
             distribution == "posterior") %>%
  mutate(
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct(),
    Sediment = if_else(
      distribution == "prior", "Prior", Sediment
    ) %>% fct(),
    Load = if_else(
      distribution == "prior", NA, Load
    )
  ) %>%
  select(starts_with("."), Species, Sediment, Load, k, t0.5, r, nu) %T>%
  print()

# Save parameter distributions
time_prior_posterior_global %>%
  write_rds(here("Carbon", "RDS", "time_prior_posterior_global.rds"))
time_prior_posterior %>%
  write_rds(here("Carbon", "RDS", "time_prior_posterior.rds"))

# 2.1.6 Prediction ####
time_prediction <- time_prior_posterior %>%
  mutate(Treatment = str_c(Species, Sediment, Load, sep = " ")) %>%
  spread_continuous(
    data = sed_sum %>%
      mutate(Treatment = str_c(Species, Sediment, Load, sep = " ")),
    predictor_name = "Day",
    group_name = "Treatment"
  ) %>%
  mutate( 
    mu = (1 - r) * exp( -k * Day ) + r,
    Proportion = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %>%
  group_by(Species, Sediment, Load, Day) %>%
  median_qi(mu, Proportion, .width = c(.5, .8, .9)) %T>%
  print()

time_prediction %>%
  write_rds(here("Carbon", "RDS", "time_prediction.rds"))

# 2.2 Depth ####
# I am envisioning a model of the form L * exp( -k * Day ) + R which can be used
# to test if there is a putative refractory fraction remaining at a certain depth.


# 3. Figures ####
# Define custom theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing.x = unit(0.6, "cm"),
                 panel.spacing.y = unit(0.6, "cm"),
                 text = element_text(family = "Futura"))

# 2.1 Figure 1 ####
require(ggh4x)
Fig_1 <- sed %>%
  filter(Depth < 10) %>%
  ggplot() +
    geom_point(aes(Day, Fraction * 100, colour = Species), 
               shape = 16, size = 2.5) +
    geom_hline(yintercept = 0, lineend = "square") +
    scale_colour_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                        guide = "none") +
    facet_nested(
      Depth ~ Species + Sediment + Load, 
      nest_line = TRUE, solo_line = TRUE,
      # axes = "x", remove_labels = "all", switch = "y", 
      scales = "free", space = "free_x",
      strip = strip_nested(
        clip = "off", by_layer_x = TRUE,
        text_x = list(element_text(face = "italic"), NULL, NULL),
        text_y = element_text(angle = 0, hjust = 0, vjust = 1)
      ),
      labeller = labeller(
        Depth = c(
          "0.5" = "0–1 cm",
          "1.5" = "1–2 cm",
          "2.5" = "2–3 cm",
          "3.5" = "3–4 cm",
          "4.5" = "4–5 cm",
          "5.5" = "5–6 cm",
          "6.5" = "6–7 cm",
          "7.5" = "7–8 cm",
          "8.5" = "8–9 cm",
          "9.5" = "9–10 cm"
        ),
        Load = c(
          "1" = "1 g",
          "7.5" = "7.5 g",
          "10" = "10 g",
          "15" = "15 g"
        )
      )
    ) +
    facetted_pos_scales(
      x = list(
        Species == "Ecklonia radiata" ~ 
          scale_x_continuous(limits = c(0, 120),
                             breaks = seq(0, 120, 30)),
        Species %>% str_detect("Laminaria") ~
          scale_x_continuous(limits = c(0, 30),
                             breaks = c(0, 30))
      ),
      y = list(
        Depth > 3 ~ 
          scale_y_continuous(limits = c(0, 5),
                             breaks = seq(0, 5, 2.5),
                             labels = scales::label_number(accuracy = c(1, 0.1, 1))),
        Depth == 2.5 ~
          scale_y_continuous(limits = c(0, 8),
                             breaks = seq(0, 8, 4)),
        Depth == 1.5 ~
          scale_y_continuous(limits = c(0, 30),
                             breaks = seq(0, 30, 15)),
        Depth == 0.5 ~
          scale_y_continuous(limits = c(0, 100),
                             breaks = seq(0, 100, 50))
      )
    ) +
    labs(x = "Time since kelp addition (days)",
         y = expression("Kelp fraction of sediment C"["org"]*" (%)")) +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 0.1, 0.2, 0.2, unit = "cm"))

Fig_1

Fig_1 %>%
  ggsave(filename = "Fig_1.pdf", path = "Figures",
         device = cairo_pdf, height = 24, width = 24, units = "cm")

# 2.2 Figure 2 ####
Fig_2 <- sed %>%
  filter(Depth < 10) %>%
  ggplot() +
    geom_point(aes(Day, Proportion * 100, colour = Species), 
               shape = 16, size = 2.5) +
    geom_hline(yintercept = 0, lineend = "square") +
    scale_colour_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                        guide = "none") +
    facet_nested(
      Depth ~ Species + Sediment + Load, 
      nest_line = TRUE, solo_line = TRUE,
      # axes = "x", remove_labels = "all", switch = "y", 
      scales = "free", space = "free_x",
      strip = strip_nested(
        clip = "off", by_layer_x = TRUE,
        text_x = list(element_text(face = "italic"), NULL, NULL),
        text_y = element_text(angle = 0, hjust = 0, vjust = 1)
      ),
      labeller = labeller(
        Depth = c(
          "0.5" = "0–1 cm",
          "1.5" = "1–2 cm",
          "2.5" = "2–3 cm",
          "3.5" = "3–4 cm",
          "4.5" = "4–5 cm",
          "5.5" = "5–6 cm",
          "6.5" = "6–7 cm",
          "7.5" = "7–8 cm",
          "8.5" = "8–9 cm",
          "9.5" = "9–10 cm"
        ),
        Load = c(
          "1" = "1 g",
          "7.5" = "7.5 g",
          "10" = "10 g",
          "15" = "15 g"
        )
      )
    ) +
    facetted_pos_scales(
      x = list(
        Species == "Ecklonia radiata" ~
          scale_x_continuous(limits = c(0, 120),
                             breaks = seq(0, 120, 30)),
        Species %>% str_detect("Laminaria") ~
          scale_x_continuous(limits = c(0, 30),
                             breaks = c(0, 30))
      ),
      y = list(
        Depth > 4 ~
          scale_y_continuous(limits = c(0, 1.2),
                             breaks = seq(0, 1.2, 0.6),
                             labels = scales::label_number(accuracy = c(1, 0.1, 0.1))),
        Depth > 2 & Depth < 4 ~
          scale_y_continuous(limits = c(0, 3),
                             breaks = seq(0, 3, 1.5),
                             labels = scales::label_number(accuracy = c(1, 0.1, 1))),
        Depth == 1.5 ~
          scale_y_continuous(limits = c(0, 5),
                             breaks = seq(0, 5, 2.5),
                             labels = scales::label_number(accuracy = c(1, 0.1, 1))),
        Depth == 0.5 ~
          scale_y_continuous(limits = c(0, 100),
                             breaks = seq(0, 100, 50))
      )
    ) +
    labs(x = "Time since kelp addition (days)",
         y = "Remaining kelp (%)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 0.1, 0.2, 0.2, unit = "cm"))

Fig_2

Fig_2 %>%
  ggsave(filename = "Fig_2.pdf", path = "Figures",
         device = cairo_pdf, height = 24, width = 24, units = "cm")

# 2.3 Figure 3 ####
# Proportion over time
Fig_3a <- time_prediction %>%
  filter(Species != "Prior") %>%
  ggplot() +
    geom_point(data = sed_sum %>% filter(Day != 0),
               aes(Day, Proportion * 100, colour = Species), 
               shape = 16, size = 2.5) +
    geom_line(aes(Day, Proportion * 100, colour = Species)) +
    geom_ribbon(aes(Day, ymin = Proportion.lower * 100, 
                    ymax = Proportion.upper * 100,
                    alpha = factor(.width),
                    fill = Species)) +
    scale_colour_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                        guide = "none") +
    scale_fill_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                      guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(limits = c(0, 100)) +
    facet_nested(
      ~ Species + Sediment + Load, 
      nest_line = TRUE, solo_line = TRUE,
      scales = "free", space = "free_x",
      strip = strip_nested(
        clip = "off", by_layer_x = TRUE,
        text_x = list(element_text(face = "italic"), NULL, NULL)
      ),
      labeller = labeller(
        Load = c(
          "1" = "1 g",
          "7.5" = "7.5 g",
          "10" = "10 g",
          "15" = "15 g"
        )
      )
    ) +
    facetted_pos_scales(
      x = list(
        Species == "Ecklonia radiata" ~
          scale_x_continuous(limits = c(0, 120),
                             breaks = seq(0, 120, 30)),
        Species %>% str_detect("Laminaria") ~
          scale_x_continuous(limits = c(0, 30),
                             breaks = c(0, 30))
      )
    ) +
    labs(x = "Time since kelp addition (days)",
         y = "Remaining kelp (%)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"))

Fig_3a

# Half-life
Fig_3b <- time_prior_posterior %>%
  filter(Species != "Prior") %>%
  ggplot() +
    geom_density(aes(t0.5, fill = Species),
                 n = 2^10, bw = 6 * 0.02, colour = NA) +
    scale_fill_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                      guide = "none") +
    scale_x_continuous(limits = c(0, 6),
                       breaks = seq(0, 6, 3),
                       oob = scales::oob_keep) +
    facet_nested(
      ~ Species + Sediment + Load, 
      nest_line = TRUE, solo_line = TRUE,
      strip = strip_nested(
        clip = "off", by_layer_x = TRUE,
        text_x = list(element_text(face = "italic"), NULL, NULL)
      ),
      labeller = labeller(
        Load = c(
          "1" = "1 g",
          "7.5" = "7.5 g",
          "10" = "10 g",
          "15" = "15 g"
        ),
        Sediment = c( # Abbreviate sand
          "Coarse sand" = "Coarse",
          "Fine sand" = "Fine"
        )
      )
    ) +
    labs(x = "Half-life of labile kelp (days)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 0.5, 0, 0.2, unit = "cm"),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())

Fig_3b

# Refractory proportion
Fig_3c <- time_prior_posterior %>%
  filter(Species != "Prior") %>%
  ggplot() +
    geom_density(aes(r * 100, fill = Species),
                 n = 2^10, bw = 20 * 0.02, colour = NA) +
    scale_fill_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                      guide = "none") +
    scale_x_continuous(limits = c(0, 20),
                       breaks = seq(0, 20, 10),
                       oob = scales::oob_keep) +
    facet_grid( ~ Species + Sediment + Load ) +
    labs(x = "Putative refractory kelp (%)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_blank())

Fig_3c

require(patchwork)
Fig_3 <- ( Fig_3a / Fig_3b / Fig_3c ) +
          plot_layout(heights = c(1, 0.3, 0.3))
Fig_3

Fig_3 %>%
  ggsave(filename = "Fig_3.pdf", path = "Figures",
         device = cairo_pdf, height = 16, width = 22, units = "cm")

# 2.4 Figure 4 ####
Fig_4 <- final %>%
  filter(Depth < 10) %>%
  ggplot() +
    geom_col(aes(Depth, C_m2, fill = Species)) +
    # geom_point(aes(Depth, C_m2), colour = "black",
    #            shape = 16, size = 2.5) +
    scale_fill_manual(values = c("#c3b300", "#627d0e", "#f1c700"),
                        guide = "none") +
    scale_x_reverse(limits = c(0, 10), breaks = 0:10) +
    facet_nested(
      ~ Species + Sediment + Load,
      nest_line = TRUE, solo_line = TRUE,
      scales = "free",
      strip = strip_nested(
        clip = "off", by_layer_x = TRUE,
        text_x = list(element_text(face = "italic"), NULL, NULL)
      ),
      labeller = labeller(
        Load = c(
          "1" = "1 g",
          "7.5" = "7.5 g",
          "10" = "10 g",
          "15" = "15 g"
        ),
        Sediment = c( # Abbreviate sand
          "Coarse sand" = "Coarse",
          "Fine sand" = "Fine"
        )
      )
    ) +
    facetted_pos_scales(
      y = list(
        Load == 1 ~
          scale_y_continuous(limits = c(0, 3),
                             breaks = seq(0, 3, 1)),
        Load == 7.5 ~ 
          scale_y_continuous(limits = c(0, 20),
                             breaks = seq(0, 20, 10)),
        Load == 10 ~
          scale_y_continuous(limits = c(0, 20),
                             breaks = seq(0, 20, 10)),
        Load == 15 ~
          scale_y_continuous(limits = c(0, 90),
                             breaks = seq(0, 90, 45))
      )
    ) +
    labs(y = expression("Final kelp carbon stock (g C m"^-2*")"),
         x = "Depth (cm)") +
    coord_flip(expand = FALSE, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"))

Fig_4

Fig_4 %>%
  ggsave(filename = "Fig_4.pdf", path = "Figures",
         device = cairo_pdf, height = 12, width = 22, units = "cm")

# Calculate corresponding carbon influx
kelp_carb %>% 
  filter(Date == "26.09.25" & Species != "Saccharina latissima") %>%
  summarise(C = mean(C/100), .by = Species) %>%
  full_join(sed %>% distinct(Species, Sediment, Load)) %>%
  mutate(Influx = Load * C / ( pi * 0.05^2 ))

# Calculate baseline carbon stocks
sed %>%
  filter(Day == 0) %>%
  mutate(C_m2 = Mass_c * Corg / ( pi * 0.05^2 )) %>%
  summarise(C_m2_mean = mean(C_m2),
            C_m2_sd = sd(C_m2),
            n = n(),
            .by = Sediment) %>%
  mutate(
    across(c(ends_with("mean"), ends_with("sd")), ~ signif(.x, 2)),
    C_m2 = glue("{C_m2_mean} ± {C_m2_sd} (n = {n})")
  )

