require(tidyverse)

nafluo1 <- read.csv("~/Desktop/PhD/Papers/Burial/Data/Pulse-chase/Bioirrigation/240317_NaFluo_1.csv", skip = 10) %>%
  mutate(Tank = c(rep(c(1, 2), each = 3 * 15), rep(NA, 6)),
         Core = c(rep(rep(seq(1, 15), each = 3), 2), rep(NA, 6)),
         Concentration = c(rep(NA, 90), c(0, 1, 2, 5, 10, 20))) %>%
  rename(Fluorescence = 3)

nafluo2 <- read.csv("~/Desktop/PhD/Papers/Burial/Data/Pulse-chase/Bioirrigation/240317_NaFluo_2.csv", skip = 10) %>%
  mutate(Tank = c(rep(c(1, 2), each = 3 * 15), rep(NA, 6)),
         Core = c(rep(rep(seq(1, 15), each = 3), 2), rep(NA, 6)),
         Concentration = c(rep(NA, 90), c(0, 1, 2, 5, 10, 20))) %>%
  rename(Fluorescence = 3)

nafluo3 <- read.csv("~/Desktop/PhD/Papers/Burial/Data/Pulse-chase/Bioirrigation/240317_NaFluo_3.csv", skip = 10) %>%
  mutate(Tank = c(rep(c(1, 2), each = 3 * 15), rep(NA, 6)),
         Core = c(rep(rep(seq(1, 15), each = 3), 2), rep(NA, 6)),
         Concentration = c(rep(NA, 90), c(0, 1, 2, 5, 10, 20))) %>%
  rename(Fluorescence = 3)

nafluo4 <- read.csv("~/Desktop/PhD/Papers/Burial/Data/Pulse-chase/Bioirrigation/240317_NaFluo_4.csv", skip = 10) %>%
  mutate(Tank = c(rep(c(1, 2), each = 3 * 15), rep(NA, 6)),
         Core = c(rep(rep(seq(1, 15), each = 3), 2), rep(NA, 6)),
         Concentration = c(rep(NA, 90), c(0, 1, 2, 5, 10, 20))) %>%
  rename(Fluorescence = 3)


standard <- bind_rows(nafluo1 %>% filter(!is.na(Concentration)),
                      nafluo2 %>% filter(!is.na(Concentration)),
                      nafluo3 %>% filter(!is.na(Concentration)),
                      nafluo4 %>% filter(!is.na(Concentration))) %>%
  select(Concentration, Fluorescence) %>%
  mutate(Plate = rep(seq(1, 4), each = 6))

standard %>%
  ggplot(aes(Concentration, Fluorescence)) +
    geom_point() +
    geom_smooth(method = "lm", se = F, colour = "black") +
    facet_grid(~Plate) +
    theme_minimal()


curve <- standard %>%
  nest_by(Plate) %>%
  mutate(mod = list(lm(Fluorescence ~ Concentration, data = data))) %>%
  reframe(broom::tidy(mod)) %>%
  select(Plate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(Intercept = "(Intercept)", Slope = Concentration)


samples <- bind_rows(nafluo1 %>% filter(is.na(Concentration)),
                     nafluo2 %>% filter(is.na(Concentration)),
                     nafluo3 %>% filter(is.na(Concentration)),
                     nafluo4 %>% filter(is.na(Concentration))) %>%
  select(Tank, Core, Fluorescence) %>%
  mutate(Plate = rep(seq(1, 4), each = 90)) %>%
  left_join(curve, by = "Plate") %>%
  mutate(Concentration = (Fluorescence - Intercept) / Slope)

samples %>%
  ggplot(aes(Plate, Concentration)) +
    geom_point() +
    geom_smooth(method = "lm", se = F, colour = "black") +
    facet_grid(~ Tank + Core) +
    theme_minimal()


