require(tidyverse)


mug <- read.csv("~/Desktop//240321_4MUG_Corewater_Kinetics.csv", skip = 11) %>%
  rename(Well = X, Sample = Time) %>%
  pivot_longer(cols = -c(Well, Sample), names_to = "Time") %>%
  mutate(min = as.numeric(str_extract(Time, pattern = "\\d+")),
         s = as.numeric(str_extract(Time, pattern = "\\b\\d{2}\\b")),
         s = ifelse(is.na(s), 0, s),
         Time = min + s / 60,
         variable = rep(rep(c("Fluorescence", "Temperature"), each = 450), 96)) %>%
  select(-c(min, s)) %>%
  pivot_wider(names_from = variable, values_from = value)
  

slopes <- mug %>% 
  slice(1:37800) %>%
  nest_by(Well) %>%
  mutate(slope = list(coef(lm(Fluorescence ~ Time, data = data)))) %>%
  unnest_wider(slope) %>%
  select(-data) %>%
  rename(x0 = `(Intercept)`, y = Time) %>%
  mutate(Concentration = rep(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500), 7),
         Sample = rep(c("1.1", "1.2", "1.3", "2.1", "2.2", "2.3", "blank"), each =  12))


ggplot(data = slopes, aes(Concentration, y, colour = Sample)) +
  geom_point() +
  geom_path()



lamc <- read.csv("~/Desktop//240321_LAMC_Corewater_Kinetics.csv", skip = 11) %>%
  rename(Well = X, Sample = Time) %>%
  pivot_longer(cols = -c(Well, Sample), names_to = "Time") %>%
  mutate(min = as.numeric(str_extract(Time, pattern = "\\d+")),
         s = as.numeric(str_extract(Time, pattern = "\\b\\d{2}\\b")),
         s = ifelse(is.na(s), 0, s),
         Time = min + s / 60,
         variable = rep(rep(c("Fluorescence", "Temperature"), each = 450), 96)) %>%
  select(-c(min, s)) %>%
  pivot_wider(names_from = variable, values_from = value)


slopes2 <- lamc %>% 
  slice(1:37800) %>%
  nest_by(Well) %>%
  mutate(slope = list(coef(lm(Fluorescence ~ Time, data = data)))) %>%
  unnest_wider(slope) %>%
  select(-data) %>%
  rename(x0 = `(Intercept)`, y = Time) %>%
  mutate(Concentration = rep(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500), 7),
         Sample = rep(c("1.1", "1.2", "1.3", "2.1", "2.2", "2.3", "blank"), each =  12))


ggplot(data = slopes2, aes(Concentration, y, colour = Sample)) +
  geom_point() +
  geom_path()







mug <- read.csv("~/Desktop//240321_4MUG_Corewater_Kinetics.csv", header = FALSE) %>%
  slice(-1:-10) %>% select(1:452) %>% `colnames<-`(.[1, ]) %>% .[-1, ] %>% `rownames<-`(NULL)
  rename(Row = 1, Column = 2, Read = 3)
  pivot_longer(cols = 4:452, values_to = "Fluorescence", names_to = "Time")
  left_join(read.csv("~/Desktop/240321_4MUG_Corewater_Kinetics.csv", header = FALSE) %>%
              slice(-1:-10) %>% select(1:3, 454:903) %>% `colnames<-`(.[1, ]) %>% .[-1, ] %>% `rownames<-`(NULL) %>%
              rename(Row = 1, Column = 2, Read = 3) %>%
              pivot_longer(cols = 4:453, values_to = "Temperature", names_to = "Time"),
            by = c("Row", "Column", "Read", "Time")) %>%
  mutate(Read = as.numeric(str_extract(Read, "(\\d)+")))

lamc <- read.csv("~/Desktop/PhD/Papers/BvF/Data/230824_LAMC.csv", header = FALSE) %>%
  slice(-1:-10) %>% select(1:453) %>% `colnames<-`(.[1, ]) %>% .[-1, ] %>% `rownames<-`(NULL) %>%
  rename(Row = 1, Column = 2, Read = 3) %>%
  pivot_longer(cols = 4:453, values_to = "Fluorescence", names_to = "Time") %>%
  left_join(read.csv("~/Desktop/PhD/Papers/BvF/Data/230824_LAMC.csv", header = FALSE) %>%
              slice(-1:-10) %>% select(1:3, 454:903) %>% `colnames<-`(.[1, ]) %>% .[-1, ] %>% `rownames<-`(NULL) %>%
              rename(Row = 1, Column = 2, Read = 3) %>%
              pivot_longer(cols = 4:453, values_to = "Temperature", names_to = "Time"),
            by = c("Row", "Column", "Read", "Time")) %>%
  mutate(Read = as.numeric(str_extract(Read, "(\\d)+")))

comb <- data.frame(rbind(mug, lamc), Date = rep("24.08.23", 86400), Substrate = rep(c("4-MUG", "LAMC"), each = 43200)) %>%
  left_join(read.csv("~/Desktop/PhD/Papers/BvF/Data/Saturation.csv"), by = c("Date", "Substrate", "Read")) %>%
  mutate(Column = as.numeric(Column),
         Time = as.numeric(Time),
         Fluorescence = as.numeric(Fluorescence),
         Temperature = as.numeric(Temperature))

standard <- comb %>%
  filter(Fluorophore == "Standard") %>%
  group_by(Substrate, Concentration, Quenching) %>%
  summarise(Fluorescence = mean(Fluorescence))

# Baly equation with intercept term (mathematically identical to Michaelis-Menten)
ggplot(data = standard, aes(Concentration, Fluorescence, colour = Quenching)) +
  geom_point() + 
  geom_smooth(data = standard %>% filter(Substrate == "LAMC" | Quenching == "No"),
              method = "nls", 
              method.args = list(formula = y ~ Fmax * b * x / (Fmax + b * x), # a < 0 here, so removed
                                 start = list(Fmax = 230000, b = 2000)),
              se = FALSE) +
  geom_smooth(data = standard %>% filter(Substrate == "4-MUG" & Quenching == "Yes"),
              method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 5000)) + # disable to see full fit
  facet_grid(~ Substrate)


# solved for x, the Baly equation without an interceot reads 
# x = Fmax * y / (b * (Fmax - y)), which is how the concentration 
# of cleaved substrate can be easily calculated from sample fluorescence

nls(Fluorescence ~ Fmax * b * Concentration / (Fmax + b * Concentration),
    start = list(Fmax = 230000, b = 2000),
    data = standard %>% filter(Substrate == "LAMC" & Quenching == "Yes"))
# y = 225276 * 1070 * x / (225276 + 1070 * x)
# x = 225276 * y / (1070 * (225276 - y))

nls(Fluorescence ~ Fmax * b * Concentration / (Fmax + b * Concentration),
    start = list(Fmax = 230000, b = 2000),
    data = standard %>% filter(Substrate == "LAMC" & Quenching == "No"))
# y = 194275 * 1668 * x / (194275 + 1668 * x)
# x = 194275 * y / (1668 * (194275 - y))

nls(Fluorescence ~ Fmax * b * Concentration / (Fmax + b * Concentration),
    start = list(Fmax = 230000, b = 2000),
    data = standard %>% filter(Substrate == "4-MUG" & Quenching == "No"))
# y = 171350 * 2155 * x / (171350 + 2155 * x)
# x = 171350 * y / (2155 * (171350 - y))

# solving the polynomial is more complicated because a polynomial technically
# cannot be uniquely solved; instead, we find the roots of a modified polynomial,
# the intercept of which is shifted down by the amount of fluorescence we
# are trying to find the equivalent concentration for

poly1 <- lm(Fluorescence ~ poly(Concentration, 3),
            data = standard %>% filter(Substrate == "4-MUG" & Quenching == "Yes"))

new <- data.frame(fit = predict(poly1, newdata = data.frame(Concentration = seq(0, 150, 2))),
                  Concentration = seq(0, 150, 2))

ggplot(data = standard %>% filter(Substrate == "4-MUG" & Quenching == "Yes"),
       aes(Concentration, Fluorescence)) +
  geom_point() +
  geom_line(data = new, aes(Concentration, fit))


polyxsolve <- function(model, order, y, length, output){ 
  
  n <- 1:length
  
  if(order == 2){
    
    rf2 <- Vectorize(function(x){
      polyroot(c(coef(model)[1] - y[x], 
                 coef(model)[2],
                 coef(model)[3]))[1]
    })
    
    root <- as.numeric(rf2(n))
    
  } else if(order == 3){
    
    rf3 <- Vectorize(function(x){
      polyroot(c(coef(model)[1] - y[x], 
                 coef(model)[2],
                 coef(model)[3],
                 coef(model)[4]))[1]
    })
    
    root <- as.numeric(rf3(n))
    
  } else {
    "Not a valid input"
  }
  
  prediction <- as.numeric(predict(model, newdata = data.frame(Concentration = root)))
  
  if(output == "x"){
    root
  } else if(output == "y"){
    prediction
  } else {
    "Not a valid output."
  }
}

polyxsolve(poly1, 3, c(2000, 3000, 4000), 3, "prediction")



nls(Fluorescence ~ Fmax * b * Concentration / (Fmax + b * Concentration) + a,
    start = list(Fmax = 230000, b = 2000, a = 100),
    data = standard %>% filter(Substrate == "4-MUG" & Quenching == "Yes"))
# y = 267100.9 * 770.4 * x / (267100.9 + 770.4 * x) + 148.0
# x = 267100.9 * (y - 148.0) / (770.4 * (148.0 + 267100.9 - y))


substrate <- comb %>% 
  filter(Fluorophore == "Substrate")

ggplot(data = substrate, aes(Time, Fluorescence, colour = Substrate)) +
  geom_point() +
  facet_grid(~ Concentration)




slopes <- comb %>% 
  nest_by(Substrate, ID, Concentration) %>%
  mutate(slope = list(coef(lm(Fluorescence ~ Time, data = data)))) %>%
  unnest_wider(slope) %>%
  select(-data) %>%
  rename(x0 = `(Intercept)`, y = Time)
  
ggplot(data = slopes %>% filter(!ID %in% c("s1", "s1+")),
       aes(Concentration, y, colour = Substrate)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_grid(~ ID)


ggplot(data = slopes %>% filter(ID %in% c("s1", "s1+")),
       aes(Concentration, x0, colour = Substrate)) +
  geom_point() +
  geom_smooth(method = "loess") +
  facet_grid(~ ID)


