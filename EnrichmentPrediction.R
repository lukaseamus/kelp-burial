#### Step 1: Calculate required 13C in seawater to enrich kelp ####

### Baseline d13C for Ecklonia radiata ###
d13C <- c(
        # Data from Stepien (2015) 10.1111/1365-2745.12451 Table S1
        -21.26, -19.26, -21.45, -18.4, -18.7, -16.7, -23.26, -23.29, -18.1,
        -16.06, -15.64, -22.15, -19.82, -18.79, -17.77, -17.69, -15.94,
        # Data from Vanderklift and Bearham (2014) 10.3354/meps10967 Table A1
        -20.7, -20.3, -18.1, -18.0, -22.0, -21.3, -22.5, -21.2, -20.9, -19.0, 
        -22.1, -21.5, -18.2, -16.9, -20.7, -19.1, -20.8, -19.2, -17.5, -17.4, 
        -17.3, -18.2, -21.2, -20.7, -17.6, -17.2, -19.3, -18.2, -19.3, -17.5, 
        -20.1, -19.5, -18.7, -18.3, -19.7, -19.1, -23.0, -20.5, -21.3, -19.7,
        -21.7, -18.7, -23.2, -20.9, -23.0, -19.0, -20.6, -21.4, -21.1, -17.7, 
        -18.5, -17.8, -19.7, -18.9, -19.4, -16.6, -21.0, -20.1, -19.7, -18.3, 
        -20.4, -19.3, -19.5, -18.7, -19.9, -19.4, -17.9, -16.7, -18.8, -18.0, 
        -18.3, -17.1, -20.1, -21.5, -21.2, -22.4, -22.6, -22.2, -23.0, -23.4, 
        -22.9, -23.8, -21.0, -19.3, -19.5, -20.4, -19.8, -18.7, -19.3, -18.1,
        -18.5, -18.4, -22.1, -20.5, -21.3, -22.4, -19.4, -19.6, -20.2, -20.5, 
        -19.0, -18.1, -19.2, -18.7, -19.4, -19.3, -18.8, -18.6, -20.4, -20.6
        )

# convert d13C to 13C atom percent
# equation: 13C at% = (d13C/1000 + 1) * 0.0112372 / 1 + ((d13C/1000 + 1) * 0.0112372)
# where 0.0112372 is the R(13C/12C)VPDB value (Table 2 in Skrzypek & Dunn 2020 10.1002/rcm.8892)
at <- function(d) (((d/1000 + 1) * 0.0112372) / (1 + ((d/1000 + 1) * 0.0112372))) * 100

at13C <- at(d13C)

mean(at13C) # 1.089566 at%
sd(at13C) # 0.001979763 at%

at13Csim <- rnorm(1e5, mean(at13C), sd(at13C))
plot(density(at13Csim))

### C assimilation by Ecklonia radiata ###
# Gross photosynthesis rates given in mg O2 g-1 dry mass h-1 for our study site from
# Staehr & Wernberg (2009) 10.1111/j.1529-8817.2008.00635.x Fig. 3
# Wernberg et al. (2016) 10.1002/lno.10362 Fig. 2c, g
# Wright et al. (2024) 10.1093/aob/mcad167 Fig. 3

# Gross rather than net rates are more applicable because enriched stable isotope
# is not respired over short timescales (Mateo et al. 2001, 10.3354/meps223157)

# Assuming a mean photosynthetic quotient (O2:CO2 molar ratio) 
# of 1.22 (Miller & Dunton 2007, 10.3354/meps329085), carbon assimilation 
# (mg C g-1 dry mass h-1) is derived by converting mg O2 to mmol O2 (dividing
# by 15.9994 * 2 = 31.9988 g mol-1, note that data from Wright et al. 2024 are
# already given in mmol O2), converting O2 to CO2 (dividing by the photosynthetic quotient) 
# and converting from mmol CO2 to mg C (multiplying by 12.011 g mol-1)
P <- c(c(1.30726204, 3.3888029, 4.2097945, 3.9021071, 0.88466563, 1.3691275, 
         2.19846, 2.8872843, 3.3773369, 3.1949424, 3.0655801, 2.9012883) / 31.9988, 
       0.084) / 1.22 * 12.011

mean(P) # 0.8372077 mg C g-1 dry mass h-1
sd(P) # 0.312305 mg C g-1 dry mass h-1

Psim <- rnorm(1e5, mean(P), sd(P))
plot(density(Psim))

### Proportional C content of Ecklonia radiata ###
C <- c(
       # Data from Staehr & Wernberg (2009) 10.1111/j.1529-8817.2008.00635.x Table 1
       0.282, 0.308, 0.303, 0.301, 
       # Data from Singh et al. (2021) 10.3389/fmars.2021.678222 Fig. 2
       0.374, 0.368, 0.352, 0.346,
       # Data from Perkins et al. (2022) 10.1007/s10533-022-00946-4 Text
       0.359,
       # Data from Zhang et al. (2020) 10.1016/j.algal.2020.102092 Table 2
       0.326
       )

mean(C) # 0.3319
sd(C) # 0.0321436

Csim <- rnorm(1e5, mean(C), sd(C))
plot(density(Csim))

### Relative C assimilation by Ecklonia radiata ###
# Since carbon assimilation is given per unit total mass, we need to divide 
# carbon assimilation by carbon content to arrive at the carbon assimilated
# by the carbon fraction of total mass (mg C g-1 C h-1)
# This is then expressed as a proportional rate (g C g-1 C h-1) by dividing by 1000
rPsim <- Psim / Csim / 1e3

mean(rPsim) # 0.002548876 proportion h-1
sd(rPsim) # 0.0009888362 proportion h-1
# these values may deviate somewhat since this is a probabilistic calculation
plot(density(rPsim))

### Calculations ###
# e = 13C enriched kelp (at%)
# b = 13C baseline kelp (at%)
# w = 13C enriched water (at%)
# p = relative C assimilation (proportion h-1)
# t = incubation time under irradiance (h)

# estimate water enrichment required to reach desired at% after t h incubation
# importantly, the added C proportional to the final carbon mass counts, which
# is calculated as p * t / (1 + p * t)

w <- function(p, t, b, e) (e - (1 - p * t / (1 + p * t)) * b) / (p * t / (1 + p * t))
# calculate for 16:8-h light:dark cycle for one week (7 d)
wsim <- w(p = rPsim, t = 7 * 16, b = at13Csim, e = 10)

median(wsim) # 41.41471 at%

plot(density(wsim))
# plot is not satisfactory because this is a quotient distribution
require(tidyverse) # use ggplot2 instead
ggplot() +
  geom_vline(xintercept = median(wsim)) +
  geom_density(aes(wsim)) +
  scale_x_continuous(limits = c(0, 100),
                     oob = scales::oob_keep) +
  xlab(expression("Required seawater "^13*"C (at%)")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

# this is the inverse function (simple mixing model)
e <- function(p, t, b, w) (1 - p * t / (1 + p * t)) * b + (p * t / (1 + p * t)) * w

# test function reversibility using the result from above
e(p = rPsim, t = 7 * 16, b = at13Csim, w = wsim) # all values converge to 10 at%
plot(density(e(p = rPsim, t = 7 * 16, b = at13Csim, w = wsim) - 10)) # visualise difference = 0
# the estimated and expected final enrichment are identical, so the function works

# calculate an example case where water has exactly 41 at%
esim <- e(p = rPsim, t = 7 * 16, b = at13Csim, w = 41)
plot(density(esim)) # most probability density lies around 10 at%

# seawater dDI13C (‰)
dDI13Cdeepsim <- rnorm(1e5, 0.88, 0.05) # Table 3 in Cheng et al. 2019 10.1002/lom3.10300
plot(density(dDI13Cdeepsim))

# however, this is for deep waters and shallow waters seem to have higher d13C
dDI13Cshallowsim <- rnorm(1e5, 1.5, 0.2) # Figure 4b for -32° in Quay et al. 2003 10.1029/2001GB001817
plot(density(dDI13Cshallowsim))

atDI13Csim <- at(dDI13Cshallowsim)
mean(atDI13Csim) # 1.112881 at%
sd(atDI13Csim) # 0.0002192243 at%

#### Step 2: Calculate required 13C sodium bicarbonate to enrich seawater ####

# calculate proportion of 13C coming from label using a further derivation of the simple 
# mixing models that are w() and e()
# e = 13C enriched seawater (at%)
# b = 13C baseline seawater (at%)
# l = 13C NaH13CO3 (at%)
# p = label proportion

p <- function(e, b, l) (e - b) / (l - b)

psim <- p(e = wsim, b = atDI13Csim, l = 98)

median(psim) # 0.4159676

ggplot() +
  geom_vline(xintercept = mean(psim), colour= "red") +
  geom_vline(xintercept = median(psim), colour= "blue") +
  geom_density(aes(psim)) +
  scale_x_continuous(limits = c(0, 1),
                     oob = scales::oob_keep) +
  xlab(expression("Required proportion of seawater at%")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())



# calculate average seawater DIC concentration (µmol l-1)
require(seacarb)
DIC <- carb(flag = 21, var1 = x2pCO2(xCO2 = 400), var2 = 8.1)$DIC * 1e6 * 1.025
DIC # 2298.611 µmol l-1

DICsim <- rnorm(1e5, DIC, DIC*0.01) # arbitrary 1% s.d. to add uncertainty
plot(density(DICsim))

# calculate required label addition
# p = proportion of carbon coming from label
# c = concentration of DIC, which is mostly HCO3- (µmol l-1)
# m = molar mass of NaH13CO3 (g mol-1)
# mgl = required label (mg l-1)

mgl <- function(p, c, m) p / (1 - p) * c * m * 1e-3

mglsim <- mgl(p = psim, c = DICsim, m = 85)

mean(mglsim) # 109.7542 mg l-1
median(mglsim) # 132.2828 mg l-1

mglp <- ggplot() +
  geom_vline(xintercept = 115) +
  geom_vline(xintercept = mean(mglsim), colour = "red") +
  geom_vline(xintercept = median(mglsim), colour = "blue") +
  geom_density(aes(mglsim)) +
  scale_x_continuous(limits = c(0, 500),
                     oob = scales::oob_keep) +
  xlab(expression("Required NaH"^13*"CO"[3]*" concentration (mg l"^-1*")")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())
mglp
# total required label for 90 l
mean(mglsim * 90 * 1e-3) # 9.877881 g

gp <- ggplot() +
  geom_vline(xintercept = mean(mglsim * 90 * 1e-3), colour = "red") +
  geom_vline(xintercept = median(mglsim * 90 * 1e-3), colour = "blue") +
  geom_vline(xintercept = 0.115 * 90) +
  geom_density(aes(mglsim * 90 * 1e-3)) +
  scale_x_continuous(limits = c(0, 20),
                     oob = scales::oob_keep) +
  xlab(expression("Required NaH"^13*"CO"[3]*" (g) for 90 l")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())
gp

require(patchwork)
mglp | gp
