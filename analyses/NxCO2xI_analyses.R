## NxCO2xI analysis script. Paths assume the root directory of this script.

##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)
library(multcomp)
library(multcompView)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

# Read in compiled data file
df <- read.csv("../data/NxCO2xI_data.csv") %>%
  mutate(n.trt = as.numeric(n.trt),
         inoc = factor(inoc, levels = c("no.inoc", "inoc")),
         co2 = factor(co2, levels = c("amb", "elv")),
         co2.inoc = str_c(co2, "_", inoc)) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nod.root.ratio < 0.05))
  ## filter all uninoculated pots that have nod biomass > 0.05 g;
  ## hard code inoc/co2 to make coefficients easier to understand

df.removed <- read.csv("../data/NxCO2xI_data.csv") %>%
  dplyr::filter(inoc == "no.inoc" & nod.root.ratio > 0.05)

##########################################################################
## Area-based leaf nitrogen content (Narea)
##########################################################################
narea <- lmer(narea ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(narea)
qqnorm(residuals(narea))
qqline(residuals(narea))
densityPlot(residuals(narea))
shapiro.test(residuals(narea))
outlierTest(narea)

# Model results
summary(narea)
Anova(narea)
r.squaredGLMM(narea)

# Post-hoc tests
test(emtrends(narea, pairwise~co2, "n.trt")) 
test(emtrends(narea, pairwise~inoc, "n.trt"))
emmeans(narea, pairwise~co2*inoc)

# Individual effects
emmeans(narea, pairwise~co2)
emmeans(narea, pairwise~inoc)
test(emtrends(narea, ~1, "n.trt"))

##########################################################################
## Mass-based leaf nitrogen content (Nmass)
##########################################################################
nmass <- lmer(nmass.focal ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
qqnorm(residuals(nmass))
qqline(residuals(nmass))
densityPlot(residuals(nmass))
shapiro.test(residuals(nmass))
outlierTest(nmass)

# Model results
summary(nmass)
Anova(nmass)
r.squaredGLMM(nmass)

# Post-hoc tests
test(emtrends(nmass, pairwise~co2, "n.trt")) 
test(emtrends(nmass, pairwise~inoc, "n.trt"))

# Individual effect of co2
emmeans(nmass, pairwise~co2)
emmeans(nmass, pairwise~inoc)
test(emtrends(nmass, ~1, "n.trt"))

##########################################################################
## Leaf mass per unit leaf area (Marea)
##########################################################################
marea <- lmer(log(marea) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
qqnorm(residuals(marea))
qqline(residuals(marea))
densityPlot(residuals(marea))
shapiro.test(residuals(marea))
outlierTest(marea)

# Model results
summary(marea)
Anova(marea)
r.squaredGLMM(marea)

# Post-hoc tests
test(emtrends(marea, pairwise~co2, "n.trt")) 
test(emtrends(marea, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(marea, pairwise~co2, type = "response")
emmeans(marea, pairwise~inoc)
test(emtrends(marea, ~1, "n.trt"))

##########################################################################
## Chlorophyll content (Chlarea)
##########################################################################
chlarea <- lmer(chl.mmolm2 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
qqnorm(residuals(chlarea))
qqline(residuals(chlarea))
densityPlot(residuals(chlarea))
shapiro.test(residuals(chlarea))
outlierTest(chlarea)

# Model results
summary(chlarea)
Anova(chlarea)
r.squaredGLMM(chlarea)

# Post-hoc tests
test(emtrends(chlarea, pairwise~co2, "n.trt")) 
test(emtrends(chlarea, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(chlarea, pairwise~co2, type = "response")
emmeans(chlarea, pairwise~inoc)
test(emtrends(chlarea, ~1, "n.trt"))

##########################################################################
## Net photosynthesis at 420ppm CO2 (Anet,420)
##########################################################################
anet <- lmer(anet ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(anet)
qqnorm(residuals(anet))
qqline(residuals(anet))
densityPlot(residuals(anet))
shapiro.test(residuals(anet))
outlierTest(anet)

# Model results
summary(anet)
Anova(anet)
r.squaredGLMM(anet)

# Pairwise comparisons
test(emtrends(anet, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(anet, pairwise~co2)
cld(emmeans(anet, pairwise~inoc))
test(emtrends(anet, ~1, "n.trt"))

##########################################################################
## Net photosynthesis at growth CO2 concentration (Anet,growth)
##########################################################################
anet.growth <- lmer(anet.growth ~ co2 * inoc * n.trt + (1|rack:co2),
                    data = df)

# Check model fit
plot(anet.growth)
qqnorm(residuals(anet.growth))
qqline(residuals(anet.growth))
densityPlot(residuals(anet.growth))
shapiro.test(residuals(anet.growth))
outlierTest(anet.growth)

# Model results
summary(anet.growth)
Anova(anet.growth)
r.squaredGLMM(anet.growth)

# Pairwise comparisons
emmeans(anet.growth, pairwise~co2*inoc)
test(emtrends(anet.growth, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(anet.growth, pairwise~co2)
cld(emmeans(anet.growth, pairwise~inoc))
test(emtrends(anet.growth, ~1, "n.trt"))

##########################################################################
## Maximum Rubisco carboxylation rate (Vcmax25)
##########################################################################
vcmax <- lmer(vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(vcmax)
qqnorm(residuals(vcmax))
qqline(residuals(vcmax))
densityPlot(residuals(vcmax))
shapiro.test(residuals(vcmax))
outlierTest(vcmax)

# Model results
summary(vcmax)
Anova(vcmax)
r.squaredGLMM(vcmax)

# Pairwise comparisons
test(emtrends(vcmax, pairwise~inoc, "n.trt"))

# Individual effect of CO2
emmeans(vcmax, pairwise~co2)
emmeans(vcmax, pairwise~inoc)
test(emtrends(vcmax, ~1, "n.trt"))

##########################################################################
## Maximum electron transport for RuBP regeneration rate (Jmax25)
##########################################################################
jmax <- lmer(jmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(jmax)
qqnorm(residuals(jmax))
qqline(residuals(jmax))
densityPlot(residuals(jmax))
shapiro.test(residuals(jmax))
outlierTest(jmax)

# Model results
summary(jmax)
Anova(jmax)
r.squaredGLMM(jmax)

# Pairwise comparisons
test(emtrends(jmax, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(jmax, pairwise~co2)
emmeans(jmax, pairwise~inoc)
test(emtrends(jmax, ~1, "n.trt"))

# Percent change for study limitation section
emmeans(jmax, ~1, "n.trt", at = list(n.trt = c(0, 630)))

##########################################################################
## Ratio of Jmax25 to Vcmax25
##########################################################################
df$jmax25.vcmax25[100] <- NA
jvmax <- lmer(jmax25.vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(jvmax)
qqnorm(residuals(jvmax))
qqline(residuals(jvmax))
densityPlot(residuals(jvmax))
shapiro.test(residuals(jvmax))
outlierTest(jvmax)

# Model results
summary(jvmax)
Anova(jvmax)
r.squaredGLMM(jvmax)

# Pairwise comparisons
test(emtrends(jvmax, pairwise~inoc, "n.trt"))
test(emtrends(jvmax, pairwise~co2, "n.trt"))

# Individual effects
emmeans(jvmax, pairwise~co2)
emmeans(jvmax, pairwise~inoc)
test(emtrends(jvmax, ~1, "n.trt"))

##########################################################################
## Dark respiration (Rd25)
##########################################################################
df$rd25[df$rd25 < 0] <- NA
df$rd25[c(29, 34, 56)] <- NA

rd25 <- lmer(rd25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(rd25)
qqnorm(residuals(rd25))
qqline(residuals(rd25))
densityPlot(residuals(rd25))
shapiro.test(residuals(rd25))
outlierTest(rd25)

# Model results
summary(rd25)
Anova(rd25)
r.squaredGLMM(rd25)

# Pairwise comparisons
test(emtrends(rd25, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(rd25, pairwise~co2)
emmeans(rd25, pairwise~inoc)
test(emtrends(rd25, ~1, "n.trt"))

##########################################################################
## Photosynthetic nitrogen-use efficiency (PNUE)
##########################################################################
df$pnue.growth[41] <- NA

pnue <- lmer(pnue.growth ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(pnue)
qqnorm(residuals(pnue))
qqline(residuals(pnue))
densityPlot(residuals(pnue))
shapiro.test(residuals(pnue))
outlierTest(pnue)

# Model results
summary(pnue)
Anova(pnue)
r.squaredGLMM(pnue)

# Pairwise comparisons
test(emtrends(pnue, pairwise~co2, "n.trt"))
test(emtrends(pnue, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(pnue, pairwise~co2)
emmeans(pnue, pairwise~inoc)

##########################################################################
## Leaf Ci:Ca (chi)
##########################################################################
chi <- lmer(chi ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(chi)
qqnorm(residuals(chi))
qqline(residuals(chi))
densityPlot(residuals(chi))
shapiro.test(residuals(chi))
outlierTest(chi)

# Model results
summary(chi)
Anova(chi)
r.squaredGLMM(chi)

# Pairwise comparisons
test(emtrends(chi, pairwise~inoc*co2, "n.trt"))
emmeans(chi, pairwise~inoc*co2)
test(emtrends(chi, pairwise~inoc, "n.trt"))
test(emtrends(chi, pairwise~co2, "n.trt"))

## Individual effect of n.trt on chi
test(emtrends(chi, ~1, "n.trt"))
emmeans(chi, pairwise~inoc)
emmeans(chi, pairwise~co2)

##########################################################################
## Total leaf area
##########################################################################
tla <- lmer(tla ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(tla)
qqnorm(residuals(tla))
qqline(residuals(tla))
densityPlot(residuals(tla))
shapiro.test(residuals(tla))
outlierTest(tla)

# Model results
summary(tla)
Anova(tla)
r.squaredGLMM(tla)

# Pairwise comparisons
test(emtrends(tla, pairwise~co2, "n.trt"))
test(emtrends(tla, pairwise~inoc, "n.trt"))
emmeans(tla, pairwise~co2*inoc)

## Individual effects
emmeans(tla, pairwise~co2)
emmeans(tla, pairwise~inoc)

## Does inoculation stimulate TLA under low soil N?
emmeans(tla, pairwise~co2*inoc, "n.trt", type = "response",
        at = list(n.trt = c(0,630)))

##########################################################################
## Total biomass
##########################################################################
tbio <- lmer(sqrt(total.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(tbio)
qqnorm(residuals(tbio))
qqline(residuals(tbio))
densityPlot(residuals(tbio))
shapiro.test(residuals(tbio))
outlierTest(tbio)

# Model results
summary(tbio)
Anova(tbio)
r.squaredGLMM(tbio)

# Pairwise comparisons
test(emtrends(tbio, pairwise~inoc, "n.trt"))
test(emtrends(tbio, pairwise~co2, "n.trt"))
emmeans(tbio, pairwise~co2*inoc)

## Individual effects
emmeans(tbio, pairwise~co2, type = "response")
emmeans(tbio, pairwise~inoc)
cld(emmeans(tbio, pairwise~co2*inoc))

##########################################################################
## Structural carbon cost to acquire nitrogen (Ncost)
##########################################################################
df$ncost[c(100, 101)] <- NA
df$ncost[c(38, 103)] <- NA
df$ncost[32] <- NA

ncost <- lmer(ncost ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(ncost)
qqnorm(residuals(ncost))
qqline(residuals(ncost))
densityPlot(residuals(ncost))
shapiro.test(residuals(ncost))
outlierTest(ncost)

# Model results
summary(ncost)
Anova(ncost)
r.squaredGLMM(ncost)

# Pairwise comparisons
test(emtrends(ncost, pairwise~co2*inoc, "n.trt"))
test(emtrends(ncost, pairwise~co2, "n.trt"))
test(emtrends(ncost, pairwise~inoc, "n.trt"))
emmeans(ncost, pairwise~co2*inoc)

## Individual effects
emmeans(ncost, pairwise~co2)
emmeans(ncost, pairwise~inoc)
test(emtrends(ncost, ~1, "n.trt"))

##########################################################################
## Belowground carbon biomass
##########################################################################
cbg <- lmer(log(cbg) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(cbg)
qqnorm(residuals(cbg))
qqline(residuals(cbg))
densityPlot(residuals(cbg))
shapiro.test(residuals(cbg))
outlierTest(cbg)

# Model results
summary(cbg)
Anova(cbg)
r.squaredGLMM(cbg)

# Pairwise comparisons
test(emtrends(cbg, pairwise~inoc, "n.trt"))
cld(emmeans(cbg, pairwise~co2*inoc, type = "response"))

## Individual effects
emmeans(cbg, pairwise~co2, type = "response")
emmeans(cbg, pairwise~inoc, type = "response")
test(emtrends(cbg, ~1, "n.trt"))

##########################################################################
## Root biomass
##########################################################################
root.bio <- lmer(root.biomass ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(root.bio)
qqnorm(residuals(root.bio))
qqline(residuals(root.bio))
densityPlot(residuals(root.bio))
shapiro.test(residuals(root.bio))
outlierTest(root.bio)

# Model results
summary(root.bio)
Anova(root.bio)
r.squaredGLMM(root.bio)

# Pairwise comparisons
emmeans(root.bio, pairwise~co2*inoc)
test(emtrends(root.bio, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(root.bio, pairwise~co2)
emmeans(root.bio, pairwise~inoc)
test(emtrends(root.bio, ~1, "n.trt"))

##########################################################################
## Whole plant nitrogen
##########################################################################
df$wpn[c(92)] <- NA
wpn <- lmer(sqrt(wpn) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(wpn)
qqnorm(residuals(wpn))
qqline(residuals(wpn))
densityPlot(residuals(wpn))
shapiro.test(residuals(wpn))
outlierTest(wpn)

# Model results
summary(wpn)
Anova(wpn)
r.squaredGLMM(wpn)

# Pairwise comparisons
test(emtrends(wpn, pairwise~inoc, "n.trt"))
test(emtrends(wpn, pairwise~co2, "n.trt"))

## Individual effects
emmeans(wpn, pairwise~co2, type = "response")
emmeans(wpn, pairwise~inoc, type = "response")
test(emtrends(wpn, ~1, "n.trt", regrid = "response"))

##########################################################################
## Root nodule biomass: root biomass
##########################################################################
df$nodule.biomass[80] <- NA

nod.bio <- lmer(sqrt(nodule.biomass) ~ co2 * inoc * n.trt + 
                  (1|rack:co2), data = df)

# Check model fit
plot(nod.bio)
qqnorm(residuals(nod.bio))
qqline(residuals(nod.bio))
densityPlot(residuals(nod.bio))
shapiro.test(residuals(nod.bio))
outlierTest(nod.bio)

# Model results
summary(nod.bio)
Anova(nod.bio)
r.squaredGLMM(nod.bio)

# Pairwise comparisons
test(emtrends(nod.bio, pairwise~inoc, "n.trt"))
test(emtrends(nod.bio, pairwise~co2, "n.trt"))

## Individual effects
emmeans(nod.bio, pairwise~co2)
emmeans(nod.bio, pairwise~inoc)
test(emtrends(nod.bio, ~1, "n.trt"))

##########################################################################
## Root nodule biomass: root biomass
##########################################################################
nod.root.ratio <- lmer(sqrt(nod.root.ratio) ~ co2 * inoc * n.trt + 
                         (1|rack:co2), data = df)

# Check model fit
plot(nod.root.ratio)
qqnorm(residuals(nod.root.ratio))
qqline(residuals(nod.root.ratio))
densityPlot(residuals(nod.root.ratio))
shapiro.test(residuals(nod.root.ratio))
outlierTest(nod.root.ratio)

# Model results
summary(nod.root.ratio)
Anova(nod.root.ratio)
r.squaredGLMM(nod.root.ratio)

# Pairwise comparisons
test(emtrends(nod.root.ratio, pairwise~inoc, "n.trt"))
test(emtrends(nod.root.ratio, pairwise~co2, "n.trt"))
emmeans(nod.root.ratio, pairwise~co2*inoc, type = "response")

## Individual effects
emmeans(nod.root.ratio, pairwise~co2)
emmeans(nod.root.ratio, pairwise~inoc)
test(emtrends(nod.root.ratio, ~1, "n.trt"))

# Percent change in inoculated pots
emmeans(nod.root.ratio, ~inoc, "n.trt", at = list(n.trt = c(0, 630)))

##########################################################################
## Percent leaf nitrogen fixed from the atmosphere (%Ndfa)
##########################################################################
df$ndfa[c(38, 85, 101, 103)] <- NA

ndfa <- lmer(sqrt(ndfa) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
 
# Check model fit
plot(ndfa)
qqnorm(residuals(ndfa))
qqline(residuals(ndfa))
densityPlot(residuals(ndfa))
shapiro.test(residuals(ndfa))
outlierTest(ndfa)

# Model results
summary(ndfa)
Anova(ndfa)
r.squaredGLMM(ndfa)

# Pairwise comparisons
test(emtrends(ndfa, ~inoc, "n.trt"))
test(emtrends(ndfa, ~co2, "n.trt"))
test(emtrends(ndfa, ~1, "n.trt"))
emmeans(ndfa, pairwise~inoc)

emmeans(ndfa, pairwise~co2, at = list(n.ppm = 0))

##########################################################################
## Table 1: Leaf N content
##########################################################################
narea.coefs <- data.frame(summary(narea)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.narea = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.narea) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.narea) %>%
  print(., row.names = FALSE)

narea.table <- data.frame(Anova(narea)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(narea.coefs) %>%
  dplyr::select(treatment, df = Df, coef.narea, 
                chisq.narea = Chisq, pval.narea = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.narea:pval.narea, round, digits = 3),
         chisq.narea = ifelse(chisq.narea <0.001 & chisq.narea >= 0, 
                              "<0.001", chisq.narea),
         pval.narea = ifelse(pval.narea <0.001 & pval.narea >= 0, 
                             "<0.001", pval.narea)) %>%
  arrange(treatment)

nmass.coefs <- data.frame(summary(nmass)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.nmass = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.nmass) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.nmass) %>%
  print(., row.names = FALSE)

nmass.table <- data.frame(Anova(nmass)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(nmass.coefs) %>%
  dplyr::select(treatment, df = Df, coef.nmass, 
                chisq.nmass = Chisq, pval.nmass = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.nmass:pval.nmass, round, digits = 3),
         chisq.nmass = ifelse(chisq.nmass <0.001 & chisq.nmass >= 0, 
                              "<0.001", chisq.nmass),
         pval.nmass = ifelse(pval.nmass <0.001 & pval.nmass >= 0, 
                             "<0.001", pval.nmass)) %>%
  arrange(treatment)

marea.coefs <- data.frame(summary(marea)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.marea = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.marea) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.marea) %>%
  print(., row.names = FALSE)

marea.table <- data.frame(Anova(marea)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(marea.coefs) %>%
  dplyr::select(treatment, df = Df, coef.marea, 
                chisq.marea = Chisq, pval.marea = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.marea:pval.marea, round, digits = 3),
         chisq.marea = ifelse(chisq.marea <0.001 & chisq.marea >= 0, 
                              "<0.001", chisq.marea),
         pval.marea = ifelse(pval.marea <0.001 & pval.marea >= 0, 
                             "<0.001", pval.marea)) %>%
  arrange(treatment)

chl.area.coefs <- data.frame(summary(chlarea)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.chl.area = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.chl.area) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.chl.area) %>%
  print(., row.names = FALSE)

chl.area.table <- data.frame(Anova(chlarea)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(chl.area.coefs) %>%
  dplyr::select(treatment, df = Df, coef.chl.area, 
                chisq.chl.area = Chisq, pval.chl.area = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.chl.area:pval.chl.area, round, digits = 3),
         chisq.chl.area = ifelse(chisq.chl.area <0.001 & chisq.chl.area >= 0, 
                              "<0.001", chisq.chl.area),
         pval.chl.area = ifelse(pval.chl.area <0.001 & pval.chl.area >= 0, 
                             "<0.001", pval.chl.area)) %>%
  arrange(treatment)

table1 <- narea.table %>% full_join(nmass.table) %>% 
  full_join(marea.table) %>% full_join(chl.area.table) %>%
  replace(is.na(.), "-")

##########################################################################
## Table 2: Gas exchange
##########################################################################
anet420.coefs <- data.frame(summary(anet)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.anet420 = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.anet420) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.anet420) %>%
  print(., row.names = FALSE)

anet420.table <- data.frame(Anova(anet)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(anet420.coefs) %>%
  dplyr::select(treatment, df = Df, coef.anet420, 
                chisq.anet420 = Chisq, pval.anet420 = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.anet420:pval.anet420, round, digits = 3),
         chisq.anet420 = ifelse(chisq.anet420 < 0.001 & chisq.anet420 >= 0, 
                              "<0.001", chisq.anet420),
         pval.anet420 = ifelse(pval.anet420 < 0.001 & pval.anet420 >= 0, 
                             "<0.001", pval.anet420)) %>%
  arrange(treatment)

anetgrowth.coefs <- data.frame(summary(anet.growth)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.anetgrowth = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.anetgrowth) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.anetgrowth) %>%
  print(., row.names = FALSE)

anetgrowth.table <- data.frame(Anova(anet.growth)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(anetgrowth.coefs) %>%
  dplyr::select(treatment, df = Df, coef.anetgrowth, 
                chisq.anetgrowth = Chisq, pval.anetgrowth = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.anetgrowth:pval.anetgrowth, round, digits = 3),
         chisq.anetgrowth = ifelse(chisq.anetgrowth < 0.001 & chisq.anetgrowth >= 0, 
                                "<0.001", chisq.anetgrowth),
         pval.anetgrowth = ifelse(pval.anetgrowth < 0.001 & pval.anetgrowth >= 0, 
                               "<0.001", pval.anetgrowth)) %>%
  arrange(treatment)

vcmax.coefs <- data.frame(summary(vcmax)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.vcmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.vcmax) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.vcmax) %>%
  print(., row.names = FALSE)

vcmax.table <- data.frame(Anova(vcmax)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(vcmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.vcmax, 
                chisq.vcmax = Chisq, pval.vcmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.vcmax:pval.vcmax, round, digits = 3),
         chisq.vcmax = ifelse(chisq.vcmax < 0.001 & chisq.vcmax >= 0, 
                              "<0.001", chisq.vcmax),
         pval.vcmax = ifelse(pval.vcmax < 0.001 & pval.vcmax >= 0, 
                             "<0.001", pval.vcmax)) %>%
  arrange(treatment)

jmax.coefs <- data.frame(summary(jmax)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.jmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.jmax) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.jmax) %>%
  print(., row.names = FALSE)

jmax.table <- data.frame(Anova(jmax)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(jmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.jmax, 
                chisq.jmax = Chisq, pval.jmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.jmax:pval.jmax, round, digits = 3),
         chisq.jmax = ifelse(chisq.jmax < 0.001 & chisq.jmax >= 0, 
                              "<0.001", chisq.jmax),
         pval.jmax = ifelse(pval.jmax < 0.001 & pval.jmax >= 0, 
                             "<0.001", pval.jmax)) %>%
  arrange(treatment)

jvmax.coefs <- data.frame(summary(jvmax)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.jvmax = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.jvmax) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.jvmax) %>%
  print(., row.names = FALSE)

jvmax.table <- data.frame(Anova(jvmax)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(jvmax.coefs) %>%
  dplyr::select(treatment, df = Df, coef.jvmax, 
                chisq.jvmax = Chisq, pval.jvmax = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.jvmax:pval.jvmax, round, digits = 3),
         chisq.jvmax = ifelse(chisq.jvmax < 0.001 & chisq.jvmax >= 0, 
                             "<0.001", chisq.jvmax),
         pval.jvmax = ifelse(pval.jvmax < 0.001 & pval.jvmax >= 0, 
                            "<0.001", pval.jvmax)) %>%
  arrange(treatment)


rd25.coefs <- data.frame(summary(rd25)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.rd25 = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.rd25) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.rd25) %>%
  print(., row.names = FALSE)

rd25.table <- data.frame(Anova(rd25)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(rd25.coefs) %>%
  dplyr::select(treatment, df = Df, coef.rd25, 
                chisq.rd25 = Chisq, pval.rd25 = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.rd25:pval.rd25, round, digits = 3),
         chisq.rd25 = ifelse(chisq.rd25 < 0.001 & chisq.rd25 >= 0, 
                            "<0.001", chisq.rd25),
         pval.rd25 = ifelse(pval.rd25 < 0.001 & pval.rd25 >= 0, 
                           "<0.001", pval.rd25)) %>%
  arrange(treatment)

table2 <- anet420.table %>% full_join(anetgrowth.table) %>% 
  full_join(vcmax.table) %>% full_join(jmax.table) %>%
  full_join(jvmax.table) %>% full_join(rd25.table) %>%
  replace(is.na(.), "-")

##########################################################################
## Table 3: PNUE/iWUE
##########################################################################
pnue.coefs <- data.frame(summary(pnue)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.pnue = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.pnue) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

pnue.table <- data.frame(Anova(pnue)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(pnue.coefs) %>%
  dplyr::select(treatment, df = Df, coef.pnue, 
                chisq.pnue = Chisq, pval.pnue = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.pnue:pval.pnue, round, digits = 3),
         chisq.pnue = ifelse(chisq.pnue < 0.001 & chisq.pnue >= 0, 
                              "<0.001", chisq.pnue),
         pval.pnue = ifelse(pval.pnue < 0.001 & pval.pnue >= 0, 
                             "<0.001", pval.pnue)) %>%
  arrange(treatment)

chi.coefs <- data.frame(summary(chi)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.chi = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.chi) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

chi.table <- data.frame(Anova(chi)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(chi.coefs) %>%
  dplyr::select(treatment, df = Df, coef.chi, 
                chisq.chi = Chisq, pval.chi = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.chi:pval.chi, round, digits = 3),
         chisq.chi = ifelse(chisq.chi < 0.001 & chisq.chi >= 0, 
                             "<0.001", chisq.chi),
         pval.chi = ifelse(pval.chi < 0.001 & pval.chi >= 0, 
                            "<0.001", pval.chi)) %>%
  arrange(treatment)

table3 <- pnue.table %>% full_join(chi.table) %>% replace(is.na(.), "-")

##########################################################################
## Table 4: Whole plant traits
##########################################################################
tla.coefs <- data.frame(summary(tla)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.tla = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.tla) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

tla.table <- data.frame(Anova(tla)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(tla.coefs) %>%
  dplyr::select(treatment, df = Df, coef.tla, 
                chisq.tla = Chisq, pval.tla = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.tla:pval.tla, round, digits = 3),
         chisq.tla = ifelse(chisq.tla <0.001 & chisq.tla >= 0, 
                            "<0.001", chisq.tla),
         pval.tla = ifelse(pval.tla <0.001 & pval.tla >= 0, 
                           "<0.001", pval.tla)) %>%
  arrange(treatment)

tbio.coefs <- data.frame(summary(tbio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.tbio = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.tbio) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

tbio.table <- data.frame(Anova(tbio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(tbio.coefs) %>%
  dplyr::select(treatment, df = Df, coef.tbio, 
                chisq.tbio = Chisq, pval.tbio = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.tbio:pval.tbio, round, digits = 3),
         chisq.tbio = ifelse(chisq.tbio < 0.001 & chisq.tbio >= 0, 
                             "<0.001", chisq.tbio),
         pval.tbio = ifelse(pval.tbio < 0.001 & pval.tbio >= 0, 
                            "<0.001", pval.tbio)) %>%
  arrange(treatment)

ncost.coefs <- data.frame(summary(ncost)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.ncost = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.ncost) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

ncost.table <- data.frame(Anova(ncost)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(ncost.coefs) %>%
  dplyr::select(treatment, df = Df, coef.ncost, 
                chisq.ncost = Chisq, pval.ncost = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.ncost:pval.ncost, round, digits = 3),
         chisq.ncost = ifelse(chisq.ncost <0.001 & chisq.ncost >= 0, 
                              "<0.001", chisq.ncost),
         pval.ncost = ifelse(pval.ncost <0.001 & pval.ncost >= 0, 
                             "<0.001", pval.ncost)) %>%
  arrange(treatment)

ndfa.coefs <- data.frame(summary(ndfa)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.ndfa = format(Estimate, scientific = TRUE, digits = 3),
         se.ndfa = round(Std..Error, digits = 3),
         t.value.ndfa = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.ndfa, se.ndfa, t.value.ndfa) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.ndfa) %>%
  print(., row.names = FALSE)

ndfa.table <- data.frame(Anova(ndfa)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(ndfa.coefs) %>%
  dplyr::select(treatment, df = Df, coef.ndfa,
                chisq.ndfa = Chisq, pval.ndfa = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment,
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.ndfa:pval.ndfa, round, digits = 3),
         chisq.ndfa = ifelse(chisq.ndfa < 0.001 & chisq.ndfa >= 0,
                             "<0.001", chisq.ndfa),
         pval.ndfa = ifelse(pval.ndfa < 0.001 & pval.ndfa >= 0,
                            "<0.001", pval.ndfa)) %>%
  arrange(treatment)

table4 <- tla.table %>% full_join(tbio.table) %>%
  full_join(ncost.table) %>% full_join(ndfa.table) %>%
  replace(is.na(.), "-")

##########################################################################
## Table S3: Components of Ncost
##########################################################################
cbg.coefs <- data.frame(summary(cbg)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.cbg = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.cbg) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

cbg.table <- data.frame(Anova(cbg)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(cbg.coefs) %>%
  dplyr::select(treatment, df = Df, coef.cbg, 
                chisq.cbg = Chisq, pval.cbg = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.cbg:pval.cbg, round, digits = 3),
         chisq.cbg = ifelse(chisq.cbg <0.001 & chisq.cbg >= 0, 
                            "<0.001", chisq.cbg),
         pval.cbg = ifelse(pval.cbg <0.001 & pval.cbg >= 0, 
                           "<0.001", pval.cbg)) %>%
  arrange(treatment)

wpn.coefs <- data.frame(summary(wpn)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.wpn = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.wpn) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

wpn.table <- data.frame(Anova(wpn)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(wpn.coefs) %>%
  dplyr::select(treatment, df = Df, coef.wpn, 
                chisq.wpn = Chisq, pval.wpn = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.wpn:pval.wpn, round, digits = 3),
         chisq.wpn = ifelse(chisq.wpn <0.001 & chisq.wpn >= 0, 
                            "<0.001", chisq.wpn),
         pval.wpn = ifelse(pval.wpn <0.001 & pval.wpn >= 0, 
                           "<0.001", pval.wpn)) %>%
  arrange(treatment)

tableS3 <- cbg.table %>% full_join(wpn.table) %>% replace(is.na(.), "-")

##########################################################################
## Table S4: Nitrogen fixation
##########################################################################
nodroot.coefs <- data.frame(summary(nod.root.ratio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.nodroot = format(Estimate, scientific = TRUE, digits = 3),
         se.nodroot = round(Std..Error, digits = 3),
         t.value.nodroot = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.nodroot, se.nodroot, t.value.nodroot) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.nodroot) %>%
  print(., row.names = FALSE)

nodroot.table <- data.frame(Anova(nod.root.ratio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(nodroot.coefs) %>%
  dplyr::select(treatment, df = Df, coef.nodroot, 
                chisq.nodroot = Chisq, pval.nodroot = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.nodroot:pval.nodroot, round, digits = 3),
         chisq.nodroot = ifelse(chisq.nodroot <0.001 & chisq.nodroot >= 0, 
                                "<0.001", chisq.nodroot),
         pval.nodroot = ifelse(pval.nodroot <0.001 & pval.nodroot >= 0, 
                               "<0.001", pval.nodroot)) %>%
  arrange(treatment)

nodbio.coefs <- data.frame(summary(nod.bio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.nodbio = format(Estimate, scientific = TRUE, digits = 3),
         se.nodbio = round(Std..Error, digits = 3),
         t.value.nodbio = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.nodbio, se.nodbio, t.value.nodbio) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.nodbio) %>%
  print(., row.names = FALSE)

nodbio.table <- data.frame(Anova(nod.bio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(nodbio.coefs) %>%
  dplyr::select(treatment, df = Df, coef.nodbio, 
                chisq.nodbio = Chisq, pval.nodbio = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.nodbio:pval.nodbio, round, digits = 3),
         chisq.nodbio = ifelse(chisq.nodbio <0.001 & chisq.nodbio >= 0, 
                               "<0.001", chisq.nodbio),
         pval.nodbio = ifelse(pval.nodbio <0.001 & pval.nodbio >= 0, 
                              "<0.001", pval.nodbio)) %>%
  arrange(treatment)

tableS4 <- nodbio.table %>% full_join(nodroot.table) %>% replace(is.na(.), "-")
