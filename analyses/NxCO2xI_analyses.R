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

# % change co2
(1.380 - 1.932) / 1.932 * 100

# % change inoc
(1.864 - 1.449) / 1.449 * 100

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

# % change CO2
(72.830 - 87.096) / 87.096 * 100

# % change inoculation
(94.700 - 65.225) / 65.225 * 100

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

# % change CO2
(137.364 - 151.834) / 151.834 * 100

# % change inoculation
(168.318 - 119.880) / 119.880 * 100

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

## Does inoculation increase positive effect of elevated CO2 on total biomass?
emmeans(tbio, pairwise~co2*inoc, at = list(n.ppm = 0)) # No

##########################################################################
## Leaf biomass
##########################################################################
leaf.bio <- lmer(sqrt(leaf.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(leaf.bio)
qqnorm(residuals(leaf.bio))
qqline(residuals(leaf.bio))
densityPlot(residuals(leaf.bio))
shapiro.test(residuals(leaf.bio))
outlierTest(leaf.bio)

# Model results
summary(leaf.bio)
Anova(leaf.bio)
r.squaredGLMM(leaf.bio)

# Pairwise comparisons
test(emtrends(leaf.bio, pairwise~co2, "n.trt"))
test(emtrends(leaf.bio, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(leaf.bio, pairwise~co2, type = "response")
emmeans(leaf.bio, pairwise~inoc)
test(emtrends(leaf.bio, ~1, "n.trt"))

##########################################################################
## Leaf biomass
##########################################################################
stem.bio <- lmer(sqrt(stem.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(stem.bio)
qqnorm(residuals(stem.bio))
qqline(residuals(stem.bio))
densityPlot(residuals(stem.bio))
shapiro.test(residuals(stem.bio))
outlierTest(stem.bio)

# Model results
summary(stem.bio)
Anova(stem.bio)
r.squaredGLMM(stem.bio)

# Pairwise comparisons
emmeans(stem.bio, pairwise~co2, type = "response")


test(emtrends(stem.bio, pairwise~co2, "n.trt"))
test(emtrends(stem.bio, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(stem.bio, pairwise~co2, type = "response")
emmeans(stem.bio, pairwise~inoc)
test(emtrends(stem.bio, ~1, "n.trt"))

##########################################################################
## Root biomass
##########################################################################
root.bio <- lmer(sqrt(root.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

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
emmeans(root.bio, pairwise~co2*inoc, type = "response")
test(emtrends(root.bio, pairwise~co2, "n.trt"))
test(emtrends(root.bio, pairwise~inoc, "n.trt"))

# Individual effects
emmeans(root.bio, pairwise~co2, type = "response")
emmeans(root.bio, pairwise~inoc)
test(emtrends(root.bio, ~1, "n.trt"))

##########################################################################
## Root nodule biomass
##########################################################################
df$nodule.biomass[80] <- NA

nod.bio <- lmer(sqrt(nodule.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

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
test(emtrends(nod.bio, pairwise~co2, "n.trt"))
test(emtrends(nod.bio, pairwise~inoc, "n.trt"))

## Individual effects
emmeans(nod.bio, pairwise~co2)
emmeans(nod.bio, pairwise~inoc)
test(emtrends(nod.bio, ~1, "n.trt"))

##########################################################################
## Root:shoot ratio
##########################################################################
df$root.shoot.ratio[101] <- NA

rootshoot <- lmer(log(root.shoot.ratio) ~ co2 * inoc * n.trt + (1|rack:co2), 
                  data = df)

# Check model fit
plot(rootshoot)
qqnorm(residuals(rootshoot))
qqline(residuals(rootshoot))
densityPlot(residuals(rootshoot))
shapiro.test(residuals(rootshoot))
outlierTest(rootshoot)

# Model results
summary(rootshoot)
Anova(rootshoot)
r.squaredGLMM(rootshoot)

# Pairwise comparisons
emmeans(rootshoot, pairwise~co2, type = "response")
cld(emmeans(rootshoot, ~co2*inoc, type = "response"))

test(emtrends(rootshoot, ~1, "n.trt"))
test(emtrends(rootshoot, ~inoc, "n.trt"))

##########################################################################
## Leaf mass fraction
##########################################################################
df$lmf[c(37, 38, 101, 102, 104)] <- NA

lmf <- lmer(lmf ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(lmf)
qqnorm(residuals(lmf))
qqline(residuals(lmf))
densityPlot(residuals(lmf))
shapiro.test(residuals(lmf))
outlierTest(lmf)

# Model results
summary(lmf)
Anova(lmf)
r.squaredGLMM(lmf)

# Pairwise comparisons
cld(emmeans(lmf, pairwise~co2*inoc))
test(emtrends(lmf, pairwise~inoc, "n.trt"))

emmeans(lmf, pairwise~co2)
emmeans(lmf, pairwise~inoc)
test(emtrends(lmf, ~1, "n.trt"))

# % change CO2
(0.398 - 0.369) / 0.369 * 100

##########################################################################
## Stem mass fraction
##########################################################################
df$smf[113] <- NA

smf <- lmer(smf ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(smf)
qqnorm(residuals(smf))
qqline(residuals(smf))
densityPlot(residuals(smf))
shapiro.test(residuals(smf))
outlierTest(smf)

# Model results
summary(smf)
Anova(smf)
r.squaredGLMM(smf)

# Pairwise comparisons
test(emtrends(smf, ~1, "n.trt"))
test(emtrends(smf, ~inoc, "n.trt"))

##########################################################################
## Root mass fraction
##########################################################################
df$rmf[c(101)] <- NA

rmf <- lmer(log(rmf) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(rmf)
qqnorm(residuals(rmf))
qqline(residuals(rmf))
densityPlot(residuals(rmf))
shapiro.test(residuals(rmf))
outlierTest(rmf)

# Model results
summary(rmf)
Anova(rmf)
r.squaredGLMM(rmf)

# Pairwise comparisons
cld(emmeans(rmf, ~co2*inoc, type = "response"))
test(emtrends(rmf, ~inoc, "n.trt"))

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
cld(emtrends(ncost, pairwise~co2*inoc, "n.trt"))
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
cbg <- lmer(cbg ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

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
nod.root.ratio <- lmer(sqrt(nod.root.ratio) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

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
cld(emmeans(nod.root.ratio, pairwise~co2*inoc, type = "response"))

## Individual effects
emmeans(nod.root.ratio, pairwise~co2)
emmeans(nod.root.ratio, pairwise~inoc)
test(emtrends(nod.root.ratio, ~1, "n.trt"))

# Percent change in inoculated pots
emmeans(nod.root.ratio, ~inoc, "n.trt", at = list(n.trt = c(0, 630)))

##########################################################################
## Biomass : pot volume
##########################################################################
bvr <- lmer(bvr ~ co2 * inoc * n.trt + (1|rack:co2), data = df)

# Check model fit
plot(bvr)
qqnorm(residuals(bvr))
qqline(residuals(bvr))
densityPlot(residuals(bvr))
shapiro.test(residuals(bvr))
outlierTest(bvr)

# Model results
summary(bvr)
Anova(bvr)
r.squaredGLMM(bvr)

# Pairwise comparisons
test(emtrends(bvr, pairwise~inoc, "n.trt"))
test(emtrends(bvr, pairwise~co2, "n.trt"))
test(emtrends(bvr, ~1, "n.trt"))
emmeans(bvr, pairwise~inoc)

##########################################################################
## Table 1: Leaf N content
##########################################################################
# Narea
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
         across(chisq.narea:pval.narea, \(x) round(x, digits = 3)),
         chisq.narea = ifelse(chisq.narea <0.001 & chisq.narea >= 0, 
                              "<0.001", chisq.narea),
         pval.narea = ifelse(pval.narea <0.001 & pval.narea >= 0, 
                             "<0.001", pval.narea)) %>%
  arrange(treatment)

# Chlarea
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
         across(chisq.chl.area:pval.chl.area, \(x) round(x, digits = 3)),
         chisq.chl.area = ifelse(chisq.chl.area <0.001 & chisq.chl.area >= 0, 
                              "<0.001", chisq.chl.area),
         pval.chl.area = ifelse(pval.chl.area <0.001 & pval.chl.area >= 0, 
                             "<0.001", pval.chl.area)) %>%
  arrange(treatment)

table1 <- narea.table %>% full_join(chl.area.table) %>%
  replace(is.na(.), "-")

##########################################################################
## Table 2: Gas exchange
##########################################################################
# Anet,420
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
         across(chisq.anet420:pval.anet420, \(x) round(x, digits = 3)),
         chisq.anet420 = ifelse(chisq.anet420 < 0.001 & chisq.anet420 >= 0, 
                              "<0.001", chisq.anet420),
         pval.anet420 = ifelse(pval.anet420 < 0.001 & pval.anet420 >= 0, 
                             "<0.001", pval.anet420)) %>%
  arrange(treatment)

# Anet,gc
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
         across(chisq.anetgrowth:pval.anetgrowth, \(x) round(x, digits = 3)),
         chisq.anetgrowth = ifelse(chisq.anetgrowth < 0.001 & chisq.anetgrowth >= 0, 
                                "<0.001", chisq.anetgrowth),
         pval.anetgrowth = ifelse(pval.anetgrowth < 0.001 & pval.anetgrowth >= 0, 
                               "<0.001", pval.anetgrowth)) %>%
  arrange(treatment)

# Vcmax25
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
         across(chisq.vcmax:pval.vcmax, \(x) round(x, digits = 3)),
         chisq.vcmax = ifelse(chisq.vcmax < 0.001 & chisq.vcmax >= 0, 
                              "<0.001", chisq.vcmax),
         pval.vcmax = ifelse(pval.vcmax < 0.001 & pval.vcmax >= 0, 
                             "<0.001", pval.vcmax)) %>%
  arrange(treatment)

# Jmax25
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
         across(chisq.jmax:pval.jmax, \(x) round(x, digits = 3)),
         chisq.jmax = ifelse(chisq.jmax < 0.001 & chisq.jmax >= 0, 
                              "<0.001", chisq.jmax),
         pval.jmax = ifelse(pval.jmax < 0.001 & pval.jmax >= 0, 
                             "<0.001", pval.jmax)) %>%
  arrange(treatment)

# Jmax25: Vcmax25
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
         across(chisq.jvmax:pval.jvmax, \(x) round(x, digits = 3)),
         chisq.jvmax = ifelse(chisq.jvmax < 0.001 & chisq.jvmax >= 0, 
                             "<0.001", chisq.jvmax),
         pval.jvmax = ifelse(pval.jvmax < 0.001 & pval.jvmax >= 0, 
                            "<0.001", pval.jvmax)) %>%
  arrange(treatment)

table2 <- anet420.table %>% full_join(anetgrowth.table) %>% 
  full_join(vcmax.table) %>% full_join(jmax.table) %>%
  full_join(jvmax.table) %>% replace(is.na(.), "-")

##########################################################################
## Table 3: Whole plant traits
##########################################################################
# Total leaf area
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
         across(chisq.tla:pval.tla, \(x) round(x, digits = 3)),
         chisq.tla = ifelse(chisq.tla <0.001 & chisq.tla >= 0, 
                            "<0.001", chisq.tla),
         pval.tla = ifelse(pval.tla <0.001 & pval.tla >= 0, 
                           "<0.001", pval.tla)) %>%
  arrange(treatment)

# Total biomass
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
         across(chisq.tbio:pval.tbio, \(x) round(x, digits = 3)),
         chisq.tbio = ifelse(chisq.tbio < 0.001 & chisq.tbio >= 0, 
                             "<0.001", chisq.tbio),
         pval.tbio = ifelse(pval.tbio < 0.001 & pval.tbio >= 0, 
                            "<0.001", pval.tbio)) %>%
  arrange(treatment)

# Leaf area ratio
lar.coefs <- data.frame(summary(lar)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.lar = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.lar) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.lar) %>%
  print(., row.names = FALSE)

lar.table <- data.frame(Anova(lar)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(lar.coefs) %>%
  dplyr::select(treatment, df = Df, coef.lar, 
                chisq.lar = Chisq, pval.lar = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.lar:pval.lar, \(x) round(x, digits = 3)),
         chisq.lar = ifelse(chisq.lar < 0.001 & chisq.lar >= 0, 
                            "<0.001", chisq.lar),
         pval.lar = ifelse(pval.lar < 0.001 & pval.lar >= 0, 
                           "<0.001", pval.lar)) %>%
  arrange(treatment)

# Carbon cost to acquire nitrogen
rootshoot.coefs <- data.frame(summary(rootshoot)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.rootshoot = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.rootshoot) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  print(., row.names = FALSE)

rootshoot.table <- data.frame(Anova(rootshoot)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(rootshoot.coefs) %>%
  dplyr::select(treatment, df = Df, coef.rootshoot, 
                chisq.rootshoot = Chisq, pval.rootshoot = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.rootshoot:pval.rootshoot, \(x) round(x, digits = 3)),
         chisq.rootshoot = ifelse(chisq.rootshoot <0.001 & chisq.rootshoot >= 0, 
                              "<0.001", chisq.rootshoot),
         pval.rootshoot = ifelse(pval.rootshoot <0.001 & pval.rootshoot >= 0, 
                             "<0.001", pval.rootshoot)) %>%
  arrange(treatment)

# Carbon cost to acquire nitrogen
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
         across(chisq.ncost:pval.ncost, \(x) round(x, digits = 3)),
         chisq.ncost = ifelse(chisq.ncost <0.001 & chisq.ncost >= 0, 
                              "<0.001", chisq.ncost),
         pval.ncost = ifelse(pval.ncost <0.001 & pval.ncost >= 0, 
                             "<0.001", pval.ncost)) %>%
  arrange(treatment)

table3 <- tla.table %>% full_join(tbio.table) %>%
  full_join(lar.table) %>%
  full_join(rootshoot.table) %>%
  full_join(ncost.table) %>%
  replace(is.na(.), "-")

##########################################################################
## Table S3: Nmass, LMA
##########################################################################
# Nmass
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
         across(chisq.nmass:pval.nmass, \(x) round(x, digits = 3)),
         chisq.nmass = ifelse(chisq.nmass <0.001 & chisq.nmass >= 0, 
                              "<0.001", chisq.nmass),
         pval.nmass = ifelse(pval.nmass <0.001 & pval.nmass >= 0, 
                             "<0.001", pval.nmass)) %>%
  arrange(treatment)

# Marea
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
         across(chisq.marea:pval.marea, \(x) round(x, digits = 3)),
         chisq.marea = ifelse(chisq.marea <0.001 & chisq.marea >= 0, 
                              "<0.001", chisq.marea),
         pval.marea = ifelse(pval.marea <0.001 & pval.marea >= 0, 
                             "<0.001", pval.marea)) %>%
  arrange(treatment)

tableS3 <- nmass.table %>% full_join(marea.table) %>%
  replace(is.na(.), "-")

##########################################################################
## Table S4: Dark respiration, PNUE
##########################################################################
# Dark respiration
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
         across(chisq.rd25:pval.rd25, \(x) round(x, digits = 3)),
         chisq.rd25 = ifelse(chisq.rd25 < 0.001 & chisq.rd25 >= 0, 
                             "<0.001", chisq.rd25),
         pval.rd25 = ifelse(pval.rd25 < 0.001 & pval.rd25 >= 0, 
                            "<0.001", pval.rd25)) %>%
  arrange(treatment)

# PNUEgc
pnue.coefs <- data.frame(summary(pnue)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.pnue = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.pnue) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.pnue) %>%
  print(., row.names = FALSE)

pnue.table <- data.frame(Anova(pnue)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(pnue.coefs) %>%
  dplyr::select(treatment, df = Df, coef.pnue, 
                chisq.pnue = Chisq, pval.pnue = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.pnue:pval.pnue, \(x) round(x, digits = 3)),
         chisq.pnue = ifelse(chisq.pnue < 0.001 & chisq.pnue >= 0, 
                             "<0.001", chisq.pnue),
         pval.pnue = ifelse(pval.pnue < 0.001 & pval.pnue >= 0, 
                            "<0.001", pval.pnue)) %>%
  arrange(treatment)

tableS4 <- rd25.table %>% full_join(pnue.table) %>%
  replace(is.na(.), "-")

##########################################################################
## Table S5: Biomass partitioning
##########################################################################
# Leaf biomass
leafbio.coefs <- data.frame(summary(leaf.bio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.leafbio = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.leafbio) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.leafbio) %>%
  print(., row.names = FALSE)

leafbio.table <- data.frame(Anova(leaf.bio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(leafbio.coefs) %>%
  dplyr::select(treatment, df = Df, coef.leafbio, 
                chisq.leafbio = Chisq, pval.leafbio = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.leafbio:pval.leafbio, \(x) round(x, digits = 3)),
         chisq.leafbio = ifelse(chisq.leafbio < 0.001 & chisq.leafbio >= 0, 
                            "<0.001", chisq.leafbio),
         pval.leafbio = ifelse(pval.leafbio < 0.001 & pval.leafbio >= 0, 
                           "<0.001", pval.leafbio)) %>%
  arrange(treatment)

# Stem biomass
stembio.coefs <- data.frame(summary(stem.bio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.stembio = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.stembio) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.stembio) %>%
  print(., row.names = FALSE)

stembio.table <- data.frame(Anova(stem.bio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(stembio.coefs) %>%
  dplyr::select(treatment, df = Df, coef.stembio, 
                chisq.stembio = Chisq, pval.stembio = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.stembio:pval.stembio, \(x) round(x, digits = 3)),
         chisq.stembio = ifelse(chisq.stembio < 0.001 & chisq.stembio >= 0, 
                                "<0.001", chisq.stembio),
         pval.stembio = ifelse(pval.stembio < 0.001 & pval.stembio >= 0, 
                               "<0.001", pval.stembio)) %>%
  arrange(treatment)

# Root biomass 
rootbio.coefs <- data.frame(summary(root.bio)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.rootbio = format(Estimate, scientific = TRUE, digits = 3),
         se.rootbio = round(Std..Error, digits = 3),
         t.value.rootbio = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.rootbio, se.rootbio, t.value.rootbio) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.rootbio) %>%
  print(., row.names = FALSE)

rootbio.table <- data.frame(Anova(root.bio)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(rootbio.coefs) %>%
  dplyr::select(treatment, df = Df, coef.rootbio, 
                chisq.rootbio = Chisq, pval.rootbio = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.rootbio:pval.rootbio, \(x) round(x, digits = 3)),
         chisq.rootbio = ifelse(chisq.rootbio <0.001 & chisq.rootbio >= 0, 
                                "<0.001", chisq.rootbio),
         pval.rootbio = ifelse(pval.rootbio <0.001 & pval.rootbio >= 0, 
                               "<0.001", pval.rootbio)) %>%
  arrange(treatment)

# Root nodule biomass
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
         across(chisq.nodbio:pval.nodbio, \(x) round(x, digits = 3)),
         chisq.nodbio = ifelse(chisq.nodbio <0.001 & chisq.nodbio >= 0, 
                               "<0.001", chisq.nodbio),
         pval.nodbio = ifelse(pval.nodbio <0.001 & pval.nodbio >= 0, 
                              "<0.001", pval.nodbio)) %>%
  arrange(treatment)

# Leaf mass fraction
lmf.coefs <- data.frame(summary(lmf)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.lmf = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.lmf) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.lmf) %>%
  print(., row.names = FALSE)

lmf.table <- data.frame(Anova(lmf)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(lmf.coefs) %>%
  dplyr::select(treatment, df = Df, coef.lmf, 
                chisq.lmf = Chisq, pval.lmf = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.lmf:pval.lmf, \(x) round(x, digits = 3)),
         chisq.lmf = ifelse(chisq.lmf < 0.001 & chisq.lmf >= 0, 
                            "<0.001", chisq.lmf),
         pval.lmf = ifelse(pval.lmf < 0.001 & pval.lmf >= 0, 
                           "<0.001", pval.lmf)) %>%
  arrange(treatment)

# Stem mass fraction
smf.coefs <- data.frame(summary(smf)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.smf = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.smf) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.smf) %>%
  print(., row.names = FALSE)

smf.table <- data.frame(Anova(smf)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(smf.coefs) %>%
  dplyr::select(treatment, df = Df, coef.smf, 
                chisq.smf = Chisq, pval.smf = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.smf:pval.smf, \(x) round(x, digits = 3)),
         chisq.smf = ifelse(chisq.smf < 0.001 & chisq.smf >= 0, 
                            "<0.001", chisq.smf),
         pval.smf = ifelse(pval.smf < 0.001 & pval.smf >= 0, 
                           "<0.001", pval.smf)) %>%
  arrange(treatment)

# Root mass fraction
rmf.coefs <- data.frame(summary(rmf)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.rmf = format(Estimate, scientific = TRUE, digits = 3)) %>%
  dplyr::select(treatment, coef.rmf) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.rmf) %>%
  print(., row.names = FALSE)

rmf.table <- data.frame(Anova(rmf)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(rmf.coefs) %>%
  dplyr::select(treatment, df = Df, coef.rmf, 
                chisq.rmf = Chisq, pval.rmf = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.rmf:pval.rmf, \(x) round(x, digits = 3)),
         chisq.rmf = ifelse(chisq.rmf < 0.001 & chisq.rmf >= 0, 
                            "<0.001", chisq.rmf),
         pval.rmf = ifelse(pval.rmf < 0.001 & pval.rmf >= 0, 
                           "<0.001", pval.rmf)) %>%
  arrange(treatment)

# Root nodule: root biomass
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
         across(chisq.nodroot:pval.nodroot, \(x) round(x, digits = 3)),
         chisq.nodroot = ifelse(chisq.nodroot <0.001 & chisq.nodroot >= 0, 
                                "<0.001", chisq.nodroot),
         pval.nodroot = ifelse(pval.nodroot <0.001 & pval.nodroot >= 0, 
                               "<0.001", pval.nodroot)) %>%
  arrange(treatment)

tableS5 <- leafbio.table %>%
  full_join(stembio.table) %>%
  full_join(rootbio.table) %>%
  full_join(nodbio.table) %>%
  full_join(lmf.table) %>%
  full_join(smf.table) %>%
  full_join(rmf.table) %>%
  full_join(nodroot.table) %>% replace(is.na(.), "-")

##########################################################################
## Table S6: Components of Ncost
##########################################################################
# Belowground carbon allocation
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
         across(chisq.cbg:pval.cbg, \(x) round(x, digits = 3)),
         chisq.cbg = ifelse(chisq.cbg <0.001 & chisq.cbg >= 0, 
                            "<0.001", chisq.cbg),
         pval.cbg = ifelse(pval.cbg <0.001 & pval.cbg >= 0, 
                           "<0.001", pval.cbg)) %>%
  arrange(treatment)

# Whole-plant nitrogen biomass
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
         across(chisq.wpn:pval.wpn, \(x) round(x, digits = 3)),
         chisq.wpn = ifelse(chisq.wpn <0.001 & chisq.wpn >= 0, 
                            "<0.001", chisq.wpn),
         pval.wpn = ifelse(pval.wpn <0.001 & pval.wpn >= 0, 
                           "<0.001", pval.wpn)) %>%
  arrange(treatment)

tableS6 <- cbg.table %>% full_join(wpn.table) %>% replace(is.na(.), "-")

##########################################################################
## Table S7: BVR
##########################################################################
# Biomass: pot volume
bvr.coefs <- data.frame(summary(bvr)$coefficient) %>%
  mutate(treatment = row.names(.),
         coef.bvr = format(Estimate, scientific = TRUE, digits = 3),
         se.bvr = round(Std..Error, digits = 3),
         t.value.bvr = round(t.value, digits = 3)) %>%
  dplyr::select(treatment, coef.bvr, se.bvr, t.value.bvr) %>%
  mutate(treatment = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")) %>%
  dplyr::select(treatment, coef.bvr) %>%
  print(., row.names = FALSE)

bvr.table <- data.frame(Anova(bvr)) %>%
  mutate(treatment = row.names(.)) %>%
  full_join(bvr.coefs) %>%
  dplyr::select(treatment, df = Df, coef.bvr, 
                chisq.bvr = Chisq, pval.bvr = Pr..Chisq.) %>%
  mutate(treatment = factor(treatment, 
                            levels = c("(Intercept)", "co2", "inoc", "n.trt", "co2:inoc",
                                       "co2:n.trt", "inoc:n.trt", "co2:inoc:n.trt")),
         across(chisq.bvr:pval.bvr, \(x) round(x, digits = 3)),
         chisq.bvr = ifelse(chisq.bvr <0.001 & chisq.bvr >= 0, 
                                "<0.001", chisq.bvr),
         pval.bvr = ifelse(pval.bvr <0.001 & pval.bvr >= 0, 
                               "<0.001", pval.bvr)) %>%
  arrange(treatment)
