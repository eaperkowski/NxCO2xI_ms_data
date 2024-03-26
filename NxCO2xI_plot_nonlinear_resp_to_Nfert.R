# Code for nitrogen WG manuscript that explores and plots nonlinear effects of 
# increasing nitrogen fertilization on indices of photosynthetic capacity
# (Vcmax25, Jmax25)
#
# Note that this script assumes that the eaperkowski/NxCO2xI_ms_data data
# repository is the root working directory
# URL: https://github.com/eaperkowski/NxCO2xI_ms_data

##########################################################################
# Load libraries and import data
##########################################################################
# Libraries
library(tidyverse)
library(nlme)
library(ggpubr)
library(emmeans)

# Read compiled data file, ensure N fertilization is numeric, then filter 
# all uninoculated plants where nodule biomass:root biomass was greater than 
# 0.05 g/g
df <- read.csv("data/NxCO2xI_data.csv") %>%
  mutate(n.trt = as.numeric(n.trt)) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nod.root.ratio < 0.05))

##########################################################################
# Vcmax25 nonlinear saturating regression model and plot
##########################################################################

# Create nonlinear saturating regression to explain nonlinear saturating
# Vcmax response to fertilization. Note that model uses data collected from 
# only uninoculated individuals grown under elevated CO2
vcmax.nls.elv <- nls(formula = vcmax25 ~ a + ((b * n.trt) / (c + n.trt)),
                     start = list(a = 24, b = 114, c = 454),
                     data = subset(df, inoc == "no.inoc" & co2 == "elv"))

# Create nonlinear saturating regression to explain nonlinear saturating
# Vcmax response to fertilization. Note that model uses data collected from 
# only uninoculated individuals grown under elevated CO2
vcmax.nls.amb <- nls(formula = vcmax25 ~ a + ((b * n.trt) / (c + n.trt)),
                     start = list(a = 3, b = 127, c = 114),
                     data = subset(df, inoc == "no.inoc" & co2 == "amb"))

# Create predicted trendlines along range in N fertilization values for 
# both CO2 treatments
vcmax.nls.elv.pred <- data.frame(
  emmeans(vcmax.nls.elv, ~1, "n.trt",
          at = list(n.trt = seq(0, 630, 1)),
          data = subset(df, inoc == "no.inoc" & co2 == "elv")))

vcmax.nls.amb.pred <- data.frame(
  emmeans(vcmax.nls.amb, ~1, "n.trt",
          at = list(n.trt = seq(0, 630, 1)),
          data = subset(df, inoc == "no.inoc" & co2 == "amb")))

# Create Vcmax plot
vcmax25.plot <- ggplot(data = subset(df, inoc == "no.inoc"), 
                       aes(x = n.trt, y = vcmax25)) +
  geom_point(aes(fill = co2), alpha = 0.75, size = 4, shape = 21) +
  geom_smooth(data = vcmax.nls.elv.pred, aes(y = emmean), color = "#b2182b",
              linewidth = 2, se = FALSE) +
  geom_smooth(data = vcmax.nls.amb.pred, aes(y = emmean), color = "#2166ac",
              linewidth = 2, se = FALSE) +
  geom_text(aes(500, 5, label=(paste(expression(" y = "*frac("114.85 * x", "454.90 + x")*" + 24.07")))),
                parse = TRUE, size = 3.5, color = "#b2182b") +
  geom_text(aes(500, 30, label=(paste(expression("y = "*frac("127.47 * x", "114.01 + x")*" +  3.29")))),
            parse = TRUE, size = 3.5, color = "#2166ac") +
  scale_fill_manual(values = c("#2166ac", "#b2182b"),
                     labels = c(expression("Ambient CO"["2"]),
                                expression("Elevated CO"["2"]))) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bold("CO"["2"]*" treatment"))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0)
vcmax25.plot

##########################################################################
## Jmax25 nonlinear saturating regression model and plot
##########################################################################

# Create nonlinear saturating regression to explain nonlinear saturating
# Jmax response to fertilization. Note that model uses data collected from 
# only uninoculated individuals grown under elevated CO2
jmax.nls.elv <- nls(formula = jmax25 ~ a + ((b * n.trt) / (c + n.trt)),
                    start = list(a = 46, b = 185, c = 347),
                    data = subset(df, inoc == "no.inoc" & co2 == "elv"))

# Create nonlinear saturating regression to explain nonlinear saturating
# Jmax response to fertilization. Note that model uses data collected from 
# only uninoculated individuals grown under ambient CO2
jmax.nls.amb <- nls(formula = jmax25 ~ a + ((b * n.trt) / (c + n.trt)),
                    start = list(a = 9, b = 207, c = 96),
                    data = subset(df, inoc == "no.inoc" & co2 == "amb"))

# Create predicted trendlines along range in N fertilization values for 
# both CO2 treatments
jmax.nls.elv.pred <- data.frame(
  emmeans(jmax.nls.elv, ~1, "n.trt",
          at = list(n.trt = seq(0, 630, 1)),
          data = subset(df, inoc == "no.inoc" & co2 == "elv")))

jmax.nls.amb.pred <- data.frame(
  emmeans(jmax.nls.amb, ~1, "n.trt",
          at = list(n.trt = seq(0, 630, 1)),
          data = subset(df, inoc == "no.inoc" & co2 == "amb")))

# Create Jmax25 plot
jmax25.plot <- ggplot(data = subset(df, inoc == "no.inoc"), 
                       aes(x = n.trt, y = jmax25)) +
  geom_point(aes(fill = co2), alpha = 0.75, size = 4, shape = 21) +
  geom_smooth(data = jmax.nls.elv.pred, aes(y = emmean), color = "#b2182b",
              linewidth = 2, se = FALSE) +
  geom_smooth(data = jmax.nls.amb.pred, aes(y = emmean), color = "#2166ac",
              linewidth = 2, se = FALSE) +
  geom_text(aes(500, 7.5, label=(paste(expression(" y = "*frac("185.98 * x", "347.68 + x")*" + 46.28")))),
            parse = TRUE, size = 3.5, color = "#b2182b") +
  geom_text(aes(500, 45, label=(paste(expression("y = "*frac("207.27 * x", "96.47 + x")*"  +  9.81")))),
            parse = TRUE, size = 3.5, color = "#2166ac") +
  scale_fill_manual(values = c("#2166ac", "#b2182b"),
                    labels = c(expression("Ambient CO"["2"]),
                               expression("Elevated CO"["2"]))) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bold("CO"["2"]*" treatment"))) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0)
jmax25.plot

##########################################################################
# Create paneled figure with Vcmax and Jmax nonlinear saturating 
# responses to N fertilization. Write file as .png with high resolution
##########################################################################
# png("[insert path here]", width = 12, 
#     height = 4.5, units = "in", res = 600)
ggarrange(vcmax25.plot, jmax25.plot, common.legend = TRUE, legend = "right",
          labels = c("a", "b"), font.label = list(size = 18),
          label.x = 0.2, label.y = 0.975)
# dev.off()
