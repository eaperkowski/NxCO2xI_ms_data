## NxCO2xI plot script. Paths assume the root directory of this script.

##########################################################################
## Load libraries and import data
##########################################################################
# Libraries
library(lme4)
library(emmeans)
library(tidyverse)
library(ggpubr)
library(multcomp)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

# Load compiled datasheet
df <- read.csv("../data/NxCO2xI_data.csv", 
               na.strings = "NA") %>%
  mutate(n.trt = as.numeric(n.trt),
         inoc = factor(inoc, levels = c("no.inoc", "inoc")),
         co2 = factor(co2, levels = c("amb", "elv"))) %>%
  filter(inoc == "inoc" | (inoc == "no.inoc" & nod.root.ratio < 0.05)) %>%
  unite(col = "co2.inoc", co2:inoc, sep = "_", remove = FALSE) %>%
  mutate(co2.inoc = factor(co2.inoc,
                           levels = c("elv_inoc", "elv_no.inoc",
                                      "amb_inoc", "amb_no.inoc")))

## Add colorblind friendly palette
co2.cols <- c("#2166ac", "#b2182b")
co2.cols <- c("#2166ac", "#b2182b")
full.cols <- c("#b2182b", "#f4a582", "#2166ac", "#92c5d3")

## Create blank plot as spacer plot
blank.plot <- ggplot() + 
  theme_bw() +
  theme(panel.background = element_rect(color = "white",
                                        fill = "white"),
        panel.border = element_rect(color = "white"))

##########################################################################
## Narea regression line prep
##########################################################################
narea <- lmer(narea ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(narea, ~inoc*co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
narea.regline <- data.frame(emmeans(narea, ~co2*inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Narea plot
##########################################################################
narea.plot <- ggplot(data = df, aes(x = n.trt, y = narea, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = narea.regline,
              aes(color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE, method = "loess") +
  geom_ribbon(data = narea.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 3.24), breaks = seq(0, 3.2, 0.8)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*" (gN m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
narea.plot

##########################################################################
## Narea regression line prep co2-by-inoculation interaction plot 
## in supplement
##########################################################################
test(emtrends(narea, pairwise~co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
narea.int.regline <- data.frame(emmeans(narea, ~co2, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response"))

##########################################################################
## Narea co2-by-inoculation interaction plot in supplement
##########################################################################
narea.int.plot <- ggplot(data = df, aes(x = n.trt, y = narea, fill = co2)) +
  geom_point(aes(shape = inoc), size = 3, alpha = 0.75) +
  geom_smooth(data = narea.int.regline,
              aes(color = co2, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = narea.int.regline,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient", "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient", "Elevated")) +
  scale_shape_manual(values = c(21, 24), labels = c("Uninoculated",
                                                    "Inoculated")) +
  scale_y_continuous(limits = c(0, 3.24), breaks = seq(0, 3.2, 0.8)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["area"]*" (gN m"^"-2"*")")),
       fill = expression(bold("CO"["2"])),
       color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))
narea.int.plot

##########################################################################
## Nmass regression line prep
##########################################################################
nmass <- lmer(nmass.focal ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nmass, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nmass.regline <- data.frame(emmeans(nmass, ~co2*inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "elv_inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, levels = c("elv_inoc","elv_no.inoc", 
                                                "amb_inoc", "amb_no.inoc")))

##########################################################################
## Nmass plot
##########################################################################
nmass.plot <- ggplot(data = df,
                     aes(x = n.trt, y = nmass.focal, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = nmass.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = nmass.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["mass"]*" (gN g"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
nmass.plot

##########################################################################
## Nmass regression line prep for co2-by-inoculation interaction plot 
## in supplement
##########################################################################
test(emtrends(nmass, pairwise~co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nmass.int.regline <- data.frame(emmeans(nmass, ~co2, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response"))

##########################################################################
## Nmass co2-by-inoculation interaction plot in supplement
##########################################################################
nmass.int.plot <- ggplot(data = df,
                         aes(x = n.trt, y = nmass.focal, fill = co2)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = nmass.int.regline,
              aes(color = co2, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = nmass.int.regline,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient", "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient", "Elevated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["mass"]*" (gN g"^"-1"*")")),
       fill = expression(bold("CO"["2"])),
       color = expression(bold("CO"["2"]))) +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))
nmass.int.plot

##########################################################################
## Marea regression line prep
##########################################################################
marea <- lmer(log(marea) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(marea, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
marea.regline <- data.frame(emmeans(marea, ~co2*inoc, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Marea plot
##########################################################################
marea.plot <- ggplot(data = df, 
                         aes(x = n.trt, y = marea, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = marea.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = marea.regline,
              aes(fill = co2.inoc, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, 15)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("M")["area"]*" (g m"^"-2"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
marea.plot

##########################################################################
## Marea regression line prep for co2-by-inoculation interaction plot 
## in supplement
##########################################################################
test(emtrends(marea, pairwise~co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
marea.int.regline <- data.frame(emmeans(marea, ~co2, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response"))

##########################################################################
## Marea co2-by-inoculation interaction plot in supplement
##########################################################################
marea.int.plot <- ggplot(data = df, aes(x = n.trt, y = marea, fill = co2)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = marea.int.regline,
              aes(color = co2, y = response), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = marea.int.regline,
              aes(fill = co2, y = response, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient", "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient", "Elevated")) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, 15)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("M")["area"]*" (g m"^"-2"*")")),
       fill = expression(bold("CO"["2"])),
       color = expression(bold("CO"["2"]))) +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))
marea.int.plot

##########################################################################
## Chlarea regression line prep
##########################################################################
chlarea <- lmer(chl.mmolm2 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(chlarea, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
chlarea.regline <- data.frame(emmeans(chlarea, ~co2*inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Chlarea plot
##########################################################################
chlarea.plot <- ggplot(data = df,
                       aes(x = n.trt, y = chl.mmolm2, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = chlarea.regline,
              aes(color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = chlarea.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.07)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("Chl")["area"]*" (mmol m"^"-2"*")")),
       fill = "Treatment", color = "Treatment", shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
chlarea.plot

##########################################################################
## Anet,420 regression line prep
##########################################################################
anet <- lmer(anet ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(anet, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
anet.regline <- data.frame(emmeans(anet, ~co2*inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Anet,420 plot
##########################################################################
anet.plot <- ggplot(data = df, 
                       aes(x = n.trt, y = anet, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = anet.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = anet.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 32), breaks = seq(0, 32, 8)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("A")["net,420"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment", shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
anet.plot

##########################################################################
## Anet,growth regression line prep
##########################################################################
anet.growth <- lmer(anet.growth ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(anet.growth, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
anet.growth.regline <- data.frame(emmeans(anet.growth, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Anet.growth plot
##########################################################################
anet.growth.plot <- ggplot(data = df, 
                    aes(x = n.trt, y = anet.growth, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = anet.growth.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = anet.growth.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("A")["net,growth"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment", shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
anet.growth.plot

##########################################################################
## Vcmax regression line prep
##########################################################################
vcmax25 <- lmer(vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(vcmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
vcmax25.regline <- data.frame(emmeans(vcmax25, ~co2*inoc, "n.trt",
                                     at = list(n.trt = seq(0, 630, 5)),
                                     type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Vcmax plot
##########################################################################
vcmax25.plot <- ggplot(data = df, 
                       aes(x = n.trt, y = vcmax25, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = vcmax25.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = vcmax25.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
vcmax25.plot

##########################################################################
## Jmax regression line prep
##########################################################################
jmax25 <- lmer(jmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(vcmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
jmax25.regline <- data.frame(emmeans(jmax25, ~co2*inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Jmax plot
##########################################################################
jmax25.plot <- ggplot(data = df, 
                      aes(x = n.trt, y = jmax25, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = jmax25.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = jmax25.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 240, 60)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
jmax25.plot

##########################################################################
## Jmax25:Vcmax25 regression line prep
##########################################################################
df$jmax25.vcmax25[100] <- NA
jvmax25 <- lmer(jmax25.vcmax25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(jvmax25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
jvmax25.regline <- data.frame(emmeans(jvmax25, ~co2*inoc, "n.trt",
                                      at = list(n.trt = seq(0, 630, 5)),
                                      type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Jmax:Vcmax plot
##########################################################################
jvmax25.plot <- ggplot(data = df,
                       aes(x = n.trt, y = jmax25.vcmax25, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = jvmax25.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = jvmax25.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(1.4, 2.2), breaks = seq(1.4, 2.2, 0.2)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"])),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
jvmax25.plot

##########################################################################
## Rd25 regression line prep
##########################################################################
df$rd25[df$rd25 < 0] <- NA
df$rd25[c(29, 34, 56)] <- NA

rd25 <- lmer(rd25 ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(rd25, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
rd25.regline <- data.frame(emmeans(rd25, ~co2*inoc, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_no.inoc", "solid", "dashed"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Rd25 plot
##########################################################################
rd25.plot <- ggplot(data = df,
                    aes(x = n.trt, y = rd25, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = rd25.regline,
              aes(linetype = linetype, color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = rd25.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 1.5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("R")["d25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
rd25.plot

##########################################################################
## PNUE regression line prep
##########################################################################
df$pnue.growth[41] <- NA

pnue <- lmer(pnue.growth ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(pnue, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
pnue.regline <- data.frame(emmeans(pnue, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_no.inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## PNUE plot
##########################################################################
pnue.plot <- ggplot(data = df,
                    aes(x = n.trt, y = pnue.growth, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = pnue.regline,
              aes(linetype = linetype, color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = pnue.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 6)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("PNUE")*" ("*mu*"mol CO"["2"]*" g"^"-1"*" N s"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
pnue.plot

##########################################################################
## PNUE regression line prep for co2-by-inoculation interaction plot 
## in supplement
##########################################################################
test(emtrends(pnue, pairwise~co2, "n.trt"))

## Emmean fxns for regression lines + error ribbons
pnue.int.regline <- data.frame(emmeans(pnue, ~co2, "n.trt",
                                        at = list(n.trt = seq(0, 630, 5)),
                                        type = "response")) %>%
  mutate(linetype = ifelse(co2 == "amb", "dashed", "solid"))

##########################################################################
## PNUE co2-by-inoculation interaction plot in supplement
##########################################################################
pnue.int.plot <- ggplot(data = df, aes(x = n.trt, y = pnue.growth, fill = co2)) +
  geom_point(aes(shape = inoc), size = 3, alpha = 0.75) +
  geom_smooth(data = pnue.int.regline,
              aes(color = co2, y = emmean, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = pnue.int.regline,
              aes(fill = co2, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_fill_manual(values = co2.cols,
                    labels = c("Ambient", "Elevated")) +
  scale_color_manual(values = co2.cols,
                     labels = c("Ambient", "Elevated")) +
  scale_shape_manual(values = c(21, 24), labels = c("Uninoculated",
                                                    "Inoculated")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 6)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("PNUE")*" ("*mu*"mol CO"["2"]*" g"^"-1"*" N s"^"-1"*")")),
       fill = expression(bold("CO"["2"])),
       color = expression(bold("CO"["2"])),
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25)) +
  guides(linetype = "none",
         fill = guide_legend(override.aes = list(shape = 21)))
pnue.int.plot

##########################################################################
## chi regression line prep
##########################################################################
chi <- lmer(chi ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(chi, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
chi.regline <- data.frame(emmeans(chi, ~co2*inoc, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "amb_inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## chi plot
##########################################################################
chi.plot <- ggplot(data = df, 
                   aes(x = n.trt, y = chi, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = chi.regline,
              aes(linetype = linetype, color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = chi.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0.49, 0.75), breaks = seq(0.5, 0.75, 0.05)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(chi*" (unitless)")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
chi.plot

##########################################################################
## Total leaf area regression line prep
##########################################################################
tla <- lmer(tla ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(tla, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
tla.regline <- data.frame(emmeans(tla, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Total leaf area plot
##########################################################################
tla.plot <- ggplot(data = df, aes(x = n.trt, y = tla, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = tla.regline,
              aes(color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = tla.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, 300)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
tla.plot

##########################################################################
## Total biomass regression line prep
##########################################################################
tbio <- lmer(total.biomass ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(tbio, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
tbio.regline <- data.frame(emmeans(tbio, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Total biomass plot
##########################################################################
tbio.plot <- ggplot(data = df,
                    aes(x = n.trt, y = total.biomass, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = tbio.regline,
              aes(color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = tbio.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Total biomass (g)",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
tbio.plot

##########################################################################
## Ncost regression line prep
##########################################################################
df$ncost[c(100, 101)] <- NA
df$ncost[c(38, 103)] <- NA
df$ncost[32] <- NA

ncost <- lmer(ncost ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(ncost, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
ncost.regline <- data.frame(emmeans(ncost, ~co2*inoc, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(co2.inoc == "elv_inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Ncost plot
##########################################################################
ncost.plot <- ggplot(data = df, aes(x = n.trt, y = ncost, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = ncost.regline,
              aes(color = co2.inoc, y = emmean, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = ncost.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["cost"]*" (gC gN"^"-1"*")")),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
ncost.plot

##########################################################################
## %Ndfa regression line prep
##########################################################################
df$ndfa[c(38, 85, 101, 103)] <- NA
ndfa <- lmer(sqrt(ndfa) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(ndfa, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
ndfa.regline <- data.frame(emmeans(ndfa, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "no.inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## %Ndfa plot
##########################################################################
ndfa.plot <- ggplot(data = df, 
                    aes(x = n.trt,  y = ndfa, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = ndfa.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = ndfa.regline,
              aes(fill = co2.inoc, y = response,
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("%N"["dfa"])),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
ndfa.plot

##########################################################################
## Belowground C allocation regression line prep for plot in supplement
##########################################################################
cbg <- lmer(log(cbg) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(cbg, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
cbg.regline <- data.frame(emmeans(cbg, ~co2*inoc, "n.trt",
                                    at = list(n.trt = seq(0, 630, 5)),
                                    type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Belowground C allocation plot in supplement
##########################################################################
cbg.plot <- ggplot(data = df, aes(x = n.trt, y = cbg, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = cbg.regline, aes(color = co2.inoc, y = response), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = cbg.regline, aes(fill = co2.inoc, y = response, 
                                      ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("C")["bg"]*" (gC)")),
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
cbg.plot

##########################################################################
## Whole-plant nitrogen biomass regression line prep for plot in supplement
##########################################################################
nwp <- lmer(sqrt(wpn) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nwp, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nwp.regline <- data.frame(emmeans(nwp, ~co2*inoc, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Nwp plot in supplement
##########################################################################
nwp.plot <- ggplot(data = df, aes(x = n.trt, y = wpn, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = nwp.regline, aes(color = co2.inoc, y = response), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = nwp.regline, aes(fill = co2.inoc, y = response, 
                                      ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold(italic("N")["wp"]*" (gN)")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
nwp.plot

##########################################################################
## Root nodule biomass regression line prep for plot in supplement
##########################################################################
df$nodule.biomass[80] <- NA

nod <- lmer(sqrt(nodule.biomass) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nod, ~co2, "n.trt"))
test(emtrends(nod, ~inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nod.regline <- data.frame(emmeans(nod, ~co2*inoc, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "no.inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Root nodule biomass plot in supplement
##########################################################################
nod.plot <- ggplot(data = df, 
                   aes(x = n.trt, y = nodule.biomass, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = nod.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = nod.regline,
              aes(fill = co2.inoc, y = response,
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(-0.01, 0.6), breaks = seq(0, 0.6, 0.15)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Nodule biomass (g)",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
nod.plot

##########################################################################
## Root nodule:root biomass regression line prep for plot in supplement
##########################################################################
nodroot <- lmer(sqrt(nod.root.ratio) ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(nodroot, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
nodroot.regline <- data.frame(emmeans(nodroot, ~co2*inoc, "n.trt",
                                   at = list(n.trt = seq(0, 630, 5)),
                                   type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         linetype = ifelse(inoc == "no.inoc", "dashed", "solid"),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Root nodule:root biomass plot in supplement
##########################################################################
nodroot.plot <- ggplot(data = df, 
                       aes(x = n.trt, y = nod.root.ratio, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = nodroot.regline,
              aes(color = co2.inoc, y = response, linetype = linetype), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = nodroot.regline,
              aes(fill = co2.inoc, y = response,
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(-0.0001, 0.4), breaks = seq(0, 0.4, 0.1)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Nodule: root biomass",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
nodroot.plot

##########################################################################
## BVR regression line prep for plot in supplement
##########################################################################
bvr <- lmer(bvr ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
shapiro.test(residuals(bvr))
test(emtrends(bvr, ~co2*inoc, "n.trt"))

## Emmean fxns for regression lines + error ribbons
bvr.regline <- data.frame(emmeans(bvr, ~co2*inoc, "n.trt",
                                  at = list(n.trt = seq(0, 630, 5)),
                                  type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## BVR plot in supplement
##########################################################################
bvr.plot <- ggplot(data = df, 
                    aes(x = n.trt,  y = bvr, fill = co2.inoc)) +
  geom_hline(yintercept = 1, col = "black", linetype = "dotted") +
  geom_hline(yintercept = 2, col = "black", linetype = "dashed") +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = bvr.regline,
              aes(color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = bvr.regline,
              aes(fill = co2.inoc, y = emmean,
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  labs(x = "Soil N fertilization (ppm)",
       y = expression(bold("BVR (g L"^"-1"*")")),
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 18) +
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
bvr.plot

##########################################################################
## Root biomass
##########################################################################
root.bio <- lmer(root.biomass ~ co2 * inoc * n.trt + (1|rack:co2), data = df)
test(emtrends(root.bio, ~co2*inoc, "n.trt"))


## Emmean fxns for regression lines + error ribbons
root.bio.regline <- data.frame(emmeans(root.bio, ~co2*inoc, "n.trt",
                                       at = list(n.trt = seq(0, 630, 5)),
                                       type = "response")) %>%
  mutate(co2.inoc = str_c(co2, "_", inoc),
         co2.inoc = factor(co2.inoc, 
                           levels = c("elv_inoc","elv_no.inoc", 
                                      "amb_inoc", "amb_no.inoc")))

##########################################################################
## Root biomass plot
##########################################################################
root.bio.plot <- ggplot(data = df, aes(x = n.trt, y = root.biomass, fill = co2.inoc)) +
  geom_jitter(aes(shape = inoc), size = 3, alpha = 0.75, width = 5) +
  geom_smooth(data = root.bio.regline,
              aes(color = co2.inoc, y = emmean), 
              linewidth = 1.5, se = FALSE) +
  geom_ribbon(data = root.bio.regline,
              aes(fill = co2.inoc, y = emmean, 
                  ymin = lower.CL, ymax = upper.CL), 
              linewidth = 1.5, alpha = 0.25) +
  scale_color_manual(values = full.cols,
                     labels = c(expression("Elevated CO"["2"]*", inoculated"),
                                expression("Elevated CO"["2"]*", uninoculated"),
                                expression("Ambient CO"["2"]*", inoculated"),
                                expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_fill_manual(values = full.cols,
                    labels = c(expression("Elevated CO"["2"]*", inoculated"),
                               expression("Elevated CO"["2"]*", uninoculated"),
                               expression("Ambient CO"["2"]*", inoculated"),
                               expression("Ambient CO"["2"]*", uninoculated"))) +
  scale_shape_manual(values = c(21, 24), 
                     labels = c("Uninoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Soil N fertilization (ppm)",
       y = "Root biomass (g)",
       fill = "Treatment", color = "Treatment",
       shape = "Inoculation") +
  theme_bw(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(linewidth = 1.25),
        legend.text.align = 0) +
  guides(linetype = "none", shape = "none",
         fill = guide_legend(override.aes = list(shape = c(24, 21, 24, 21))))
root.bio.plot

##########################################################################
## Figure 1: leaf N plots
##########################################################################
# png("[insert path here]",
#     height = 8, width = 12, units = "in", res = 600)
ggarrange(narea.plot, nmass.plot, marea.plot, chlarea.plot, 
          ncol = 2, nrow = 2, align = "hv", legend = "right",
          common.legend = TRUE, labels = c("(a)", "(b)", "(c)", "(d)"),
          font.label = list(size = 18))
# dev.off()

##########################################################################
## Figure 2: leaf physiology plots
##########################################################################
# png("[insert path here]",
#     height = 12, width = 12, units = "in", res = 600)
ggarrange(anet.plot, anet.growth.plot, vcmax25.plot, jmax25.plot,
          jvmax25.plot, rd25.plot,
          ncol = 2, nrow = 3, align = "hv", legend = "right",
          common.legend = TRUE, font.label = list(size = 18), 
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))
# dev.off()

##########################################################################
## Figure 3: nitrogen and water use efficiency plots
##########################################################################
# png("insert path here",
#     height = 4, width = 12, units = "in", res = 600)
ggarrange(pnue.plot, chi.plot,
          ncol = 2, nrow = 1, align = "hv", legend = "right",
          common.legend = TRUE, font.label = list(size = 18), 
          labels = c("(a)", "(b)"))
dev.off()

##########################################################################
## Figure 4: whole plant plots
##########################################################################
# png("[insert path here]",
#     height = 8, width = 12, units = "in", res = 600)
ggarrange(tla.plot, tbio.plot, ncost.plot, ndfa.plot,
          ncol = 2, nrow = 2, align = "hv", legend = "right",
          labels = c("(a)", "(b)", "(c)", "(d)"), common.legend = TRUE,
          font.label = list(size = 18))
# dev.off()

##########################################################################
## Figure S1: figure showing NxCO2 interaction for leaf N content
##########################################################################
# png("[insert path here]",
#     height = 7, width = 10, units = "in", res = 600)
ggarrange(narea.int.plot, nmass.int.plot, marea.int.plot, 
          ncol = 2, nrow = 2, align = "hv", legend = "right",
          common.legend = TRUE, labels = c("(a)", "(b)", "(c)"),
          font.label = list(size = 18))
# dev.off()

##########################################################################
## Figure S2: figure showing NxCO2 interaction for PNUE,growth
##########################################################################
# png("[insert path here]",
#     height = 4.5, width = 8, units = "in", res = 600)
pnue.int.plot
# dev.off()

##########################################################################
## Figure S3: figure showing components of Ncost (belowground C and whole
## plant N)
##########################################################################
# png("[insert path here]",
#    height = 4, width = 12, units = "in", res = 600)
ggarrange(cbg.plot, nwp.plot,
          align = "hv", common.legend = TRUE,
          nrow = 1, ncol = 2,
          legend = "right", labels = c("(a)", "(b)"), 
          font.label = list(size = 18))
# dev.off()

##########################################################################
## Figure S4: nitrogen fixation plots
##########################################################################
# png("[insert path here]",
#     height = 4, width = 12, units = "in", res = 600)
ggarrange(nod.plot, nodroot.plot,
          align = "hv", common.legend = TRUE,
          nrow = 1, ncol = 2,
          legend = "right", labels = c("(a)", "(b)"), 
          font.label = list(size = 18))
# dev.off()

##########################################################################
## Figure S5: BVR + root biomass
##########################################################################
png("../../2022_NxCO2xI/working_drafts/figs/NxCO2xI_figS5_bvr.png",
    height = 4.5, width = 12, units = "in", res = 600)
ggarrange(bvr.plot, root.bio.plot,
          align = "hv", common.legend = TRUE,
          nrow = 1, ncol = 2,
          legend = "right", labels = c("(a)", "(b)"), 
          font.label = list(size = 18))
dev.off()

