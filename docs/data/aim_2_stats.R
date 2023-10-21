library(ggplot2)
library(ggbreak)
library(ggprism)
library(ARTool)

pyr_labs <- c("Pyruvate", "Pyruvate in PC (0.04 mM)")
names(pyr_labs) <- c("Pyruvate", "Pyruvate in PC")

pyr_labs_full <- c("PC", "Pyruvate", "Pyruvate in PC (0.04 mM)")
names(pyr_labs_full) <- c("PC", "Pyruvate", "Pyruvate in PC")

source("data/Aim_2_Data_Gastroc_3Param.R")
source("data/Aim_2_Data_Soleus_3Param.R")


aim_2_combined_3Param <- rbind(Gast_palm_kinetic_values, Gast_pyr_kinetic_values, Gast_palm_pyr_kinetic_values, Sol_pyr_kinetic_values, Sol_palm_kinetic_values, Sol_palm_pyr_kinetic_values)

aim_2_combined_3Param$Smoke <- factor(aim_2_combined_3Param$Smoke, levels = unique(aim_2_combined_3Param$Smoke))
aim_2_combined_3Param$Tissue <- factor(aim_2_combined_3Param$Tissue, levels = unique(aim_2_combined_3Param$Tissue))
aim_2_combined_3Param$Substrate <- factor(aim_2_combined_3Param$Substrate, levels = unique(aim_2_combined_3Param$Substrate))
aim_2_combined_3Param$Subject <- factor(aim_2_combined_3Param$Subject, levels = unique(aim_2_combined_3Param$Subject))

# for(x in unique(aim_2_combined_3Param$Smoke)){
#   for(y in unique(aim_2_combined_3Param$Tissue)){
#     for(z in unique(aim_2_combined_3Param$Substrate)){
#       print(
#         identify_outliers(subset(subset(subset(aim_2_combined_3Param, Smoke %in% paste(x)), Tissue %in% paste(y)), Substrate %in% paste(z)), Km)
#         
#       )
#       
#       
#     }
#   }
# }
identify_outliers(subset(subset(subset(aim_2_combined_3Param, Smoke %in% "Control"), Tissue %in% "Soleus"), Substrate %in% "PC"), Km)

aim_2_combined_3Param$Km[aim_2_combined_3Param$Km == Gast_palm_pyr_kinetic_values[8,2]] <- NA
aim_2_combined_3Param$Vmax[aim_2_combined_3Param$Vmax == Gast_palm_pyr_kinetic_values[8,1]] <- NA
#aim_2_combined_3Param$Vmax[aim_2_combined_3Param$Vmax == Sol_palm_pyr_kinetic_values[14,1]] <- NA


aim_2_combined_3Param$Km[aim_2_combined_3Param$Km == Gast_palm_kinetic_values[5,2]] <- NA
#aim_2_combined_3Param$Vmax[aim_2_combined_3Param$Vmax == Gast_palm_kinetic_values[8,1]] <- NA

aim_2_combined_3Param$Km[aim_2_combined_3Param$Km == Sol_pyr_kinetic_values[10,2]] <- NA
#aim_2_combined_3Param$Km[aim_2_combined_3Param$Km == Sol_palm_pyr_kinetic_values[17,2]] <- NA

aim_2_combined_3Param$Km[aim_2_combined_3Param$Km == Sol_palm_kinetic_values[4, 2]] <- NA
aim_2_combined_3Param$Km[aim_2_combined_3Param$Km == Sol_palm_kinetic_values[5, 2]] <- NA
aim_2_combined_3Param$Km[aim_2_combined_3Param$Km == Sol_palm_kinetic_values[9, 2]] <- NA
aim_2_combined_3Param$Km[aim_2_combined_3Param$Km == Sol_palm_kinetic_values[17, 2]] <- NA

#aim_2_combined_3Param$Vmax[aim_2_combined_3Param$Vmax == Sol_palm_kinetic_values[9,1]] <- NA

## Normality & homogeneity ----
aim_2_combined_3Param |> 
  group_by(Tissue, Substrate) |> 
  levene_test(Vmax ~ Smoke)

as.data.table(aim_2_combined_3Param)[,.(Statistic = shapiro.test(Vmax)$statistic, 
                                 P.value = shapiro.test(Vmax)$p.value),
                              by = .(Tissue, Substrate, Smoke)]

aim_2_combined_3Param |> 
  group_by(Tissue, Substrate) |> 
  levene_test(Km ~ Smoke)

as.data.table(aim_2_combined_3Param)[,.(Statistic = shapiro.test(Km)$statistic, 
                                 P.value = shapiro.test(Km)$p.value),
                              by = .(Tissue, Substrate, Smoke)]

aim_2_combined_3Param |> 
  group_by(Tissue, Substrate) |> 
  levene_test(C ~ Smoke)

as.data.table(aim_2_combined_3Param)[,.(Statistic = shapiro.test(C)$statistic, 
                                        P.value = shapiro.test(C)$p.value),
                                     by = .(Tissue, Substrate, Smoke)]

# 
# effs_est <- rbind(Gast_pyr_kinetic_values, Gast_palm_pyr_kinetic_values, Sol_pyr_kinetic_values, Sol_palm_pyr_kinetic_values)
# 
# eff_results <- effectsize::cohens_f(lm(effs_est$Vmax ~ effs_est$Smoke * effs_est$Tissue * effs_est$Substrate))
# 
# eff_results1 <- effectsize::cohens_f(lm(effs_est1$Vmax ~ effs_est1$Smoke * effs_est1$Tissue * effs_est1$Substrate))
# 
# pwr::pwr.anova.test(k = 8, f = 0.41, sig.level = 0.05, power = 0.95)

# anova_km_3Param <- aov(lm(Km ~ Smoke * Tissue * Substrate, data = aim_2_combined_3Param))
# 
# anova_km_3Param2 <- aov(lm(Km ~ Smoke * Tissue * Substrate, data = subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "Pyruvate in PC"))))
# 
# anova_km_3ParamPC <- aov(lm(Km ~ Smoke * Tissue, data = subset(aim_2_combined_3Param, Substrate %in% c("PC"))))
# 
# anova_vmax_3Param <- aov(lm(Vmax ~ Smoke * Tissue * Substrate, data = aim_2_combined_3Param))
# 
# anova_vmax_3Param2 <- aov(lm(Vmax ~ Smoke * Tissue * Substrate, data = subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "Pyruvate in PC"))))
# 
# anova_vmax_3ParamPC <- aov(lm(Vmax ~ Smoke * Tissue, data = subset(aim_2_combined_3Param, Substrate %in% c("PC"))))
# 
# anova_summary_km_3Param <- summary(anova_km_3Param2)
# 
# anova_summary_vmax_3Param <- summary(anova_vmax_3Param2)
# 
# anova_PC_summary_km_3Param <- summary(anova_km_3ParamPC)
# 
# anova_PC_summary_vmax_3Param <- summary(anova_vmax_3ParamPC)
# 
# ART_PC_Vmax <- anova(ARTool::art(Vmax ~ Smoke * Tissue, data = na.omit(subset(aim_2_combined_3Param, Substrate %in% c("PC")))))
# 
# ART_Pyr_Km <- anova(ARTool::art(Km ~ Smoke * Tissue * Substrate, data = na.omit(subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "Pyruvate in PC")))))

#ANOVAs ----
## 3-way ANOVAs ----
### Km ----
ART_all_Km <- anova(art(Km ~ Smoke * Tissue * Substrate, data = na.omit(aim_2_combined_3Param)))

### Vmax ----
ART_all_vmax <- anova(art(Vmax ~ Smoke * Tissue * Substrate, data = na.omit(aim_2_combined_3Param)))

### C ----
ART_all_C <- anova(art(C ~ Smoke * Tissue * Substrate, data = na.omit(subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "PC")))))

## 2-way ANOVAS ----
### Km ----
all_km2way <- data.frame()
for(y in unique(aim_2_combined_3Param$Substrate)){
  km_2way <- anova(art(Km ~ Smoke * Tissue, data = na.omit(subset(aim_2_combined_3Param, Substrate %in% paste(y)))))
  km_2way$Substrate <- paste(y)
  
  all_km2way <- rbind(all_km2way, km_2way)
}

### Vmax ----
all_vmax2way <- data.frame()
for(y in c("Pyruvate", "Pyruvate in PC")){
  vmax_2way <- anova(lm(Vmax ~ Smoke * Tissue, data = na.omit(subset(aim_2_combined_3Param, Substrate %in% paste(y)))))
  vmax_2way$Substrate <- paste(y)
  
  all_vmax2way <- rbind(all_vmax2way, vmax_2way)
}

PCvmax_2way <- anova(art(Vmax ~ Smoke * Tissue, data = na.omit(subset(aim_2_combined_3Param, Substrate %in% "PC"))))
PCvmax_2way$Substrate <- "PC"

### C ----
all_c2way <- data.frame()
for(y in unique(aim_2_combined_3Param$Substrate)){
  c_2way <- anova(art(C ~ Smoke * Tissue, data = na.omit(subset(aim_2_combined_3Param, Substrate %in% paste(y)))))
  c_2way$Substrate <- paste(y)
  
  all_c2way <- rbind(all_c2way, c_2way)
}

# Graphs ----

# all_vmax_3Param <- ggplot(aim_2_combined_3Param, aes(x = Tissue, y = Vmax, fill = Smoke)) +
#   facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs_full)) +
#   stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
#   theme_prism() +
#   scale_fill_prism() +
#   coord_cartesian(ylim = c(0, 80)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 80, 10)) +
#   theme(axis.title.x = element_blank())


all_km_3Param <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "Pyruvate in PC")), aes(x = Tissue, y = Km, fill = Smoke)) +
  facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs_full)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 0.7)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .7, .1)) +
  theme(axis.title.x = element_blank())



# pyruvate_graphs_3Param <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "Pyruvate in PC")), aes(x = Tissue, y = Km, fill = Smoke)) +
#   facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs)) +
#   stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
#   theme_prism() +
#   scale_fill_prism() +
#   coord_cartesian(ylim = c(0, 0.7)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, .7, .1)) +
#   theme(axis.title.x = element_blank())
# 
# 
# palmitate_graph_3Param <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("PC")), aes(x = Tissue, y = Km, fill = Smoke)) +
#   facet_wrap(~Substrate) +
#   stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
#   theme_prism() +
#   scale_fill_prism() +
#   coord_cartesian(ylim = c(0, 0.08)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, .8, .01)) +
#   theme(axis.title.x = element_blank())


all_vmax_3Param <- ggplot(aim_2_combined_3Param, aes(x = Tissue, y = Vmax, color = Smoke)) +
  facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs_full)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), size = 1, fill = NA) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 90)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 90, 10)) +
  theme(axis.title.x = element_blank())

pyruvate_graphs_3Param_vmax <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "Pyruvate in PC")), aes(x = Tissue, y = Vmax, color = Smoke)) +
  facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), size = 1, fill = NA) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 100, 10)) +
  theme(axis.title.x = element_blank())

palmitate_graph_3Param_vmax <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("PC")), aes(x = Tissue, y = Vmax, color = Smoke)) +
  facet_wrap(~Substrate) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), size = 1, fill = NA) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 50)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 50, 5)) +
  theme(axis.title.x = element_blank())

pyruvate_graphs_3Param <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "Pyruvate in PC")), aes(x = Tissue, y = Km, color = Smoke)) +
  facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), size = 1, fill = NA) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.7)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .7, .1)) +
  theme(axis.title.x = element_blank())


palmitate_graph_3Param <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("PC")), aes(x = Tissue, y = Km, color = Smoke)) +
  facet_wrap(~Substrate) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), size = 1, fill = NA) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.06)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .8, .01)) +
  theme(axis.title.x = element_blank())

pyruvate_graphs_3Param_C <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("Pyruvate", "Pyruvate in PC")), aes(x = Tissue, y = Km, color = Smoke)) +
  facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), size = 1, fill = NA) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.7)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .7, .1)) +
  theme(axis.title.x = element_blank())


graph_3Param_C <- ggplot(aim_2_combined_3Param, aes(x = Tissue, y = C, color = Smoke)) +
  facet_wrap(~Substrate) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), size = 1, fill = NA) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 30)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 30, 5)) +
  theme(axis.title.x = element_blank()) +
  ylab("Initial Rate (C)")


# Percents ----  
library(scales)
percent_combined <- rbind(Gast_percent, Sol_percent)

percent_combined$Smoke <- factor(percent_combined$Smoke, levels = unique(percent_combined$Smoke))
percent_combined$Tissue <- factor(percent_combined$Tissue, levels = unique(percent_combined$Tissue))
percent_combined$Subject <- factor(percent_combined$Subject, levels = unique(percent_combined$Subject))
# identify_outliers(subset(subset(percent_combined, Smoke %in% "Control"), Tissue %in% "Soleus"), Vmax)


percent_combined |> 
  group_by(Tissue) |> 
  levene_test(Vmax ~ Smoke)

as.data.table(percent_combined)[,.(Statistic = shapiro.test(Vmax)$statistic, 
                                        P.value = shapiro.test(Vmax)$p.value),
                                     by = .(Tissue, Smoke)]

percent_combined |> 
  group_by(Tissue) |> 
  levene_test(Km ~ Smoke)

as.data.table(percent_combined)[,.(Statistic = shapiro.test(Km)$statistic, 
                                        P.value = shapiro.test(Km)$p.value),
                                     by = .(Tissue, Smoke)]

# Deltas ----  
delta_combined <- rbind(Gast_delta, Sol_delta)

delta_combined$Smoke <- factor(delta_combined$Smoke, levels = unique(delta_combined$Smoke))
delta_combined$Tissue <- factor(delta_combined$Tissue, levels = unique(delta_combined$Tissue))
delta_combined$Subject <- factor(delta_combined$Subject, levels = unique(delta_combined$Subject))
#identify_outliers(subset(subset(delta_combined, Smoke %in% "Control"), Tissue %in% "Gastrocnemius"), Km)


delta_combined |> 
  group_by(Tissue) |> 
  levene_test(Vmax ~ Smoke)

as.data.table(delta_combined)[,.(Statistic = shapiro.test(Vmax)$statistic, 
                                   P.value = shapiro.test(Vmax)$p.value),
                                by = .(Tissue, Smoke)]

delta_combined |> 
  group_by(Tissue) |> 
  levene_test(Km ~ Smoke)

as.data.table(delta_combined)[,.(Statistic = shapiro.test(Km)$statistic, 
                                   P.value = shapiro.test(Km)$p.value),
                                by = .(Tissue, Smoke)]





contributions_combined <- rbind(Gast_contributions, Sol_contributions)

contributions_combined$Smoke <- factor(contributions_combined$Smoke, levels = unique(contributions_combined$Smoke))
contributions_combined$Tissue <- factor(contributions_combined$Tissue, levels = unique(contributions_combined$Tissue))
contributions_combined$Subject <- factor(contributions_combined$Subject, levels = unique(contributions_combined$Subject))

contributions_longer <- pivot_longer(contributions_combined, 1:2, names_to = "Substrate", values_to = "Contribution")

contributions_longer$Condition <- paste(contributions_longer$Tissue, contributions_longer$Smoke)

contributions_longer$Substrate <- factor(contributions_longer$Substrate, levels = c("Pyruvate", "PC"))
# 2 parameter models ----

# source("data/Aim_2_Data_Gastroc.R")
# source("data/Aim_2_Data_Soleus.R")
# aim_2_combined <- rbind(Gast_palm_kinetic_values, Gast_pyr_kinetic_values, Gast_palm_pyr_kinetic_values, Sol_pyr_kinetic_values, Sol_palm_kinetic_values, Sol_palm_pyr_kinetic_values)
# 
# ## Normality & homogeneity ----
# aim_2_combined |> 
#   group_by(Tissue, Substrate) |> 
#   levene_test(Vmax ~ Smoke)
# 
# as.data.table(aim_2_combined)[,.(Statistic = shapiro.test(Vmax)$statistic, 
#                                  P.value = shapiro.test(Vmax)$p.value),
#                               by = .(Tissue, Substrate, Smoke)]
# 
# aim_2_combined |> 
#   group_by(Tissue, Substrate) |> 
#   levene_test(Km ~ Smoke)
# 
# as.data.table(aim_2_combined)[,.(Statistic = shapiro.test(Km)$statistic, 
#                                  P.value = shapiro.test(Km)$p.value),
#                               by = .(Tissue, Substrate, Smoke)]
# 
# # 
# # effs_est <- rbind(Gast_pyr_kinetic_values, Gast_palm_pyr_kinetic_values, Sol_pyr_kinetic_values, Sol_palm_pyr_kinetic_values)
# # 
# # eff_results <- effectsize::cohens_f(lm(effs_est$Vmax ~ effs_est$Smoke * effs_est$Tissue * effs_est$Substrate))
# # 
# # eff_results1 <- effectsize::cohens_f(lm(effs_est1$Vmax ~ effs_est1$Smoke * effs_est1$Tissue * effs_est1$Substrate))
# # 
# # pwr::pwr.anova.test(k = 8, f = 0.41, sig.level = 0.05, power = 0.95)
# 
# anova_km <- aov(lm(Km ~ Smoke * Tissue * Substrate, data = aim_2_combined))
# 
# anova_vmax <- aov(lm(Vmax ~ Smoke * Tissue * Substrate, data = aim_2_combined))
# 
# anova_summary_km <- summary(anova_km)
# 
# anova_summary_vmax <- summary(anova_vmax)
# 
# 
# 
# all_vmax <- ggplot(aim_2_combined, aes(x = Tissue, y = Vmax, fill = Smoke)) +
#   facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs_full)) +
#   stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
#   theme_prism() +
#   scale_fill_prism() +
#   coord_cartesian(ylim = c(0, 80)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 80, 10)) +
#   theme(axis.title.x = element_blank())
# 
# all_km <- ggplot(aim_2_combined, aes(x = Tissue, y = Km, fill = Smoke)) +
#   facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs_full)) +
#   stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
#   theme_prism() +
#   scale_fill_prism() +
#   scale_y_cut(c( 0.0005,0.003), scales = "free", space = 0.25) +
#   expand_limits(y = c(0, 0.3001)) +
#   theme(axis.title.x = element_blank())
# 
# pyruvate_graphs <- ggplot(subset(aim_2_combined, Substrate %in% c("Pyruvate", "Pyruvate in PC")), aes(x = Tissue, y = Km, fill = Smoke)) +
#   facet_wrap(~Substrate, labeller = labeller(Substrate = pyr_labs)) +
#   stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
#   theme_prism() +
#   scale_fill_prism() +
#   coord_cartesian(ylim = c(0, 0.3)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, .3, .05)) +
#   theme(axis.title.x = element_blank())
# 
# 
# palmitate_graph <- ggplot(subset(aim_2_combined, Substrate %in% c("PC")), aes(x = Tissue, y = Km, fill = Smoke)) +
#   facet_wrap(~Substrate) +
#   stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
#   theme_prism() +
#   scale_fill_prism() +
#   coord_cartesian(ylim = c(0, 0.0016)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, .002, .0002)) +
#   theme(axis.title.x = element_blank())