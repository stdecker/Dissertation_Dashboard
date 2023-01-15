# Final Graphs & Stats ----
library(patchwork)
library(dunn.test)
library(ggbeeswarm)

source("data/aim_2_stats.R")

Sidak <- function(vecP)
  #
  # This function corrects a vector of probabilities for multiple testing
  # using the Bonferroni (1935) and Sidak (1967) corrections.
  #
  # References: Bonferroni (1935), Sidak (1967), Wright (1992).
  #
  # Bonferroni, C. E. 1935. Il calcolo delle assicurazioni su gruppi di teste. 
  # Pp. 13-60 in: Studi in onore del Professore Salvatore Ortu Carboni. Roma.
  #
  # Sidak, Z. 1967. Rectangular confidence regions for the means of multivariate 
  # normal distributions. Journal of the American Statistical Association 62:626-633.
#
# Wright, S. P. 1992. Adjusted P-values for simultaneous inference. 
# Biometrics 48: 1005-1013. 
#
#                  Pierre Legendre, May 2007
{
  vecP = sort(vecP, decreasing = F)  
  k = length(vecP)
  
  vecPB = 0
  vecPS = 0
  
  for(i in 1:k) {
    bonf = vecP[i]*k
    if(bonf > 1) bonf=1
    vecPB = c(vecPB, bonf)
    vecPS = c(vecPS, (1-(1-vecP[i])^(k-(i-1))))
  }
  #
  return(list(OriginalP=vecP, BonfP=vecPB[-1], SidakP=vecPS[-1]))
}

aim_2_combined_3Param$Condition <- paste(aim_2_combined_3Param$Tissue, aim_2_combined_3Param$Smoke)

#Vmax ----

all_vmax2way <- data.frame()
for(y in unique(aim_2_combined_3Param$Substrate)){
  vmax_2way <- anova(art(Vmax ~ Smoke * Tissue, data = na.omit(subset(aim_2_combined_3Param[-c(2)], Substrate %in% paste(y)))))
  vmax_2way$Substrate <- paste(y)
  
  print(effectsize::eta_squared(vmax_2way))
  all_vmax2way <- rbind(all_vmax2way, vmax_2way)
}

Pyruvate_dunn_vmax <- na.omit(subset(aim_2_combined_3Param[-c(2)], Substrate %in% "Pyruvate"))
Pyruvate_dunn_vmax$Condition <- paste(Pyruvate_dunn_vmax$Tissue, Pyruvate_dunn_vmax$Smoke)
dunn.test(Pyruvate_dunn_vmax$Vmax, Pyruvate_dunn_vmax$Condition, method = "hs", list = T)

wide_pyr_vmax <- pivot_wider(Pyruvate_dunn_vmax, names_from = Condition, id_cols = Subject, values_from = c(Vmax))
effectsize::cohens_d(wide_pyr_vmax$`Soleus Control`, wide_pyr_vmax$`Soleus Smoke`, pooled_sd = F)

pyruvate_vmax <- ggplot(subset(aim_2_combined_3Param, Substrate %in% "Pyruvate"), aes(x = Tissue, y = Vmax, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = 'white') +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 90)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 90, 10)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(V[max]~(pmol[O['2']]/sec/mg[wt])))) +
  # ggtitle("Pyruvate") +
  # geom_point(data = subset(aim_2_combined_3Param, Substrate %in% "Pyruvate"), aes(x = Tissue, y = Vmax, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = subset(aim_2_combined_3Param, Substrate %in% "Pyruvate"), aes(x = Tissue, y = Vmax, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.45, y = 88, label = 'Smoke Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 84, label = 'Tissue Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 80, label = 'Interaction Effect: p = 0.012', hjust = 0, fontface = 2) +  
  # annotate('text', x = 1.15, y = 35, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 1.3, y = 37, label = '†', fontface = 2, hjust = 0.5, size = 3.5) +  
  # annotate('text', x = 1.23, y = 36.5, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  geom_bracket(xmin = c(.755), xmax = c(1.245), y.position = 70, inherit.aes = F, label = "p = 0.026", size = 1, fontface = 2)# +
  #  # geom_vline(xintercept = 1.5, size = 1, alpha = 0.5, linetype = "longdash") #+   # +
  # geom_bracket(xmin = c(1.245), xmax = c(1.755), y.position = 90, inherit.aes = F, label = "p < 0.001", size = 1, fontface = 2) +
  # geom_bracket(xmin = c(1.245), xmax = c(2.255), y.position = 97, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2)

pyrPC_dunn_vmax <- na.omit(subset(aim_2_combined_3Param[-c(2)], Substrate %in% "Pyruvate in PC"))
pyrPC_dunn_vmax$Condition <- paste(pyrPC_dunn_vmax$Tissue, pyrPC_dunn_vmax$Smoke)
pyrPC_dunn_results_vmax <- dunn.test(pyrPC_dunn_vmax$Vmax, pyrPC_dunn_vmax$Condition, method = "none", list = T)
Sidak(c(pyrPC_dunn_results_vmax$P[2], pyrPC_dunn_results_vmax$P[5]))
Sidak(c(pyrPC_dunn_results_vmax$P[1], pyrPC_dunn_results_vmax$P[6]))

wide_pyrPC_vmax <- pivot_wider(pyrPC_dunn_vmax, names_from = Condition, id_cols = Subject, values_from = c(Vmax))
effectsize::cohens_d(wide_pyrPC_vmax$`Soleus Smoke`, wide_pyrPC_vmax$`Soleus Control`, pooled_sd = F)

pyruvatePC_vmax <- ggplot(subset(aim_2_combined_3Param, Substrate %in% "Pyruvate in PC"), aes(x = Tissue, y = Vmax, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 90)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 90, 10)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())  +
  labs(y=expression(bold(V[max]~(pmol[O['2']]/sec/mg[wt])))) +
  # ggtitle("Pyruvate + PC (0.04 mM)") +
  # geom_point(data = subset(aim_2_combined_3Param, Substrate %in% "Pyruvate in PC"), aes(x = Tissue, y = Vmax, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = subset(aim_2_combined_3Param, Substrate %in% "Pyruvate in PC"), aes(x = Tissue, y = Vmax, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.45, y = 88, label = 'Smoke Main Effect: p = 0.035', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 84, label = 'Tissue Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 80, label = 'Interaction Effect: p = 0.276', hjust = 0, fontface = 2)# +
   # geom_vline(xintercept = 1.5, size = 1, alpha = 0.5, linetype = "longdash") #+   #+  
  # annotate('text', x = 1.78, y = 90, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 2.23, y = 87, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  # geom_bracket(xmin = c(.775), xmax = c(1.225), y.position = 58, inherit.aes = F, label = "p = 0.104", size = 1, fontface = 2) +
  # geom_bracket(xmin = c(.755), xmax = c(1.755), y.position = 92, inherit.aes = F, label = "p = 0.005", size = 1, fontface = 2) +
  # geom_bracket(xmin = c(1.245), xmax = c(2.255), y.position = 101, inherit.aes = F, label = "p < 0.001", size = 1, fontface = 2)

PC_dunn_vmax <- na.omit(subset(aim_2_combined_3Param[-c(2)], Substrate %in% "PC"))
PC_dunn_vmax$Condition <- paste(PC_dunn_vmax$Tissue, PC_dunn_vmax$Smoke)
PC_dunn_results_vmax <- dunn.test(PC_dunn_vmax$Vmax, PC_dunn_vmax$Condition, method = "none", list = T)
Sidak(c(PC_dunn_results_vmax$P[2], PC_dunn_results_vmax$P[5]))

wide_PC_vmax <- pivot_wider(PC_dunn_vmax, names_from = Condition, id_cols = Subject, values_from = c(Vmax))
effectsize::cohens_d(wide_PC_vmax$`Gastrocnemius Smoke`, wide_PC_vmax$`Soleus Smoke`, pooled_sd = F)

palmitate_vmax <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("PC")), aes(x = Tissue, y = Vmax, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 45)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 45, 5)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank())  +
  labs(y=expression(bold(V[max]~(pmol[O['2']]/sec/mg[wt])))) +
  # ggtitle("PC") +
  # geom_point(data = subset(aim_2_combined_3Param, Substrate %in% "PC"), aes(x = Tissue, y = Vmax, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = subset(aim_2_combined_3Param, Substrate %in% "PC"), aes(x = Tissue, y = Vmax, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.45, y = 44, label = 'Smoke Main Effect: p = 0.240', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 42, label = 'Tissue Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 40, label = 'Interaction Effect: p = 0.952', hjust = 0, fontface = 2)# +
   # geom_vline(xintercept = 1.5, size = 1, alpha = 0.5, linetype = "longdash") #+   #+  
  # annotate('text', x = 1.78, y = 41.5, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 2.23, y = 45, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  # geom_bracket(xmin = c(1.245), xmax = c(2.245), y.position = 49, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2) +
  # geom_bracket(xmin = c(.755), xmax = c(1.755), y.position = 43, inherit.aes = F, label = "p = 0.010", size = 1, fontface = 2)


#Km ----
all_km2way <- data.frame()
for(y in unique(aim_2_combined_3Param$Substrate)){
  km_2way <- anova(art(Km ~ Smoke * Tissue, data = na.omit(subset(aim_2_combined_3Param[-c(1)], Substrate %in% paste(y)))))
  km_2way$Substrate <- paste(y)
  
  print(effectsize::eta_squared(km_2way))
  all_km2way <- rbind(all_km2way, km_2way)
}

Pyruvate_dunn_km <- na.omit(subset(aim_2_combined_3Param[-c(1)], Substrate %in% "Pyruvate"))
Pyruvate_dunn_km$Condition <- paste(Pyruvate_dunn_km$Tissue, Pyruvate_dunn_km$Smoke)
dunn.test(Pyruvate_dunn_km$Km, Pyruvate_dunn_km$Condition, method = "hs", list = T)

wide_pyr_km <- pivot_wider(Pyruvate_dunn_km, names_from = Condition, id_cols = Subject, values_from = c(Km))
effectsize::cohens_d(wide_pyr_km$`Soleus Smoke`, wide_pyr_km$`Soleus Control`, pooled_sd = F)

pyruvate_km <- ggplot(subset(aim_2_combined_3Param, Substrate %in% "Pyruvate"), aes(x = Tissue, y = Km, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, .9)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .9, .1)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(K[m]~'([Pyruvate]'~'mM)'))) +
  # ggtitle("Pyruvate") +
  # geom_point(data = subset(aim_2_combined_3Param, Substrate %in% "Pyruvate"), aes(x = Tissue, y = Km, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = subset(aim_2_combined_3Param, Substrate %in% "Pyruvate"), aes(x = Tissue, y = Km, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.45, y = 0.88, label = 'Smoke Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 0.84, label = 'Tissue Main Effect: p = 0.008', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 0.8, label = 'Interaction Effect: p = 0.008', hjust = 0, fontface = 2) +  
  # annotate('text', x = 1.2, y = .71, label = '*', fontface = 2, hjust = 0.5, size = 6) +  
  # annotate('text', x = 1.26, y = .725, label = '#', fontface = 2, hjust = 0.5, size = 4) +  
  # annotate('text', x = 2.2, y = .565, label = '*', fontface = 2, hjust = 0.5, size = 6) +  
  # annotate('text', x = 2.26, y = .58, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  geom_bracket(xmin = c(.755), xmax = c(1.23), y.position = 0.72, inherit.aes = F, label = "p < 0.001", size = 1, fontface = 2) +
  #geom_bracket(xmin = c(1.27), xmax = c(1.74), y.position = 0.72, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2) +
  geom_bracket(xmin = c(1.78), xmax = c(2.245), y.position = 0.6, inherit.aes = F, label = "p = 0.080", size = 1, fontface = 2)# +
   # geom_vline(xintercept = 1.5, size = 1, alpha = 0.5, linetype = "longdash") #+  # +
  #geom_bracket(xmin = c(.755), xmax = c(2.245), y.position = 0.81, inherit.aes = F, label = "p = 0.014", size = 1, fontface = 2)

pyrPC_dunn_km <- na.omit(subset(aim_2_combined_3Param[-c(1)], Substrate %in% "Pyruvate in PC"))
pyrPC_dunn_km$Condition <- paste(pyrPC_dunn_km$Tissue, pyrPC_dunn_km$Smoke)
pyrPC_dunn_results_km <- dunn.test(pyrPC_dunn_km$Km, pyrPC_dunn_km$Condition, method = "none", list = T)
Sidak(c(pyrPC_dunn_results_km$P[1], pyrPC_dunn_results_km$P[6]))
Sidak(c(pyrPC_dunn_results_km$P[2], pyrPC_dunn_results_km$P[5]))

wide_pyrPC_km <- pivot_wider(pyrPC_dunn_km, names_from = Condition, id_cols = Subject, values_from = c(Km))
effectsize::cohens_d(wide_pyrPC_km$`Gastrocnemius Smoke`, wide_pyrPC_km$`Soleus Smoke`, pooled_sd = F)

pyruvatePC_km <- ggplot(subset(aim_2_combined_3Param, Substrate %in% "Pyruvate in PC"), aes(x = Tissue, y = Km, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .5, .1)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(K[m]~'([Pyruvate]'~'mM)'))) +
  # ggtitle("Pyruvate + PC (0.04 mM)") +
  # geom_point(data = subset(aim_2_combined_3Param, Substrate %in% "Pyruvate in PC"), aes(x = Tissue, y = Km, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = subset(aim_2_combined_3Param, Substrate %in% "Pyruvate in PC"), aes(x = Tissue, y = Km, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.45, y = .49, label = 'Smoke Main Effect: p = 0.038', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = .465, label = 'Tissue Main Effect: p = 0.048', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = .44, label = 'Interaction Effect: p = 0.344', hjust = 0, fontface = 2)# +
   # geom_vline(xintercept = 1.5, size = 1, alpha = 0.5, linetype = "longdash") #+   #+   
  # annotate('text', x = 2.23, y = .51, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  # geom_bracket(xmin = c(.775), xmax = c(2.225), y.position = 0.81, inherit.aes = F, label = "p = 0.014", size = 1, fontface = 2) +
  # geom_bracket(xmin = c(.775), xmax = c(2.225), y.position = 0.81, inherit.aes = F, label = "p = 0.014", size = 1, fontface = 2)

PC_dunn_km <- na.omit(subset(aim_2_combined_3Param[-c(1)], Substrate %in% "PC"))
PC_dunn_km$Condition <- paste(PC_dunn_km$Tissue, PC_dunn_km$Smoke)
PC_dunn_results_km <- dunn.test(PC_dunn_km$Km, PC_dunn_km$Condition, method = "none", list = T)
Sidak(c(PC_dunn_results_km$P[2], PC_dunn_results_km$P[5]))

wide_PC_km <- pivot_wider(PC_dunn_km, names_from = Condition, id_cols = Subject, values_from = c(Km))
effectsize::cohens_d(wide_PC_km$`Soleus Smoke`, wide_PC_km$`Gastrocnemius Smoke`, pooled_sd = F)

palmitate_km <- ggplot(subset(aim_2_combined_3Param, Substrate %in% c("PC")), aes(x = Tissue, y = Km, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 0.05)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .05, .005)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(K[m]~'([PC]'~'mM)'))) +
  # ggtitle("PC") +
  # geom_point(data = subset(aim_2_combined_3Param, Substrate %in% "PC"), aes(x = Tissue, y = Km, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = subset(aim_2_combined_3Param, Substrate %in% "PC"), aes(x = Tissue, y = Km, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.45, y = .049, label = 'Smoke Main Effect: p = 0.220', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = .0465, label = 'Tissue Main Effect: p = 0.002', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = .044, label = 'Interaction Effect: p = 0.220', hjust = 0, fontface = 2)# +
   # geom_vline(xintercept = 1.5, size = 1, alpha = 0.5, linetype = "longdash") #+   #+  
  #annotate('text', x = 1.78, y = 0.017, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 2.23, y = 0.051, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  # geom_bracket(xmin = c(.75), xmax = c(1.75), y.position = 0.016, inherit.aes = F, label = "p = 0.014", size = 1, fontface = 2) +
  # geom_bracket(xmin = c(1.245), xmax = c(2.25), y.position = 0.049, inherit.aes = F, label = "p = 0.009", size = 1, fontface = 2)




#Percents ----

percent_vmax2way <- anova(art(Vmax ~ Smoke * Tissue, data = percent_combined))
percent_km2way <- anova(art(Km ~ Smoke * Tissue, data = percent_combined[-c(18,34),]))

effectsize::eta_squared(art(Vmax ~ Smoke * Tissue, data = percent_combined))
effectsize::eta_squared(art(Km ~ Smoke * Tissue, data = percent_combined[-c(18,34),]))

percent_combined$Condition <- paste(percent_combined$Tissue, percent_combined$Smoke)
percent_vmax_results<- dunn.test(percent_combined$Vmax, percent_combined$Condition, method = "none", list = T)
Sidak(c(percent_vmax_results$P[1], percent_vmax_results$P[6]))

wide_per_vmax <- pivot_wider(percent_combined, names_from = Condition, id_cols = Subject, values_from = c(Vmax))
effectsize::cohens_d(wide_per_vmax$`Gastrocnemius Control`, wide_per_vmax$`Gastrocnemius Smoke`, pooled_sd = F)

percent_km_results <- dunn.test(percent_combined$Km[-c(18,34)], percent_combined$Condition[-c(18,34)], method = "none", list = T)
Sidak(c(percent_km_results$P[1], percent_km_results$P[6]))

wide_per_km <- pivot_wider(percent_combined[-c(18,34),], names_from = Condition, id_cols = Subject, values_from = c(Km))
effectsize::cohens_d(wide_per_km$`Soleus Control`, wide_per_km$`Soleus Smoke`)

percent_combined$Condition <- paste(percent_combined$Tissue, percent_combined$Smoke)

pyr_percent_km <- ggplot(percent_combined[-c(18,34),], aes(x = Tissue, y = Km, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 3.5)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 3.5, .5), labels = label_percent()) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  # labs(y=expression(bold(paste(Change~'in'~K[m]('(Pyruvate + PC)/Pyruvate'))))) +
  labs(y=expression(bold(paste(over(Pyruvate + PC~K[m],Pyruvate~K[m]))))) +
  # geom_point(data = percent_combined[-c(34),], aes(x = Tissue, y = Km, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = percent_combined[-c(34),], aes(x = Tissue, y = Km, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.45, y = 3.4, label = 'Smoke Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 3.25, label = 'Tissue Main Effect: p = 0.102', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 3.1, label = 'Interaction Effect: p = 0.143', hjust = 0, fontface = 2) +  
  # annotate('text', x = 1.17, y = 3.35, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 1.29, y = 3.41, label = '†', fontface = 2, hjust = 0.5, size = 3.5) +  
  # annotate('text', x = 1.23, y = 3.4, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  geom_hline(yintercept = 1, linetype = "longdash", alpha = 0.5, size = 1) +  
  #annotate('text', x = 1.23, y = .9, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_bracket(xmin = c(.755), xmax = c(1.245), y.position = 2.72, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = .02) +
  geom_bracket(xmin = c(1.75), xmax = c(2.245), y.position = 1.8, inherit.aes = F, label = "p = 0.053", size = 1, fontface = 2, tip.length = .02)


pyr_percent_vmax <- ggplot(percent_combined, aes(x = Tissue, y = Vmax, fill = Smoke)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, color = "black") +
  theme_prism() +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(0, 2.5)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 2.5, .5), labels = label_percent()) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  # labs(y=expression(bold(Change~'in'~V[max]~('(Pyruvate + PC)/Pyruvate')))) +
  labs(y=expression(bold(paste(over(Pyruvate + PC~V[max],Pyruvate~V[max]))))) +
  new_scale_fill() +
  # geom_point(data = percent_combined, aes(x = Tissue, y = Vmax, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = percent_combined, aes(x = Tissue, y = Vmax, fill = Smoke), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("white", "white")) +
  annotate('text', x = 0.45, y = 2.42, label = 'Smoke Main Effect: p = 0.027', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 2.32, label = 'Tissue Main Effect: p = 0.843', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 2.22, label = 'Interaction Effect: p = 0.095', hjust = 0, fontface = 2) +  
  # annotate('text', x = 1.17, y = 1.17, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 1.29, y = 1.21, label = '†', fontface = 2, hjust = 0.5, size = 3.5) +  
  # annotate('text', x = 1.23, y = 1.2, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  # annotate('text', x = 1.775, y = 1.6, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  # annotate('text', x = 2.23, y = 1.2, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_hline(yintercept = 1, linetype = "longdash", alpha = 0.5, size = 1) +  
  #annotate('text', x = 1.23, y = 2, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_bracket(xmin = c(.755), xmax = c(1.245), y.position = 1.97, inherit.aes = F, label = "p = 0.003", size = 1, fontface = 2, tip.length = 0.015)



percent_combined$Km2 <- percent_combined$Km -1
percent_combined$Vmax2 <- percent_combined$Vmax -1

pyr_percent_km2 <- ggplot(percent_combined[-c(18,34),], aes(x = Tissue, y = Km2, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(-1, 2)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1, 2, .5), labels = label_percent()) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  # labs(y=expression(bold(Change~'in'~K[m]~('(Pyruvate + PC)/Pyruvate')))) +
  labs(y=expression(bold(paste(over(Pyruvate + PC~K[m],Pyruvate~K[m]),~'(',Percent~Change,')')))) +
  # geom_point(data = percent_combined[-c(34),], aes(x = Tissue, y = Km, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = percent_combined[-c(34),], aes(x = Tissue, y = Km2, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.55, y = 1.92, label = 'Smoke Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 1.81, label = 'Tissue Main Effect: p = 0.102', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 1.7, label = 'Interaction Effect: p = 0.143', hjust = 0, fontface = 2) +  
  # annotate('text', x = 1.17, y = 3.35, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 1.29, y = 3.41, label = '†', fontface = 2, hjust = 0.5, size = 3.5) +  
  # annotate('text', x = 1.23, y = 3.4, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  geom_hline(yintercept = 0, size = 1) +  
  #annotate('text', x = 1.23, y = .9, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_bracket(xmin = c(.755), xmax = c(1.245), y.position = 1.72, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = .02) +
  geom_bracket(xmin = c(1.75), xmax = c(2.245), y.position = .8, inherit.aes = F, label = "p = 0.053", size = 1, fontface = 2, tip.length = .02)# +
   # geom_vline(xintercept = 1.5, size = 1, alpha = 0.5, linetype = "longdash") #+  


pyr_percent_vmax2 <- ggplot(percent_combined, aes(x = Tissue, y = Vmax2, color = Condition)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(-1, 1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1, 1, .25), labels = label_percent()) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  # labs(y=expression(bold(Change~'in'~V[max]~('(Pyruvate + PC)/Pyruvate')))) +
  labs(y=expression(bold(paste(over(Pyruvate + PC~V[max],Pyruvate~V[max]),~'(',Percent~Change,')')))) +
  # geom_point(data = percent_combined, aes(x = Tissue, y = Vmax, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = percent_combined, aes(x = Tissue, y = Vmax2, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  annotate('text', x = 0.55, y = 0.95, label = 'Smoke Main Effect: p = 0.027', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 0.88, label = 'Tissue Main Effect: p = 0.843', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 0.81, label = 'Interaction Effect: p = 0.095', hjust = 0, fontface = 2) +  
  # annotate('text', x = 1.17, y = 1.17, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 1.29, y = 1.21, label = '†', fontface = 2, hjust = 0.5, size = 3.5) +  
  # annotate('text', x = 1.23, y = 1.2, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  # annotate('text', x = 1.775, y = 1.6, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  # annotate('text', x = 2.23, y = 1.2, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_hline(yintercept = 0, size = 1) +  
  #annotate('text', x = 1.23, y = 2, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_bracket(xmin = c(.755), xmax = c(1.245), y.position = .97, inherit.aes = F, label = "p = 0.003", size = 1, fontface = 2, tip.length = 0.015)# +
   # geom_vline(xintercept = 1.5, size = 1, alpha = 0.5, linetype = "longdash") #+  





percent_combined$Km_log <- log(percent_combined$Km, 2)
percent_combined$Vmax_log <- log(percent_combined$Vmax, 2)

pyr_percent_km_log <- ggplot(percent_combined[-c(18,34),], aes(x = Tissue, y = Km_log, fill = Smoke)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, color = "black") +
  theme_prism() +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(-3, 2.5), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-3, 2.5, .5)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(Change~'in'~K[m]~'('*log['2']*'((Pyruvate + PC)/Pyruvate))'))) +
  new_scale_fill() +
  # geom_point(data = percent_combined[-c(34),], aes(x = Tissue, y = Km, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = percent_combined[-c(34),], aes(x = Tissue, y = Km_log, fill = Smoke), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("white", "white")) +
  annotate('text', x = 0.45, y = 2.45, label = 'Smoke Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 2.25, label = 'Tissue Main Effect: p = 0.102', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 2.05, label = 'Interaction Effect: p = 0.143', hjust = 0, fontface = 2) +  
  # annotate('text', x = 1.17, y = 3.35, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 1.29, y = 3.41, label = '†', fontface = 2, hjust = 0.5, size = 3.5) +  
  # annotate('text', x = 1.23, y = 3.4, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  geom_hline(yintercept = 0, size = 1) +  
  #annotate('text', x = 1.23, y = .9, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_bracket(xmin = c(.755), xmax = c(1.245), y.position = 1.72, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = .02) +
  geom_bracket(xmin = c(1.75), xmax = c(2.245), y.position = 2.1, inherit.aes = F, label = "p = 0.053", size = 1, fontface = 2, tip.length = .02)


pyr_percent_vmax_log <- ggplot(percent_combined, aes(x = Tissue, y = Vmax_log, fill = Smoke)) +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, color = "black") +
  theme_prism() +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(-2, 1.5), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-2, 1.5, .5)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(Change~'in'~V[max]~'('*log['2']*'((Pyruvate + PC)/Pyruvate))'))) +
  new_scale_fill() +
  # geom_point(data = percent_combined, aes(x = Tissue, y = Vmax, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = percent_combined, aes(x = Tissue, y = Vmax_log, fill = Smoke), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("white", "white")) +
  annotate('text', x = 0.45, y = 1.48, label = 'Smoke Main Effect: p = 0.027', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 1.36, label = 'Tissue Main Effect: p = 0.843', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 1.24, label = 'Interaction Effect: p = 0.095', hjust = 0, fontface = 2) +  
  # annotate('text', x = 1.17, y = 1.17, label = '*', fontface = 2, hjust = 0.5, size = 6) + 
  # annotate('text', x = 1.29, y = 1.21, label = '†', fontface = 2, hjust = 0.5, size = 3.5) +  
  # annotate('text', x = 1.23, y = 1.2, label = '#', fontface = 2, hjust = 0.5, size = 4) +
  # annotate('text', x = 1.775, y = 1.6, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  # annotate('text', x = 2.23, y = 1.2, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_hline(yintercept = 0, size = 1) +  
  #annotate('text', x = 1.23, y = 2, label = '*', fontface = 2, hjust = 0.5, size = 6) +
  geom_bracket(xmin = c(.755), xmax = c(1.245), y.position = 1.05, inherit.aes = F, label = "p = 0.003", size = 1, fontface = 2, tip.length = 0.015)




#means ----
library(plotrix)

vmax_means <- aggregate(x = Vmax ~ Smoke + Tissue + Substrate, data = aim_2_combined_3Param, FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
km_means <- aggregate(x = Km ~ Smoke + Tissue + Substrate, data = aim_2_combined_3Param, FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))


per_vmax_means <- aggregate(x = Vmax2 ~ Smoke + Tissue, data = percent_combined, FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
per_km_means <- aggregate(x = Km2 ~ Smoke + Tissue, data = percent_combined[-c(18,34),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))



# Deltas ----

delta_vmax2way <- anova(art(Vmax ~ Smoke * Tissue, data = delta_combined))
delta_km2way <- anova(art(Km ~ Smoke * Tissue, data = delta_combined))

pyr_delta_km <- ggplot(delta_combined, aes(x = Tissue, y = Km, fill = Smoke)) +
  geom_hline(yintercept = 0, size = 1) +  
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, color = "black") +
  theme_prism() +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(-.80, .20)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-.80, .20, .10)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(Delta~K[m]))) +
  new_scale_fill() +
  # geom_point(data = delta_combined, aes(x = Tissue, y = Km, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = delta_combined, aes(x = Tissue, y = Km, fill = Smoke), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("white", "white"))


pyr_delta_vmax <- ggplot(delta_combined, aes(x = Tissue, y = Vmax, fill = Smoke)) +
  geom_hline(yintercept = 0, size = 1) +  
  #stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, color = "black") +
  theme_prism() +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(-50, 50)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-50, 50, 10)) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(Delta~V[max]))) +
  new_scale_fill() +
  # geom_point(data = delta_combined, aes(x = Tissue, y = Vmax, fill = Smoke), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = delta_combined, aes(x = Tissue, y = Vmax, fill = Smoke), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("white", "white"))




## Contributions ----
contributions_km2way <- anova(art(PC ~ Smoke * Tissue, data = contributions_combined))
effectsize::eta_squared(art(PC ~ Smoke * Tissue, data = contributions_combined))

contributions_plot <- ggplot(contributions_longer, aes(x = Smoke, y = Contribution, fill = Substrate)) +
  stat_summary(geom = "bar", size = 1, color = "black") +
  #stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.3, size = 1, na.rm = TRUE) +
  theme_prism() +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, .10), labels = label_percent()) +
  theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
  labs(y=expression(bold(Contrubution~to~Maximal~Respiration))) +
  facet_grid(~Tissue) +
  new_scale_fill() +
  # geom_point(data = subset(contributions_longer, Substrate %in% "PC"), aes(x = Smoke, y = Contribution, fill = Substrate), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(cex = 3, data = subset(contributions_longer, Substrate %in% "PC"), aes(x = Smoke, y = Contribution, fill = Substrate), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("white", "white"))

# ggplot(contributions_longer, aes(x = Smoke, y = Contribution, fill = Substrate)) +
#   geom_bar_pattern(stat = "summary", color = "black", size = 1, position = "identity", aes(pattern = Substrate)) +
#   theme_prism() +
#   scale_fill_manual(values = c("white", "grey25")) +
#   coord_cartesian(ylim = c(0, 1), clip = "off") +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, .10), labels = label_percent()) +
#   theme(axis.title.x = element_blank(),         axis.text.x = element_blank(),         axis.ticks.x = element_blank()) +
#   labs(y=expression(bold(Contrubution~to~Maximal~Respiration))) +
#   facet_grid(~Tissue) +
#   new_scale_fill() +
#   stat_summary(geom = "bar", size = 1, color = "black", aes(x = Smoke, y = Contribution, fill = Substrate), show.legend = F) +
#   scale_fill_manual(values = alpha(c("white", "grey25"), c(0.01, 0.99))) +
#   new_scale_fill() +
#   geom_bar_pattern(data = subset(contributions_longer, Substrate %in% "PC"), aes(x = Smoke, y = Contribution),
#                    stat = "summary", color = "black", size = 1, position = "identity", aes(pattern = Substrate), fill = 'grey25', pattern_angle = 0) +
#   new_scale_fill() +
#   # geom_point(data = subset(contributions_longer, Substrate %in% "PC"), aes(x = Smoke, y = Contribution, fill = Substrate), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   geom_beeswarm(cex = 3, data = subset(contributions_longer, Substrate %in% "PC"), aes(x = Smoke, y = Contribution, fill = Substrate), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   scale_fill_manual(values = c("white", "white"))
  





