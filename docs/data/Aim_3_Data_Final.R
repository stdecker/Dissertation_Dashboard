library(dplyr)
library(ggpubr)
library(rstatix)
library(ggnewscale)
library(scales)
library(ggprism)
library(flextable)
library(ggbeeswarm)

source("data/Aim_3_Data.R")

adp_errors <- adp_mm[-c(1,3,5,6,7,8)] |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

adp_errors$min <- adp_errors$mean - (adp_errors$Smoke == "Smoke")*adp_errors$se
adp_errors$max <- adp_errors$mean + (adp_errors$Smoke == "Control")*adp_errors$se


adp_mml$Condition <- paste(adp_mml$Tissue, adp_mml$Smoke)
adp_errors$Condition <- paste(adp_errors$Tissue, adp_errors$Smoke)

adp_kinetics <- ggplot(data = adp_mm, aes(x = Concentration, y = Rate, color = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = adp_errors, aes(x = Concentration, y = mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, show.legend = F) +
  coord_cartesian(ylim = c(0, 80), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,80,10)) +
  labs(x = expression(bold(paste("[ADP] (", mu, "M)"))), y = expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  geom_line(data = individ_adp_mml, aes(x = S, y = v, linetype = Condition), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c(1,3,1,3)) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("blue", "darkblue", "red", "darkred")) +
  geom_point(data = adp_mm, aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)


cat_errors <- cat_mm[-c(1,3,5,6,7,8)] |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs = funs(mean, se=sd(.)/sqrt(n())))

cat_errors$min <- cat_errors$mean + (cat_errors$Smoke == "Smoke")*cat_errors$se
cat_errors$max <- cat_errors$mean + (cat_errors$Smoke == "Control")*cat_errors$se

cat_errors$Condition <- paste(cat_errors$Tissue, cat_errors$Smoke)

cat_kinetics <- ggplot(data = cat_mm, aes(x = Concentration, y = Rate, color = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = cat_errors, aes(x = Concentration, y = mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, fill = "white", show.legend = F) +
  coord_cartesian(ylim = c(0, 120), xlim = c(-0.25, 5.25), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 120, 20)) +
  labs(x = expression(bold(paste("[CAT] (", mu, "M)"))), y = expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  geom_line(data = individ_cat_mml, aes(x = S, y = v, linetype = Condition), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c(1,3,1,3)) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("blue", "darkblue", "red", "darkred")) +
  geom_point(data = cat_mm, aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)


# Graphs ----
## ADP ----
# identify_outliers(subset(subset(individ_adp_kinetics, Smoke %in% "Control"), Tissue %in% "Gastrocnemius"), Vmax)

adp_km <- ggplot(data = individ_adp_kinetics, aes(x = Condition, y = Km, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 900)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 900, 100)) +
  labs(y= expression(bold(paste("[ADP] (", mu, "M)")))) +
  # ggtitle(expression(bold(K[m]))) + 
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_adp_kinetics, aes(x = Condition, y = Km, color = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 880, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 840, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 800, label = 'Interaction: p < 0.001', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 790, inherit.aes = F, label = "p = 0.042", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Soleus Control", y.position = 660, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Soleus Control", y.position = 730, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.01)

adp_vmax <- ggplot(data = individ_adp_kinetics, aes(x = Condition, y = Vmax, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 160)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 160, 20)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  # ggtitle(expression(bold(V[max]))) +
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Vmax, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_adp_kinetics, aes(x = Condition, y = Vmax, color = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 157, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 150, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 143, label = 'Interaction: p = 0.813', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 145, inherit.aes = F, label = "p = 0.013", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Soleus Control", y.position = 120, inherit.aes = F, label = "p = 0.027", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Soleus Smoke", y.position = 133, inherit.aes = F, label = "p = 0.031", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Gastrocnemius Control", y.position = 66, inherit.aes = F, label = "p = 0.014", size = 1, fontface = 2, tip.length = 0.02)

## CAT ----

cat_imax <- ggplot(data = individ_cat_kinetics[-c(28,36),], aes(x = Condition, y = Imax, color =  Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 80)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 80, 10)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  # ggtitle(expression(bold(I[max]))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = Imax, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_cat_kinetics[-c(28,36),], aes(x = Condition, y = Imax, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 79, label = 'Smoke: p = 0.056', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 75, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 71, label = 'Interaction: p = 0.506', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 75, inherit.aes = F, label = "p = 0.368", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = 50, inherit.aes = F, label = "p = 0.236", size = 1, fontface = 2, tip.length = 0.02) #+
  # geom_bracket(xmin = "Gastrocnemius Control", xmax = "Soleus Control", y.position = 102, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.02) +
  # geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Soleus Smoke", y.position = 117, inherit.aes = F, label = "p = 0.008", size = 1, fontface = 2, tip.length = 0.02) 


cat_km <- ggplot(data = individ_cat_kinetics[-c(36, 7),], aes(x = Condition, y = IC50, color =  Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 5, 0.5)) +
  labs(y= expression(bold('[CAT] (μM)'))) +
  # ggtitle(expression(bold(IC['50']))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = IC50, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_cat_kinetics[-c(36, 7),], aes(x = Condition, y = IC50, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 4.9, label = 'Smoke: p = 0.679', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 4.7, label = 'Tissue: p = 0.481', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 4.5, label = 'Interaction: p = 0.541', hjust = 0, fontface = 2)

cat_percent <- ggplot(data = individ_cat_kinetics[-c(7, 24),], aes(x = Condition, y = Inhibition, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(-1.2, 0)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1.3, 0, 0.1), labels = scales::percent_format()) +
  labs(y=expression(bold(Inhibition~('%'~V[max])))) +
  # ggtitle(expression(bold(CAT~Inhibition))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = Inhibition, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_cat_kinetics[-c(7, 24),], aes(x = Condition, y = Inhibition, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_x_discrete(position = "top") +
  annotate('text', x = 4.45, y = -1, label = 'Smoke: p < 0.001', hjust = 1, fontface = 2) +
  annotate('text', x = 4.45, y = -1.05, label = 'Tissue: p = 0.048', hjust = 1, fontface = 2) +
  annotate('text', x = 4.45, y = -1.1, label = 'Interaction: p = 0.004', hjust = 1, fontface = 2) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = -1.02, inherit.aes = F, label = "p < 0.001", size = 1, fontface = 2, tip.length = -0.02, vjust = 2.2) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = -.8, inherit.aes = F, label = "p = 0.077", size = 1, fontface = 2, tip.length = -0.02, vjust = 2.2) 


cat_vmax <- ggplot(data = individ_cat_kinetics, aes(x = Condition, y = C, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 200), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 200, 25)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  # ggtitle(expression(bold(Initial~Rate~"(C)"))) +
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = C, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_cat_kinetics, aes(x = Condition, y = C, color = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 198, label = 'Smoke: p = 0.004', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 190, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 182, label = 'Interaction: p = 0.717', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 166, inherit.aes = F, label = "p = 0.073", size = 1, fontface = 2, tip.length = 0.02) +
  # geom_bracket(xmin = "Gastrocnemius Control", xmax = "Soleus Control", y.position = 164, inherit.aes = F, label = "p = 0.006", size = 1, fontface = 2, tip.length = 0.02) +
  # geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Soleus Smoke", y.position = 177, inherit.aes = F, label = "p = 0.008", size = 1, fontface = 2, tip.length = 0.02) #+
  geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Gastrocnemius Control", y.position = 82, inherit.aes = F, label = "p = 0.113", size = 1, fontface = 2, tip.length = 0.02)

## Main Plot ----
mainplot_all <- ggplot(data = subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
                                                            , 'Rot')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  #geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 200), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 200, 25)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), linetype = "longdash", color = "grey50", size = .6) + 
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
                                                  , 'Rot')), 
                aes(x = State, y = Rate, color = Condition), 
             dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  scale_x_discrete(labels = c(expression(bold(GM)), expression(bold(GMD['5000 μM'])), expression(bold(GMDS)), expression(bold(CAT['5.0 μM'])), expression(bold(FCCP[Peak])), expression(bold(Rot)), expression(bold(AmA~'&'~Omy)))) +
  #Basal
  # annotate('text', x = 0.45, y = 198, label = 'Smoke: p = 0.408', hjust = 0, fontface = 2) +
  # annotate('text', x = 0.45, y = 191, label = 'Tissue: p = 0.283', hjust = 0, fontface = 2) +
  # annotate('text', x = 0.45, y = 184, label = 'Interaction: p = 0.829', hjust = 0, fontface = 2) +
  #GM
  annotate('text', x = .51, y = 198, label = 'Smoke: p = 0.150', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 184, label = 'Interaction: p = 0.175', hjust = 0, fontface = 2) +
  # geom_bracket(xmin = c(1.665), xmax = c(2.115), y.position = 55, inherit.aes = F, label = "p = 0.048", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(1.885), xmax = c(2.345), y.position = 70, inherit.aes = F, label = "p = 0.019", size = 1, fontface = 2, tip.length = 0.01) +
  #GMD
  annotate('text', x = 1.51, y = 198, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 191, label = 'Tissue: p = 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 184, label = 'Interaction: p = 0.604', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(1.625), xmax = c(1.885), y.position = 70, inherit.aes = F, label = "p = 0.014", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(2.115), xmax = c(2.35), y.position = 125, inherit.aes = F, label = "p = 0.038", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(1.665), xmax = c(2.115), y.position = 125, inherit.aes = F, label = "p = 0.048", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(1.885), xmax = c(2.345), y.position = 135, inherit.aes = F, label = "p = 0.019", size = 1, fontface = 2, tip.length = 0.01) +
  #GMDS
  annotate('text', x = 2.51, y = 198, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 184, label = 'Interaction: p = 0.604', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(2.625), xmax = c(2.885), y.position = 75, inherit.aes = F, label = "p = 0.052", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(3.115), xmax = c(3.345), y.position = 160, inherit.aes = F, label = "p = 0.080", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(2.665), xmax = c(3.115), y.position = 150, inherit.aes = F, label = "p = 0.012", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(2.885), xmax = c(3.345), y.position = 125, inherit.aes = F, label = "p = 0.008", size = 1, fontface = 2, tip.length = 0.01) +
  #CAT
  annotate('text', x = 3.51, y = 198, label = 'Smoke: p = 1.000', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 184, label = 'Interaction: p = 0.814', hjust = 0, fontface = 2) +
  # geom_bracket(xmin = c(3.665), xmax = c(4.115), y.position = 130, inherit.aes = F, label = "p < 0.001", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(3.885), xmax = c(4.345), y.position = 140, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.01) +
  #FCCP
  annotate('text', x = 4.51, y = 198, label = 'Smoke: p = 0.004', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 184, label = 'Interaction: p = 0.857', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(4.625), xmax = c(4.885), y.position = 95, inherit.aes = F, label = "p = 0.122", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(5.115), xmax = c(5.345), y.position = 165, inherit.aes = F, label = "p = 0.101", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(4.665), xmax = c(5.115), y.position = 162, inherit.aes = F, label = "p = 0.008", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(4.885), xmax = c(5.345), y.position = 172, inherit.aes = F, label = "p = 0.015", size = 1, fontface = 2, tip.length = 0.01) +
  #Rot
  annotate('text', x = 5.51, y = 198, label = 'Smoke: p = 0.022', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 184, label = 'Interaction: p = 0.468', hjust = 0, fontface = 2)# +
  # geom_bracket(xmin = c(5.665), xmax = c(6.115), y.position = 100, inherit.aes = F, label = "p = 0.002", size = 1, fontface = 2, tip.length = 0.01) +
  # geom_bracket(xmin = c(5.885), xmax = c(6.345), y.position = 110, inherit.aes = F, label = "p = 0.015", size = 1, fontface = 2, tip.length = 0.01) +
  #Ama & Omy
  # annotate('text', x = 6.51, y = 198, label = 'Smoke: p = 0.072', hjust = 0, fontface = 2) +
  # annotate('text', x = 6.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  # annotate('text', x = 6.51, y = 184, label = 'Interaction: p = 0.535', hjust = 0, fontface = 2) #+
# geom_bracket(xmin = c(7.665), xmax = c(8.115), y.position = 60, inherit.aes = F, label = "p = 0.008", size = 1, fontface = 2, tip.length = 0.01) +
# geom_bracket(xmin = c(7.885), xmax = c(8.345), y.position = 70, inherit.aes = F, label = "p = 0.010", size = 1, fontface = 2, tip.length = 0.01)

# RCR ----
cat_RCR <- ggplot(data = cat_data, aes(x = Condition, y = RCR, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 8), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 8, 1)) +
  labs(y=expression(bold(RCR~('GMDS/GM')))) +
  ggtitle("Respiratory Control Ratio") +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = Inhibition, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = cat_data, aes(x = Condition, y = RCR, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', y = 7.9, x = .55, label = 'Smoke: p = 0.075', hjust = 0, fontface = 2) +
  annotate('text', y = 7.55, x = .55, label = 'Tissue: p = 0.754', hjust = 0, fontface = 2) +
  annotate('text', y = 7.2, x = .55, label = 'Interaction: p = 0.835', hjust = 0, fontface = 2)

RCR_anova <- anova(art(data = cat_data, RCR ~ Smoke * Tissue))
effectsize::eta_squared(RCR_anova)

# RCR ----
cat_Thermodynamic_Coupling <- ggplot(data = cat_data[-c(34),], aes(x = Condition, y = Thermodynamic_Coupling, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 1.1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.1, .1)) +
  labs(y=expression(bold('q-value'))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = Inhibition, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = cat_data[-c(34),], aes(x = Condition, y = Thermodynamic_Coupling, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = .55, y = 1.09, label = 'Smoke: p = 0.071', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = 1.04, label = 'Tissue: p = 0.544', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = .99, label = 'Interaction: p = 0.235', hjust = 0, fontface = 2)

# levene_test(data = cat_data[-c(34),], Thermodynamic_Coupling ~ Smoke * Tissue)
# 
# as.data.table(cat_data[-c(34),])[,.(Statistic = shapiro.test(Thermodynamic_Coupling)$statistic,
#                                     P.value = shapiro.test(Thermodynamic_Coupling)$p.value),
#                                 by = .(Tissue, Smoke)]

TC_anova <- anova(art(data = cat_data[-c(27, 24),], Thermodynamic_Coupling ~ Smoke * Tissue))
effectsize::eta_squared(TC_anova)


## Specific Rate

cat_Specific_rate <- ggplot(data = cat_data[-c(25, 27),], aes(x = Condition, y = SpecificRate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 1.1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.1, .1)) +
  labs(y=expression(bold('Complex I/II Specific Rate (GM/Rot)'))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = SpecificRate, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = cat_data[-c(25, 27),], aes(x = Condition, y = SpecificRate, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = .55, y = 1.09, label = 'Smoke: p = 0.071', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = 1.04, label = 'Tissue: p = 0.544', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = .99, label = 'Interaction: p = 0.235', hjust = 0, fontface = 2)

# levene_test(data = cat_data[-c(34),], Thermodynamic_Coupling ~ Smoke * Tissue)
# 
# as.data.table(cat_data[-c(34),])[,.(Statistic = shapiro.test(Thermodynamic_Coupling)$statistic,
#                                     P.value = shapiro.test(Thermodynamic_Coupling)$p.value),
#                                 by = .(Tissue, Smoke)]


SpecificRate_anova <- anova(art(data = cat_data[-c(25,27, 41:44),], SpecificRate ~ Smoke * Tissue))
effectsize::eta_squared(SpecificRate_anova)

## Rot/FCCP

cat_RotFCCP <- ggplot(data = cat_data[-c(9),], aes(x = Condition, y = CIIPercent, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 1.1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.1, .1), labels = label_percent()) +
  labs(y=expression(bold(Complex~II~Contribution~('Rot/FCCP'['Peak'])))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = CIIPercent, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = cat_data[-c(9),], aes(x = Condition, y = CIIPercent, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = .55, y = 1.09, label = 'Smoke: p = 0.015', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = 1.04, label = 'Tissue: p = 0.466', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = .99, label = 'Interaction: p = 0.048', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = 0.88, inherit.aes = F, label = "p = 0.006", size = 1, fontface = 2, tip.length = 0.02)

# levene_test(data = cat_data[-c(34),], Thermodynamic_Coupling ~ Smoke * Tissue)
# 
# as.data.table(cat_data[-c(34),])[,.(Statistic = shapiro.test(Thermodynamic_Coupling)$statistic,
#                                     P.value = shapiro.test(Thermodynamic_Coupling)$p.value),
#                                 by = .(Tissue, Smoke)]

identify_outliers(subset(cat_data, Condition %in% "Soleus Smoke"), CIIPercent)

RotFCCP_anova <- anova(art(data = cat_data[-c(9, 41:44),], CIIPercent ~ Smoke * Tissue))
effectsize::eta_squared(RotFCCP_anova)

dunn.test::dunn.test(cat_data$CIIPercent, cat_data$Condition, method = "hs", list = T)
rstatix::cohens_d(data = subset(cat_data[-c(9),], Tissue %in% "Gastrocnemius"), formula = CIIPercent ~ Smoke)
rstatix::cohens_d(data = subset(cat_data[-c(9),], Tissue %in% "Soleus"), formula = CIIPercent ~ Smoke)
rstatix::cohens_d(data = subset(cat_data[-c(9),], Smoke %in% "Control"), formula = CIIPercent ~ Tissue)
rstatix::cohens_d(data = subset(cat_data[-c(9),], Smoke %in% "Smoke"), formula = CIIPercent ~ Tissue)
rstatix::cohens_d(data = subset(cat_data[-c(9),], Condition %in% c("Gastrocnemius Smoke", "Soleus Control")), formula = CIIPercent ~ Condition)
rstatix::cohens_d(data = subset(cat_data[-c(9),], Condition %in% c("Gastrocnemius Control", "Soleus Smoke")), formula = CIIPercent ~ Condition)

## FCCP/CAT

cat_FCCPCAT <- ggplot(data = cat_data[-c(24),], aes(x = Condition, y = FCCP_per_CAT, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 3), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 3, .5)#, labels = label_percent()
                     ) +
  labs(y=expression(bold(FCCP['Peak']*'/'*CAT['5.0']))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = CIIPercent, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = cat_data[-c(24),], aes(x = Condition, y = FCCP_per_CAT, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 2.55, y = 2.97, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 2.55, y = 2.85, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 2.55, y = 2.73, label = 'Interaction: p = 0.001', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = 2.83, inherit.aes = F, label = "p = 0.007", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 2.43, inherit.aes = F, label = "p = 0.047", size = 1, fontface = 2, tip.length = 0.02)

# levene_test(data = cat_data[-c(34),], Thermodynamic_Coupling ~ Smoke * Tissue)
# 
# as.data.table(cat_data[-c(34),])[,.(Statistic = shapiro.test(Thermodynamic_Coupling)$statistic,
#                                     P.value = shapiro.test(Thermodynamic_Coupling)$p.value),
#                                 by = .(Tissue, Smoke)]

identify_outliers(subset(cat_data, Condition %in% "Gastrocnemius Control"), FCCP_per_CAT)

FCCPCAT_anova <- anova(art(data = cat_data[-c(24),], FCCP_per_CAT ~ Smoke * Tissue))
effectsize::eta_squared(FCCPCAT_anova)

dunn.test::dunn.test(cat_data$FCCP_per_CAT[-c(24)], cat_data$Condition[-c(24)], method = "hs", list = T)
rstatix::cohens_d(data = subset(cat_data, Tissue %in% "Gastrocnemius"), formula = FCCP_per_CAT ~ Smoke)
rstatix::cohens_d(data = subset(cat_data, Tissue %in% "Soleus"), formula = FCCP_per_CAT ~ Smoke)
rstatix::cohens_d(data = subset(cat_data, Smoke %in% "Control"), formula = FCCP_per_CAT ~ Tissue)
rstatix::cohens_d(data = subset(cat_data, Smoke %in% "Smoke"), formula = FCCP_per_CAT ~ Tissue)
rstatix::cohens_d(data = subset(cat_data, Condition %in% c("Gastrocnemius Smoke", "Soleus Control")), formula = FCCP_per_CAT ~ Condition)
rstatix::cohens_d(data = subset(cat_data, Condition %in% c("Gastrocnemius Control", "Soleus Smoke")), formula = FCCP_per_CAT ~ Condition)

## GMDS/CAT

cat_GMDSCAT <- ggplot(data = cat_data[-c(9),], aes(x = Condition, y = GMDS_per_CAT, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 3.5), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 3.5, .5)#, labels = label_percent()
  ) +
  labs(y=expression(bold(GMDS*'/'*CAT['5.0']))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = CIIPercent, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = cat_data[-c(9),], aes(x = Condition, y = GMDS_per_CAT, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 2.55, y = 3.45, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 2.55, y = 3.32, label = 'Tissue: p = 0.004', hjust = 0, fontface = 2) +
  annotate('text', x = 2.55, y = 3.19, label = 'Interaction: p = 0.003', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = 3.37, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 2.5, inherit.aes = F, label = "p = 0.076", size = 1, fontface = 2, tip.length = 0.02)

# levene_test(data = cat_data[-c(34),], Thermodynamic_Coupling ~ Smoke * Tissue)
# 
# as.data.table(cat_data[-c(34),])[,.(Statistic = shapiro.test(Thermodynamic_Coupling)$statistic,
#                                     P.value = shapiro.test(Thermodynamic_Coupling)$p.value),
#                                 by = .(Tissue, Smoke)]

identify_outliers(subset(cat_data, Condition %in% "Soleus Smoke"), GMDS_per_CAT)

GMDSCAT_anova <- anova(art(data = cat_data[-c(9),], GMDS_per_CAT ~ Smoke * Tissue))
effectsize::eta_squared(GMDSCAT_anova)

dunn.test::dunn.test(cat_data$GMDS_per_CAT[-c(9)], cat_data$Condition[-c(9)], method = "hs", list = T)
rstatix::cohens_d(data = subset(cat_data, Tissue %in% "Gastrocnemius"), formula = GMDS_per_CAT ~ Smoke)
rstatix::cohens_d(data = subset(cat_data, Tissue %in% "Soleus"), formula = GMDS_per_CAT ~ Smoke)
rstatix::cohens_d(data = subset(cat_data, Smoke %in% "Control"), formula = GMDS_per_CAT ~ Tissue)
rstatix::cohens_d(data = subset(cat_data, Smoke %in% "Smoke"), formula = GMDS_per_CAT ~ Tissue)
rstatix::cohens_d(data = subset(cat_data, Condition %in% c("Gastrocnemius Smoke", "Soleus Control")), formula = GMDS_per_CAT ~ Condition)
rstatix::cohens_d(data = subset(cat_data, Condition %in% c("Gastrocnemius Control", "Soleus Smoke")), formula = GMDS_per_CAT ~ Condition)

aggregate(GMDS_per_CAT ~ Smoke * Tissue, data = cat_data[-c(9),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

## GMDS/FCCP

cat_GMDSFCCP <- ggplot(data = cat_data, aes(x = Condition, y = GMDS_per_FCCP, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 1.75), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.75, .25)#, labels = label_percent()
  ) +
  labs(y=expression(bold(GMDS*'/'*FCCP[Peak]))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = CIIPercent, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = cat_data[-c(9),], aes(x = Condition, y = GMDS_per_FCCP, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = .55, y = 1.73, label = 'Smoke: p = 0.210', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = 1.66, label = 'Tissue: p = 0.569', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = 1.59, label = 'Interaction: p = 0.349', hjust = 0, fontface = 2) #+
  # geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = 3.37, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.02) +
  # geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 2.5, inherit.aes = F, label = "p = 0.075", size = 1, fontface = 2, tip.length = 0.02)

# levene_test(data = cat_data[-c(34),], Thermodynamic_Coupling ~ Smoke * Tissue)
# 
# as.data.table(cat_data[-c(34),])[,.(Statistic = shapiro.test(Thermodynamic_Coupling)$statistic,
#                                     P.value = shapiro.test(Thermodynamic_Coupling)$p.value),
#                                 by = .(Tissue, Smoke)]

identify_outliers(subset(cat_data, Condition %in% "Soleus Control"), GMDS_per_FCCP)

GMDSFCCP_anova <- anova(art(data = cat_data, GMDS_per_FCCP ~ Smoke * Tissue))
effectsize::eta_squared(GMDSFCCP_anova)

dunn.test::dunn.test(cat_data$GMDS_per_FCCP, cat_data$Condition, method = "hs", list = T)
rstatix::cohens_d(data = subset(cat_data, Tissue %in% "Gastrocnemius"), formula = GMDS_per_FCCP ~ Smoke)
rstatix::cohens_d(data = subset(cat_data, Tissue %in% "Soleus"), formula = GMDS_per_FCCP ~ Smoke)
rstatix::cohens_d(data = subset(cat_data, Smoke %in% "Control"), formula = GMDS_per_FCCP ~ Tissue)
rstatix::cohens_d(data = subset(cat_data, Smoke %in% "Smoke"), formula = GMDS_per_FCCP ~ Tissue)
rstatix::cohens_d(data = subset(cat_data, Condition %in% c("Gastrocnemius Smoke", "Soleus Control")), formula = GMDS_per_FCCP ~ Condition)
rstatix::cohens_d(data = subset(cat_data, Condition %in% c("Gastrocnemius Control", "Soleus Smoke")), formula = GMDS_per_FCCP ~ Condition)

## CAT Contribution

cat_CAT_contribution <- ggplot(data = cat_data, aes(x = Condition, y = CAT_cont, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("blue", "darkblue", "red", "darkred")) +
  coord_cartesian(ylim = c(0, 50), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 50, 5)#, labels = label_percent()
  ) +
  labs(y=expression(bold(Delta~GMDS*'-'*CAT['5.0']))) +
  ggtitle(expression(bold(Delta~GMDS*'-'*CAT['5.0']))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = CIIPercent, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = cat_data[-c(9),], aes(x = Condition, y = CAT_cont, color =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = .55, y = 48, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = 45, label = 'Tissue: p = 0.938', hjust = 0, fontface = 2) +
  annotate('text', x = .55, y = 42, label = 'Interaction: p = 0.716', hjust = 0, fontface = 2)# +
  # geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = 3.37, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.02) +
  # geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 2.5, inherit.aes = F, label = "p = 0.075", size = 1, fontface = 2, tip.length = 0.02)

# levene_test(data = cat_data[-c(34),], Thermodynamic_Coupling ~ Smoke * Tissue)
# 
# as.data.table(cat_data[-c(34),])[,.(Statistic = shapiro.test(Thermodynamic_Coupling)$statistic,
#                                     P.value = shapiro.test(Thermodynamic_Coupling)$p.value),
#                                 by = .(Tissue, Smoke)]

identify_outliers(subset(cat_data, Condition %in% "Gastrocnemius Smoke"), CAT_cont)

CAT_cont_anova <- anova(art(data = cat_data, CAT_cont ~ Smoke * Tissue))
effectsize::eta_squared(CAT_cont_anova)

dunn.test::dunn.test(cat_data$CAT_cont, cat_data$Condition, method = "hs", list = T)
rstatix::cohens_d(data = subset(cat_data, Tissue %in% "Gastrocnemius"), formula = CAT_cont ~ Smoke)
rstatix::cohens_d(data = subset(cat_data, Tissue %in% "Soleus"), formula = CAT_cont ~ Smoke)
rstatix::cohens_d(data = subset(cat_data, Smoke %in% "Control"), formula = CAT_cont ~ Tissue)
rstatix::cohens_d(data = subset(cat_data, Smoke %in% "Smoke"), formula = CAT_cont ~ Tissue)
rstatix::cohens_d(data = subset(cat_data, Condition %in% c("Gastrocnemius Smoke", "Soleus Control")), formula = CAT_cont ~ Condition)
rstatix::cohens_d(data = subset(cat_data, Condition %in% c("Gastrocnemius Control", "Soleus Smoke")), formula = CAT_cont ~ Condition)

aggregate(CIIPercent ~ Smoke * Tissue, data = cat_data, FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

