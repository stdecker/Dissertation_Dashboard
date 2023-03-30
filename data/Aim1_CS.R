# CS Activity Analysis --------

library(readxl)
library(tidyr)
library(ggpubr)
library(ggnewscale)
library(ggprism)
library(flextable)
library(psych)
library(ggbeeswarm)
library(ARTool)
library(emmeans)
library(stringr)
library(plotrix)

#Load data ----
CS_data <- read_excel("C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Freezer Studies\\20 min incubations.xlsx", sheet = "CS Activity")

# Remove all NA
CS_data <- na.omit(CS_data)

# Long format for analysis and Tidy data
CS_data_long <- pivot_longer(CS_data, 5:12, names_to = "Group", values_to = "CS")

CS_data_long[c('Smoke', 'Tissue')] <- str_split_fixed(CS_data_long$Group, ' ', 2)

CS_data_long$Smoke <- factor(CS_data_long$Smoke, levels = unique(CS_data_long$Smoke))

CS_data_long$Tissue[CS_data_long$Tissue == "Gastroc"] <- "Gastrocnemius"

# Order variables
CS_data_long$Tissue <- factor(CS_data_long$Tissue, levels = c("Heart", "Soleus", "Gastrocnemius", "Aorta"))

CS_data_long$Date <- as.factor(CS_data_long$Date)

# Stats (ANOVA, effect size, means and SD for reported analysis)
CS_anova <- anova(art(CS ~ Smoke * Tissue + Error(Date), data = CS_data_long))

effectsize::eta_squared(CS_anova)

aggregate(CS ~ Smoke * Tissue, data = CS_data_long, FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

# Create new frame for post hoc tests
for(y in unique(CS_data_long$Tissue)){
  dunn <- subset(CS_data_long, Tissue %in% paste(y))

  assign(paste("dunn", paste(y), sep = "_"), dunn)

}

# Dunn post hoc test
Gast_dunn <- dunn.test::dunn.test(dunn_Gastrocnemius$CS, dunn_Gastrocnemius$Smoke, list = T, method = "none")
Sol_dunn <- dunn.test::dunn.test(dunn_Soleus$CS, dunn_Soleus$Smoke, list = T, method = "none")
Aorta_dunn <- dunn.test::dunn.test(dunn_Aorta$CS, dunn_Aorta$Smoke, list = T, method = "none")
Heart_dunn <- dunn.test::dunn.test(dunn_Heart$CS, dunn_Heart$Smoke, list = T, method = "none")

rstatix::cohens_d(data = subset(CS_data_long, Tissue %in% "Gastrocnemius"), formula = CS ~ Smoke, paired = T)
rstatix::cohens_d(data = subset(CS_data_long, Tissue %in% "Soleus"), formula = CS ~ Smoke, paired = T)
rstatix::cohens_d(data = subset(CS_data_long, Tissue %in% "Aorta"), formula = CS ~ Smoke, paired = T)
rstatix::cohens_d(data = subset(CS_data_long, Tissue %in% "Heart"), formula = CS ~ Smoke, paired = T)

CS_data_long$Condition <- paste(CS_data_long$Tissue, CS_data_long$Smoke)

# Plot CS Activities
CS_activity <- ggplot(data = CS_data_long, aes(x = Tissue, y = CS, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), size = 1, linetype = "longdash", color = "grey", alpha = 0.5) +
  scale_color_manual(values = c("purple", "purple4", "blue", "blue4", "green", "green4",  "red", "red4")) +
  coord_cartesian(ylim = c(0, 90), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 90, 10)) +
  labs(y=expression(bold('Citrate Synthase Activity (AU)'))) +
  geom_beeswarm(data = CS_data_long, aes(x = Tissue, y = CS, color = Condition, pch = Tissue), dodge.width = 0.99, size = 2, stroke = 1.5, show.legend = FALSE, cex = 1.75) +
  scale_shape_manual(values = c(24, 22, 21, 23)) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  annotate('text', x = 3.55, y = 89, label = 'Smoke: p = 0.007', hjust = 0, fontface = 2) +
  annotate('text', x = 3.55, y = 86, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 3.55, y = 83, label = 'Interaction: p = 0.276', hjust = 0, fontface = 2) +
  geom_bracket(xmin = .75, xmax = 1.25, y.position = 86, inherit.aes = F, label = "bold(paste('p = 0.475;', ~bolditalic(d),' = 0.25'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0) +
  geom_bracket(xmin = 1.75, xmax = 2.25, y.position = 47, inherit.aes = F, label = "bold(paste('p = 0.241;', ~bolditalic(d),' = 0.41'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0) +
  geom_bracket(xmin = 2.75, xmax = 3.25, y.position = 20, inherit.aes = F, label = "bold(paste('p = 0.024;', ~bolditalic(d),' = 1.07'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0) +
  geom_bracket(xmin = 3.75, xmax = 4.25, y.position = 70, inherit.aes = F, label = "bold(paste('p = 0.055;', ~bolditalic(d),' = 0.93'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0)

# Smoke had a significant effect on CS activity, not reflective of actual changes in mito content due to acute exposure.
# Redo analysis, but only for control tissues

effectsize::eta_squared(anova(art(CS ~ Tissue + Error(Date), data = subset(CS_data_long, Smoke %in% "Control"))))

#subset only control data
cs_sub <- subset(CS_data_long, Smoke %in% "Control")
pairwise.wilcox.test(cs_sub$CS, cs_sub$Tissue, p.adjust.method = "holm")
dunn.test::dunn.test(cs_sub$CS, cs_sub$Tissue, method = "hs", list = T)
m_CS <- lm(CS~Tissue, data = cs_sub)
es_cs <- eff_size(emmeans(m_CS, ~Tissue), sigma = sigma(m_CS), edf = df.residual(m_CS))
describeBy(cs_sub$CS, cs_sub$Tissue)

# Graph for publication
# con_CS_activity <- ggplot(data = subset(CS_data_long, Smoke %in% "Control"), aes(x = Tissue, y = CS, fill = Tissue)) +
#   stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, color = "black") +
#   theme_prism() +
#   scale_fill_manual(values = c("white", "grey33", "grey66", "black")) +
#   coord_cartesian(ylim = c(0, 90), clip = "off") +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 90, 10)) +
#   labs(y=expression(bold('Citrate Synthase Activity (AU)'))) +
#   # ggtitle("Citrate Synthase Activity") +
#   new_scale_fill() +
#   # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = CIIPercent, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   geom_beeswarm(data = subset(CS_data_long, Smoke %in% "Control"), aes(x = Tissue, y = CS, fill = Tissue, pch = Tissue), dodge.width = 0.99, size = 2, stroke = 1.5, show.legend = FALSE, cex = 1.75) +
#   scale_fill_manual(values = c("white", "white", "white", "white")) +
#   scale_shape_manual(values = c(24, 22, 21, 23)) +
#   theme(axis.title.x = element_blank(),
#         legend.position = "none",
#         axis.title.y = element_text(size = 18)) +
#   annotate('text', x = 3.55, y = 89, label = 'Main Effect: p < 0.001', hjust = 0, fontface = 2) +
#   geom_bracket(xmin = 1, xmax = 2, y.position = 86, inherit.aes = F, label = "bold(paste('p < 0.001'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0) +
#   geom_bracket(xmin = 1, xmax = 3, y.position = 81, inherit.aes = F, label = "bold(paste('p < 0.001'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0) +
#   geom_bracket(xmin = 1, xmax = 4, y.position = 76, inherit.aes = F, label = "bold(paste('p = 0.025'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0) +
#   geom_bracket(xmin = 2, xmax = 3, y.position = 35, inherit.aes = F, label = "bold(paste('p = 0.025'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0) +
#   geom_bracket(xmin = 2, xmax = 4, y.position = 62, inherit.aes = F, label = "bold(paste('p < 0.001'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0) +
#   geom_bracket(xmin = 3, xmax = 4, y.position = 57, inherit.aes = F, label = "bold(paste('p = 0.027'))", type = 'expression', size = 1, fontface = 2, tip.length = 0.0)

## Graph for Dashboard
con_CS_activity <- ggplot(data = subset(CS_data_long, Smoke %in% "Control"), aes(x = Tissue, y = CS, color = Tissue)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 0.75, fill = "white") +
  geom_beeswarm(data = subset(CS_data_long, Smoke %in% "Control"), aes(x = Tissue, y = CS, color = Tissue, pch = Tissue), dodge.width = 0.99, size = 2.5, stroke = .75, show.legend = FALSE, cex = 1.75, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("green", "red", "blue", "purple")) +
  coord_cartesian(ylim = c(0, 90), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 90, 10)) +
  labs(y=expression(bold('Citrate Synthase Activity (AU)'))) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  scale_shape_manual(values = c(24, 22, 21, 23)) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18)) +
  annotate('text', x = 3.55, y = 89, label = 'Main Effect: p < 0.001', hjust = 0, fontface = 2)


## Mean CS activities to be used for normalizing respiration

mean_gastroc <- mean(CS_data$`Control Gastroc`)
mean_soleus <- mean(CS_data$`Control Soleus`)
mean_aorta <- mean(CS_data$`Control Aorta`)
mean_heart <- mean(CS_data$`Control Heart`)

