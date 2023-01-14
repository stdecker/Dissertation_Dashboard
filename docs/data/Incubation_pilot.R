library(readxl)
library(dplyr)
library(tidyr)
library(ggpubr)
library(reshape2)
library(data.table)
library(rstatix)
library(plyr)
library(drc)
library(ggnewscale)
library(scales)
library(RColorBrewer)
library(EnvStats)
library(ggprism)
library(kableExtra)
library(tidyverse)
library(rvg)
library(flextable)
library(officer)
library(officedown)
library(knitr)
library(gridExtra)
library(ggbrace)
library(reshape2)
library(segmented)


#Load data ----
data2 <- read_excel('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Project 2\\Analysis File2test.xlsx', sheet = "Final Data", col_names = TRUE)


data2 <- data2[-c(10, 11, 12:15)]

#colnames(data2)[colnames(data2) == "GMDS...9"] <- "GMDS"

data2$Incubation <- factor(data2$Incubation, levels = c("CON", "1-hour", "3-hour", "6-hour"))

data2_long <- pivot_longer(data2, 6:9, names_to = "State", values_to = "Rate")

data2_long <- na.omit(data2_long)

IncubationPlot <- ggplot(subset(data2_long, State %in% c("GMDS")), aes(x = Incubation, y = Rate, fill = Incubation)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge2(width = 0.9), color = "black", size = 1, na.rm = TRUE) +
  geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,150), breaks = seq(0,150,10)) +
  ggtitle("State III Respiration Under Different Incubations in Gastroc") +
  ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  xlab("") +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11)
  ) +
  stat_compare_means(label.y = 140) +
  stat_compare_means(ref.group = "CON", label.y = 130, label = "p.format") +
  scale_fill_grey(start = 1, end = 0.25)

## Graph full protocol
fulltest <- subset(data2, Incubation != "CON")

full_long <- pivot_longer(fulltest, 6:9, names_to = "State", values_to = "Rate")

FullTest_Plot <- ggplot(data = full_long, aes(x = State, y = Rate, fill = Incubation)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge2(width = 0.9), color = "black", size = 1) +
  geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,120), breaks = seq(0,120,10)) +
  ggtitle("Respiration Rates Under Different Incubations in Gastroc") +
  ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  xlab("") +
  theme(
    plot.title = element_text(size=11)
  ) +
  scale_fill_grey(start = 1, end = 0.25) +
  guides(fill = guide_legend(override.aes = list(shape = NA)))


## Percents ----
data2 <- data2[!is.na(data2$Subject),]

data2 <- separate(data2, Subject, "date", " ")

melted <- na.omit(pivot_longer(data2, 6:9))

melted <- subset(melted, name == c("GMDS"))

melted <- melted[-c(2,9),]

percent <- pivot_wider(melted, date, names_from = Incubation, values_from = value)

percent$`1-Hour` <- (percent$`1-hour`/percent$CON)-1

percent$`3-Hour` <- (percent$`3-hour`/percent$CON)-1

percent$`6-Hour` <- (percent$`6-hour`/percent$CON)-1


percent_long <- pivot_longer(percent, c(2,6:8))


Percent_Incubation <- ggplot(subset(percent_long, name %in% c(paste(names(percent[6:8])))), aes(x = name, y = value, fill = name)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge2(width = 0.9), color = "black", size = 1, na.rm = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0), limits = c(-1,0), breaks = seq(-1,0,.10)) +
  ggtitle("Percent Inhibition Under Incubations (Gastroc)") +
  ylab(bquote(atop(bold("Inhibition of State III Respiration"), bold("relative to CON")))) +
  theme_prism() +
  xlab("") +
  theme(
    plot.title = element_text(size=11),
    legend.position =  "none"
  ) +
  scale_x_discrete(position = "top") +
  stat_summary(aes(label = round(..y.., 2)), fun.y=mean, geom = "text", vjust = 7.5) +
  scale_fill_grey(start = 1, end = 0) +
  new_scale_fill() +
  geom_jitter(data = subset(percent_long, name %in% c(paste(names(percent[6:8])))), aes(x = name, y = value, fill = name), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), shape = 21, size = 2, stroke = 1.5) +
  scale_fill_manual(values = c("white", "white", "white")) +
  guides(fill = guide_legend(override.aes = list(shape = NA)))


# Percent 2 ----
percent_long2 <- percent_long[-c(6),]

Percent_Incubation2 <- ggplot(subset(percent_long2, name %in% c(paste(names(percent[6:8])))), aes(x = name, y = value, fill = name)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge2(width = 0.9), color = "black", size = 1, na.rm = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0), limits = c(-1,0), breaks = seq(-1,0,.10)) +
  ggtitle("Percent Inhibition Under Incubations (Gastroc)") +
  ylab(bquote(atop(bold("Inhibition of State III Respiration"), bold("relative to CON")))) +
  theme_prism() +
  xlab("") +
  theme(
    plot.title = element_text(size=11),
    legend.position =  "none"
  ) +
  scale_x_discrete(position = "top") +
  geom_point() +
  stat_summary(aes(label = round(..y.., 2)), fun.y=mean, geom = "text", vjust = 7.5) +
  scale_fill_grey(start = 1, end = 0.25)


# Soleus ----
data3 <- read_excel('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Project 2\\Analysis File2Soleus.xlsx', sheet = "Final Data", col_names = TRUE)


data3 <- data3[-c(10:15)]

#colnames(data2)[colnames(data2) == "GMDS...9"] <- "GMDS"

data3$Incubation <- factor(data3$Incubation, levels = c("CON", "1-Hour", "3-Hour", "6-Hour"))

data3_long <- pivot_longer(data3, 6:9, names_to = "State", values_to = "Rate")

data3_long <- na.omit(data3_long)

SolIncubationPlot <- ggplot(subset(data3_long, State %in% c("GMDS")), aes(x = Incubation, y = Rate, fill = Incubation)) +
  stat_summary(fun.data = "mean_se",geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge2(width = 0.9), color = "black", size = 1, na.rm = TRUE) +
  geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,250), breaks = seq(0,250,25)) +
  ggtitle("State III Respiration Under Different Incubations in Soleus") +
  ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  xlab("") +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11)
  ) +
  stat_compare_means(label.y = 240) +
  stat_compare_means(ref.group = "CON", label.y = 230, label = "p.format") +
  scale_fill_grey(start = 1, end = 0.25)

# Soleus Percent ----
Sol_inc_percent <- aggregate(data3$GMDS, list(data3$Incubation), mean)
Sol_inc_percent <- pivot_wider(Sol_inc_percent, names_from = Group.1, values_from = x)

Sol_inc_percent$`1-Hour Percent` <- Sol_inc_percent$`1-Hour`/Sol_inc_percent$CON
Sol_inc_percent$`3-Hour Percent` <- Sol_inc_percent$`3-Hour`/Sol_inc_percent$CON
Sol_inc_percent$`6-Hour Percent` <- Sol_inc_percent$`6-Hour`/Sol_inc_percent$CON

graph_sol <- melt(Sol_inc_percent[5:7])



#AGAIN!!  Load data ----
INCdata <- read_excel('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Project 2\\Analysis File2.xlsx', sheet = "Final Data", col_names = TRUE)

INCdata <- INCdata[-c(11:15)]

INCdata <- na.omit(INCdata)

INCdata$RCR <- INCdata$GMDS/INCdata$GM

INCdata$Concentration <- factor(INCdata$Concentration, levels = c("0", "0.004", "0.04", "0.4"))

INCdata <- separate(INCdata, Subject, "Subject", " ")

INCdata_long <- pivot_longer(INCdata, 7:10, names_to = "State", values_to = "Rate")



#INCdata_long$Concentration <- factor(INCdata_long$Concentration, levels = c("0", "0.004", "0.04", "0.4"))

data2_long <- separate(data2_long, Subject, "Subject", " ")

data2_long$Incubation <- factor(data2_long$Incubation, levels = unique(data2_long$Incubation))

summary(aov(Rate ~ Incubation * State + Error(factor(Subject)), data2_long, type = "III"))

dependent <- as.list(data2[6:9])

output <- vector("list")
for(y in names(dependent)){
  output[y] <- summary(aov(dependent[[y]] ~ data2$Incubation))
  cat(paste('\nDependent var:', y, '\n'))
}

ANOVA_pvalues <- data.frame("Pr(>F)" = sapply(output, getElement, name = "Pr(>F)"))

# anova_test(data = data2_long, dv = Rate, wid = Subject, within = c(Incubation, State))

inc_sol_percent <- ggplot(graph_sol, aes(x = variable, y = value, fill = variable)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge2(width = 0.9), color = "black", size = 1, na.rm = TRUE) +
  scale_fill_grey(start = 0.8, end = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),limits = c(0,1), expand = c(0,0)) +
  ylab("Respiration Rate (Relative to CON)") +
  xlab("") +
  theme(legend.position = "none") +
  ggtitle("State III Respiration of Soleus Relative to CON") +
  theme_prism()

incubation_time_plot <- ggplot(data = INCdata_long, aes(x = Time, y = Rate, fill = Concentration)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, color = "black", size = 1) +
  geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,10)) +
  ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  xlab("Incubation Time") +
  theme(
    plot.title = element_text(size=11),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 12)
  ) +
  # scale_fill_grey(start = 1, end = 0.25, name = "Concentration (μg/mL)") +
  scale_color_prism() +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  facet_grid(~State)

incubation_concentration_plot <- ggplot(data = INCdata_long, aes(x = Concentration, y = Rate, fill = Time)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, color = "black", size = 1) +
  geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,10)) +
  ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  xlab("Concentration (μg/mL)") +
  theme(
    plot.title = element_text(size=11),
    strip.background = element_blank()
  ) +
  scale_fill_brewer(palette = "Paired") +
  #scale_fill_grey(start = 1, end = 0.25, name = "Incubation Time") +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  facet_grid(~State) 

incubation_rcr_plot <- ggplot(INCdata, aes(x = Time, y = RCR, fill = Concentration)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, color = "black", size = 1) +
  geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,10), breaks = seq(0,10,1)) +
  ylab("Respiratory Control Ratio") +
  theme_prism() +
  xlab("Incubation Time") +
  theme(
    plot.title = element_text(size=11)
  ) +
  scale_fill_brewer(palette = "Paired") +
  #scale_fill_grey(start = 1, end = 0.25, name = "Concentration (μg/mL)") +
  guides(fill = guide_legend(override.aes = list(shape = NA)))

incubation_GMDS_plot <- ggplot(data = subset(INCdata_long, State %in% "GMDS"), aes(x = Concentration, y = Rate, fill = Time)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, color = "black", size = 1) +
  geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,10)) +
  ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  xlab("Concentration (μg/mL)") +
  theme(
    plot.title = element_text(size=11),
    strip.background = element_blank()
  ) +
  scale_fill_brewer(palette = "Paired") +
  #scale_fill_grey(start = 1, end = 0.25, name = "Incubation Time") +
  guides(fill = guide_legend(override.aes = list(shape = NA)))

Highdata <- read_excel('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Project 2\\Analysis File2High.xlsx', sheet = "Final Data", col_names = TRUE)

Highdata <- Highdata[-c(13:16)]

Highdata <- na.omit(Highdata)

Highdata$RCR <- Highdata$GMDS/Highdata$GM

Highdata$Concentration <- factor(Highdata$Concentration, levels = c("0", "400", "800"))

Highdata <- separate(Highdata, Subject, "Subject", " ")

Highdata_long <- pivot_longer(Highdata, 8:11, names_to = "State", values_to = "Rate")



for(i in unique(Highdata$Time)){
  for(y in unique(Highdata$Tissue)){
    highplot <- ggplot(data = subset(subset(Highdata_long, Tissue %in% paste(y)), Time %in% paste(i)), aes(x = Concentration, y = Rate, fill = Concentration)) +
      stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
      stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, color = "black", size = 1) +
      geom_point(aes(shape = Subject), size = 1.75) +
      geom_path(aes(group = Subject)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,10)) +
      ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
      theme_prism() +
      xlab("CSC Concentration (μg/mL)") +
      theme(
        plot.title = element_text(size=11),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 12)#,
        #legend.position = "none"
      ) +
      # scale_fill_grey(start = 1, end = 0.25, name = "Concentration (μg/mL)") +
      guides(fill = guide_legend(override.aes = list(shape = NA))) +
      facet_grid(~State) 
    
    assign(paste("High", y, i, sep = "_"), highplot)
  }
}

Gas_1HrP <- `High_Gastroc_1-Hour P` +
  stat_compare_means(method = "wilcox.test", ref.group = "0", label.y = 90, label = "p.signif", paired = F) +
  stat_compare_means(label.y = 95, paired = TRUE, label = "p.signif") + 
  ggtitle("Permeabilized First")

Gas_1Hr <- `High_Gastroc_1-Hour` +
  stat_compare_means(method = "wilcox.test", ref.group = "0", label.y = 90, label = "p.signif", paired = F) +
  stat_compare_means(label.y = 95, paired = TRUE, label = "p.signif") + 
  ggtitle("Incubated First")


gastroc_higher_concentrations <- cowplot::plot_grid(Gas_1Hr, Gas_1HrP)

Sol_1HrP <- `High_Soleus_1-Hour P` + 
  coord_cartesian(ylim = c(0,200), expand = T) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,250,20)) +
  stat_compare_means(method = "wilcox", ref.group = "0", label.y = 170, label = "p.signif", paired = T)+
  stat_compare_means(label.y = 180, paired = T, label = "p.signif") + 
  ggtitle("Permeabilized First")

Sol_1Hr <- `High_Soleus_1-Hour` + 
  coord_cartesian(ylim = c(0,200), expand = T) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,250,20)) +
  stat_compare_means(method = "wilcox", ref.group = "0", label.y = 170, label = "p.signif", paired = T)+
  stat_compare_means(label.y = 180, paired = T, label = "p.signif")+ 
  ggtitle("Incubated First")

soleus_higher_concentrations <- cowplot::plot_grid(Sol_1Hr, Sol_1HrP)


Aor_1Hr <- `High_Aorta_1-Hour` + 
  coord_cartesian(ylim = c(0,50), expand = T) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,250,5)) +
  stat_compare_means(method = "wilcox.test", ref.group = "0", label.y = 43, label = "p.format", paired = F) +
  stat_compare_means(label.y = 45, paired = F) + 
  ggtitle("Incubated First")


Aor_1HrP <- `High_Aorta_1-Hour P` + 
  coord_cartesian(ylim = c(0,50), expand = T) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,250,5)) +
  #stat_compare_means(method = "wilcox.test", ref.group = "0", label.y = 43, label = "p.format", paired = F) +
  #stat_compare_means(label.y = 45, paired = F) + 
  ggtitle("Permeabilized First")

aorta_higher_concentrations <- cowplot::plot_grid(Aor_1Hr, Aor_1HrP)

high_rcr <- ggplot(Highdata_long, aes(x = Tissue, y = RCR, fill = Concentration)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, color = "black", size = 1) +
  geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,10), breaks = seq(0,10,1)) +
  ylab("Respiratory Control Ratio") +
  theme_prism() +
  xlab("Incubation Time") +
  theme(
    plot.title = element_text(size=11)
  ) +
  #scale_fill_brewer(palette = "Paired") +
  #scale_fill_grey(start = 1, end = 0.25, name = "Concentration (μg/mL)") +
  guides(fill = guide_legend(override.aes = list(shape = NA)))

for(y in unique(Highdata$Tissue)){
  highplot <- ggplot(data = subset(Highdata_long, Tissue %in% paste(y)), aes(x = Concentration, y = Rate, fill = Concentration, colour = Subject)) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
    stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, color = "black", size = 1) +
    geom_point(size = 1.75) +
    scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,10)) +
    ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
    theme_prism() +
    xlab("CSC Concentration (μg/mL)") +
    theme(
      plot.title = element_text(size=11),
      strip.background = element_blank(),
      axis.text.x = element_text(size = 12)#,
      #legend.position = "none"
    ) +
    scale_fill_grey(start = 1, end = 0.25, name = "Concentration (μg/mL)") +
    guides(fill = guide_legend(override.aes = list(shape = NA))) +
    facet_grid(~State) 
  
  assign(paste("High", y, sep = "_"), highplot)
}


High_Soleus <- High_Soleus + 
  coord_cartesian(ylim = c(0,200), expand = T) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,250,20))


whenwillthisend <- readxl::read_xlsx('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Project 2\\Analysis File2 0-5.xlsx', sheet = "Final Data", col_names = TRUE)

whenwillthisend <- whenwillthisend[-c(13:16)]

whenwillthisend <- na.omit(whenwillthisend)

whenwillthisend$RCR <- whenwillthisend$GMDS/whenwillthisend$GM

whenwillthisend$Concentration <- scales::percent(whenwillthisend$Concentration)

whenwillthisend$Concentration <- factor(whenwillthisend$Concentration, levels = unique(whenwillthisend$Concentration))

whenwillthisend <- separate(whenwillthisend, Subject, "Subject", " ")

whenwillthisend$Response <- whenwillthisend$CytC/whenwillthisend$GMDS

whenwillthisend_long <- pivot_longer(whenwillthisend, 8:12, names_to = "State", values_to = "Rate")

whenwillthisend_long$State <- factor(whenwillthisend_long$State, levels = unique(whenwillthisend_long$State))

for(y in unique(whenwillthisend$Tissue)){
  whenwillthisendplot <- ggplot(data = subset(subset(whenwillthisend_long, Tissue %in% paste(y)), State %in% c("GM","GMD", "GMDS")), aes(x = Concentration, y = Rate, fill = Concentration)) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
    stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, size = 1, color = "black") +
    geom_path(aes(group = Subject, color = Subject), size = 1.25) +
    geom_point(position = position_dodge(0.9)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,120), breaks = seq(0,120,10)) +
    ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
    theme_prism() +
    theme(
      plot.title = element_text(size=11),
      strip.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.placement = "outside",
      axis.title.x = element_blank()
      #legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    guides(fill = guide_legend(override.aes = list(shape = NA))) +
    new_scale_color() +
    geom_point(data = subset(subset(filter(whenwillthisend_long, Response >= 1.1), Tissue %in% paste(y)), State %in% c("GM","GMD", "GMDS")), position = position_dodge(0.9), color = "yellow") +
    facet_grid(~State, switch = 'x') +
    labs(caption = "Yellow indicates CytC Response > 1.1")
  
  assign(paste("0-5", y, sep = "_"), whenwillthisendplot)
}

`0-5_Gastroc` <- `0-5_Gastroc` +
  stat_compare_means(paired = T, label.y = 115, label.x = 1.75, method = 'anova') +
  stat_compare_means(method = "t.test", paired = F, label.y = 110, ref.group = '0.0%', label = "p.format",                      hide.ns = F, size = 2)


incubation_es <- ggplot(data = subset(subset(whenwillthisend_long, Tissue %in% "Gastroc"), State %in% c("GMDS")), aes(x = Concentration, y = Rate, fill = Concentration)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, size = 1, color = "black") +
  #geom_path(aes(group = Subject, color = Subject), size = 1.25) +
  #geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,90), breaks = seq(0,120,10)) +
  ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  theme(
    #plot.title = element_text(size=11),
    strip.background = element_blank(),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ggtitle("Inhibition of Maximal Mitochondrial Respiration by CSC in Gastrocnemius")  +
  geom_bracket(xmin = 1, xmax = 5, y.position = 80, label = "Inhibition = ~25%; Cohen's d = 1.56", size = 1.25, fontface = 2, inherit.aes = F, label.size = 5)

#ggsave('incubation_es.png', incubation_es, width = 10)


## Gastroc Incubation Changes

incubation_percent <- pivot_wider(subset(whenwillthisend_long, State %in% "GMDS"), id_cols = c("Subject", "Tissue"), names_from = Concentration, values_from = Rate)

incubation_percent$`1%` <- -(1-(incubation_percent$`1%`/incubation_percent$`0%`)) * 100
incubation_percent$`2%` <- -(1-(incubation_percent$`2%`/incubation_percent$`0%`)) * 100
incubation_percent$`3%` <- -(1-(incubation_percent$`3%`/incubation_percent$`0%`)) * 100
incubation_percent$`4%` <- -(1-(incubation_percent$`4%`/incubation_percent$`0%`)) * 100
incubation_percent$`5%` <- -(1-(incubation_percent$`5%`/incubation_percent$`0%`)) * 100
incubation_percent$`0%` <- -(1-(incubation_percent$`0%`/incubation_percent$`0%`)) * 100

incubation_change <- pivot_longer(data = incubation_percent, 3:8, names_to = "Change", values_to = "Percent Change")



gastroc_incubation_percent_change <- ggplot(subset(incubation_change, Tissue %in% "Gastroc"), aes(x = Change, y = `Percent Change`, group = 1)) +
  geom_line(aes(group = Subject, color = Subject), alpha = 0.5) +
  geom_point(aes(group = Subject, color = Subject), alpha = 0.5) +
  geom_path(stat = "summary", fun.y = "mean", size = 1) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), 
               width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "point", fun.y = "mean", position = position_dodge(0.9), 
               na.rm = T, size = 2, color = "black") +
  theme_prism() +
  scale_y_continuous(expand = c(0,0), limits = c(-100, 25), breaks = seq(-100, 25, 10)) +
  xlab("Concentration of CSC") +
  ggtitle('Percent Change with SEM') +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = 0, size = 1, linetype = "longdash") +
  stat_summary(aes(label = round(..y.., 2)), fun.y=mean, geom = "text", vjust = 10)

soleus_incubation_es <- ggplot(data = subset(subset(whenwillthisend_long, Tissue %in% "Soleus"), State %in% c("GMDS")), aes(x = Concentration, y = Rate, fill = Concentration)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9), na.rm = T, size = 1, color = "black") +
  #geom_path(aes(group = Subject, color = Subject), size = 1.25) +
  #geom_point(position = position_dodge(0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,120,10)) +
  ylab(expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  theme(
    #plot.title = element_text(size=11),
    strip.background = element_blank(),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ggtitle("Inhibition of Maximal Mitochondrial Respiration by CSC in Soleus")  +
  geom_bracket(xmin = 1, xmax = 4, y.position = 95, label = "Inhibition = ~10%", size = 1.25, fontface = 2, inherit.aes = F, label.size = 5)

gastroc_95CI <- ggplot(subset(incubation_change, Tissue %in% "Gastroc"), aes(x = Change, y = `Percent Change`, group = 1)) +
  geom_line(aes(group = Subject, color = Subject), alpha = 0.5) +
  geom_point(aes(group = Subject, color = Subject), alpha = 0.5) +
  geom_path(stat = "summary", fun.y = "mean", size = 1) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_normal, position = position_dodge(0.9), 
               width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "point", fun.y = "mean", position = position_dodge(0.9), 
               na.rm = T, size = 2, color = "black") +
  theme_prism() +
  scale_y_continuous(expand = c(0,0), limits = c(-100, 25), breaks = seq(-100, 25, 10)) +
  xlab("Concentration of CSC") +
  ggtitle('Percent Change with 95% CI') +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = 0, size = 1, linetype = "longdash") +
  stat_summary(aes(label = round(..y.., 2)), fun.y=mean, geom = "text", vjust = 10)

`0-5_Soleus` <- `0-5_Soleus` +
  #stat_compare_means(paired = T, label.y = 175, label.x = 1.75, method = 'anova') +
  #stat_compare_means(method = "t.test", paired = F, label.y = 165, ref.group = '0.0%', label = "p.format",                      hide.ns = F, size = 2) +
  scale_y_continuous(expand = c(0,0), limits = c(0,200), breaks = seq(0,200,25)) +
  stat_compare_means(label.y = 180, paired = TRUE) +
  stat_compare_means(label.y = 170, ref.group = "0.0%", paired = TRUE, label = "p.format")

soleus_incubation_percent_change <- ggplot(subset(na.omit(incubation_change), Tissue %in% "Soleus"), aes(x = Change, y = `Percent Change`, group = 1)) +
  geom_line(aes(group = Subject, color = Subject), alpha = 0.5, na.rm = TRUE) +
  geom_point(aes(group = Subject, color = Subject), alpha = 0.5, na.rm = TRUE) +
  geom_path(stat = "summary", fun.y = "mean", size = 1, na.rm = TRUE) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge(0.9), 
               width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "point", fun.y = "mean", position = position_dodge(0.9), 
               na.rm = T, size = 2, color = "black") +
  theme_prism() +
  scale_y_continuous(expand = c(0,0), limits = c(-100, 100), breaks = seq(-100, 100, 25)) +
  xlab("Concentration of CSC") +
  ggtitle('Percent Change with SEM') +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = 0, size = 1, linetype = "longdash") +
  stat_summary(aes(label = round(..y.., 2)), fun.y=mean, geom = "text", vjust = 10)

soleus_95CI <- ggplot(subset(incubation_change, Tissue %in% "Soleus"), aes(x = Change, y = `Percent Change`, group = 1)) +
  geom_line(aes(group = Subject, color = Subject), alpha = 0.5) +
  geom_point(aes(group = Subject, color = Subject), alpha = 0.5) +
  geom_path(stat = "summary", fun.y = "mean", size = 1) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_normal, position = position_dodge(0.9), 
               width = 0.3, size = 1, na.rm = T) +
  stat_summary(geom = "point", fun.y = "mean", position = position_dodge(0.9), 
               na.rm = T, size = 2, color = "black") +
  theme_prism() +
  coord_cartesian(ylim = c(-75, 75)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-75, 75, 25)) +
  ggtitle('Percent Change with 95% CI') +
  xlab("Concentration of CSC") +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = 0, size = 1, linetype = "longdash") +
  stat_summary(aes(label = round(..y.., 2)), fun.y=mean, geom = "text", vjust = 10)
