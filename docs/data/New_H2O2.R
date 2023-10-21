library(readxl)
library(dplyr)
library(tidyr)
library(ggpubr)
library(drc)
library(ggnewscale)
library(EnvStats)
library(ggprism)
library(flextable)
library(reshape2)
library(ggbeeswarm)
library(ARTool)
library(ggpattern)
library(rstatix)
library(data.table)
library(plotrix)
library(dunn.test)

smoke_comparisons <- list(c("Control", "Smoke"))

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

# Import H2O2 data ----
h2o2_data <- read_excel('C:\\Users\\u1159489\\Box\\LayecLab\\CS-THNR\\O2K Analysis Files\\H2O2\\H2O2.xlsx', sheet = "Final H2O2 Data")

h2o2_data$Tissue[h2o2_data$Tissue == "Gastroc"] <- "Gastrocnemius"

h2o2_data[h2o2_data == 0] <- NA

h2o2_data <- na.omit(h2o2_data)

h2o2_data$Group[ h2o2_data$Group == "Con"] <- "Control"

h2o2_data$Group[ h2o2_data$Group == "Smo"] <- "Smoke"

h2o2_data <- separate( h2o2_data, Subject, into = "Subject", sep = " ")

h2o2_data$Smoke <- factor(h2o2_data$Group, levels = unique(h2o2_data$Group))

h2o2_data$Tissue <- factor(h2o2_data$Tissue, levels = unique(h2o2_data$Tissue))

h2o2_data$Sex <- factor(h2o2_data$Sex, levels = unique(h2o2_data$Sex))

h2o2_data <- h2o2_data |> 
  mutate_at(vars(9:14), list(percent = ~./GMS))


h2o2_long <- pivot_longer(h2o2_data[1:14], cols = 7:14, names_to = "State", values_to = "Rate")

h2o2_long$State <- factor(h2o2_long$State, levels = unique( h2o2_long$State))

h2o2_long_percent <- pivot_longer(h2o2_data[-c(7:14)], cols = 8:13, names_to = "State", values_to = "Rate")

h2o2_long_percent <- separate(h2o2_long_percent, State, into = c("State", "fake"), sep = "_")

h2o2_long_percent$State <- factor(h2o2_long_percent$State, levels = unique(h2o2_long_percent$State))

h2o2_data$RCR <- h2o2_data$`ADP 5000`/h2o2_data$GMS

h2o2_data$Condition <- paste(h2o2_data$Tissue, h2o2_data$Group)


# H2O2 Kinetics -----
h2o2_mm <- pivot_longer(h2o2_data, 9:14, values_to = "Rate", names_to = "Concentration")

h2o2_mm <- separate(h2o2_mm, Concentration, into = c("ADP", "Concentration"), sep = " ")

h2o2_mm$Concentration[is.na(h2o2_mm$Concentration)] <- 0

h2o2_mm$Concentration <- as.numeric( h2o2_mm$Concentration)


## Individual H2O2 ----
individ_h2o2_kinetics <- data.frame()
individ_h2o2_mml <- data.frame()

for(z in unique(h2o2_mm$Subject)){
  for(y in unique(h2o2_mm$Group)){
    for(x in unique(h2o2_mm$Tissue)){
      skip_to_next <- FALSE
      
      # Note that print(b) fails since b doesn't exist
      
      tryCatch({
        h2o2_model <- drm(Rate ~ Concentration, data = subset(subset(subset(h2o2_mm, Subject %in% paste(z)), Smoke %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower Limit", "Upper Limit", "IC50")))
        
        h2o2_model$Tissue <- paste(x)
        
        h2o2_model$Smoke <- paste(y)
        
        h2o2_model$Subject <- paste(z)
        
        coefs <- data.frame("C" = c(h2o2_model$coefficients[1]), "Imax" = c(h2o2_model$coefficients[2]), "IC50" = c(h2o2_model$coefficients[3]))
        
        coefs$Inhibition <- -(1-coefs$Imax/coefs$C) 
        
        coefs$Tissue <- paste(x)
        
        coefs$Smoke <- paste(y)
        
        coefs$Subject <- paste(z)
        
        individ_h2o2_kinetics <- rbind(individ_h2o2_kinetics, coefs)
        
        mml <- data.frame(S = seq(0, max(h2o2_mm$Concentration), length.out = 100))
        
        mml$v <- predict(h2o2_model, newdata =  mml)
        
        mml$Tissue <- paste(x) 
        
        mml$Smoke <- paste(y)
        
        mml$Subject <- paste(z)
        
        individ_h2o2_mml <- rbind(individ_h2o2_mml, mml)
        
        
      }, error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { next }
    }
  }
}

individ_h2o2_kinetics$Condition <- paste(individ_h2o2_kinetics$Tissue, individ_h2o2_kinetics$Smoke)
individ_h2o2_kinetics$Smoke <- factor(individ_h2o2_kinetics$Smoke, levels = unique(individ_h2o2_kinetics$Smoke))
individ_h2o2_kinetics$Smoke <- factor(individ_h2o2_kinetics$Smoke, levels = unique(individ_h2o2_kinetics$Smoke))

individ_h2o2_mml$Condition <- paste(individ_h2o2_mml$Tissue, individ_h2o2_mml$Smoke)

individ_h2o2_kinetics$Smoke <- factor(individ_h2o2_kinetics$Smoke, levels = unique(individ_h2o2_kinetics$Smoke))
individ_h2o2_kinetics$Tissue <- factor(individ_h2o2_kinetics$Tissue, levels = unique(individ_h2o2_kinetics$Tissue))

h2o2_errors <- h2o2_mm[-c(1,3,5,6,7,8,9)] |> 
  group_by(Tissue, Group, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

h2o2_errors$min <- h2o2_errors$Rate_mean - (h2o2_errors$Group == "Smoke")*h2o2_errors$Rate_se
h2o2_errors$max <- h2o2_errors$Rate_mean + (h2o2_errors$Group == "Control")*h2o2_errors$Rate_se

h2o2_long$Condition <- paste(h2o2_long$Tissue, h2o2_long$Group)

## Plots -------------------
identify_outliers(data = subset(individ_h2o2_kinetics, Condition %in% "Soleus Smoke"), Imax)
imax_anova <- anova(art(data = individ_h2o2_kinetics[-c(35, 36),], Imax ~ Smoke * Tissue))
effectsize::eta_squared(imax_anova)

aggregate(Imax ~ Smoke * Tissue, data = individ_h2o2_kinetics[-c(35, 36),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

h2o2_imax <- ggplot(data = individ_h2o2_kinetics[-c(35, 36),], aes(x = Tissue, y = Imax, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 0.25), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.25, 0.05)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +  
  geom_beeswarm(data = individ_h2o2_kinetics[-c(35, 36),], aes(x = Tissue, y = Imax, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 0.25*0.99, label = 'Smoke: p = 0.877', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.25*0.94, label = 'Tissue: p = 0.188', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.25*0.89, label = 'Interaction: p = 0.787', hjust = 0, fontface = 2)

identify_outliers(data = subset(individ_h2o2_kinetics, Condition %in% "Gastrocnemius Control"), IC50)
IC50_anova <- anova(art(data = individ_h2o2_kinetics[-c(33, 34, 36, 41, 42),], IC50 ~ Smoke * Tissue))
effectsize::eta_squared(IC50_anova)

aggregate(IC50 ~ Smoke * Tissue, data = individ_h2o2_kinetics[-c(33, 34, 36, 41, 42),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

h2o2_IC50_dunn <- dunn.test(individ_h2o2_kinetics$IC50[-c(33, 34, 36, 41, 42)], individ_h2o2_kinetics$Condition[-c(33, 34, 36, 41, 42)], list = T, method = "none", kw = F)
Sidak(c(h2o2_IC50_dunn$P[2], h2o2_IC50_dunn$P[5]))
Sidak(c(h2o2_IC50_dunn$P[1], h2o2_IC50_dunn$P[6]))
rstatix::cohens_d(data = subset(individ_h2o2_kinetics[-c(33, 34, 36, 41, 42),], Tissue %in% "Gastrocnemius"), formula = IC50 ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_kinetics[-c(33, 34, 36, 41, 42),], Tissue %in% "Soleus"), formula = IC50 ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_kinetics[-c(33, 34, 36, 41, 42),], Smoke %in% "Control"), formula = IC50 ~ Tissue)
rstatix::cohens_d(data = subset(individ_h2o2_kinetics[-c(33, 34, 36, 41, 42),], Smoke %in% "Smoke"), formula = IC50 ~ Tissue)

h2o2_km <- ggplot(data = individ_h2o2_kinetics[-c(33, 34, 36, 41, 42),], aes(x = Tissue, y = IC50, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 120), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +  
  geom_beeswarm(data = individ_h2o2_kinetics[-c(33, 34, 36, 41, 42),], aes(x = Tissue, y = IC50, color = Condition), dodge.width = .95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 120*0.99, label = 'Smoke: p = 0.037', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 120*0.94, label = 'Tissue: p = 0.002', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 120*0.89, label = 'Interaction: p = 0.223', hjust = 0, fontface = 2) +
  geom_bracket(xmin = 0.755, xmax = 1.245, y.position = 40, inherit.aes = F, label = "p = 0.043", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = 1.755, xmax = 2.245, y.position = 117, inherit.aes = F, label = "p = 0.052", size = 1, fontface = 2, tip.length = 0.01)


identify_outliers(data = subset(individ_h2o2_kinetics, Condition %in% "Soleus Control"), C)
C_anova <- anova(art(data = individ_h2o2_kinetics, C ~ Smoke * Tissue))
effectsize::eta_squared(C_anova)

h2o2_vmax <- ggplot(data = individ_h2o2_kinetics, aes(x = Tissue, y = C, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.2)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +  
  geom_beeswarm(data = individ_h2o2_kinetics, aes(x = Tissue, y = C, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 1*0.99, label = 'Smoke: p = 0.468', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 1*0.94, label = 'Tissue: p = 0.249', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 1*0.89, label = 'Interaction: p = 0.427', hjust = 0, fontface = 2)


identify_outliers(data = subset(individ_h2o2_kinetics, Condition %in% "Soleus Control"), Inhibition)
Inhibition_anova <- anova(art(data = individ_h2o2_kinetics[-c(42),], Inhibition ~ Smoke * Tissue))
effectsize::eta_squared(Inhibition_anova)

h2o2_percent <- ggplot(data = individ_h2o2_kinetics[-c(42),], aes(x = Tissue, y = Inhibition, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(-1, 0), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1, 0, 0.1), labels = scales::percent_format()) +
  labs(y=expression(bold(Inhibition~('%'~V[max])))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  scale_x_discrete(position = "top") + 
  geom_beeswarm(data = individ_h2o2_kinetics[-c(42),], aes(x = Tissue, y = Inhibition, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = -0.89, label = 'Smoke: p = 0.169', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = -0.94, label = 'Tissue: p = 0.379', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = -0.99, label = 'Interaction: p = 0.652', hjust = 0, fontface = 2)

for(y in unique(h2o2_mm$Concentration)){
  print(identify_outliers(subset(h2o2_mm, Concentration %in% y), Rate))
}

h2o2_errors$Condition <- paste(h2o2_errors$Tissue, h2o2_errors$Group)

h2o2_kinetics <- ggplot(data = h2o2_mm[-c(212:216),], aes(x = Concentration, y = Rate, color = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = h2o2_errors, aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
  coord_cartesian(ylim = c(0, .8), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.8, 0.1)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  geom_line(data = individ_h2o2_mml, aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  geom_point(data = h2o2_mm[-c(212:216),], aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)

h2o2_kinetics_Gast <- ggplot(data = subset(h2o2_mm[-c(212:216),], Tissue %in% "Gastrocnemius"), aes(x = Concentration, y = Rate, fill = Smoke)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = subset(h2o2_errors, Tissue %in% "Gastrocnemius"), aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Smoke), inherit.aes = F, size = 1, width = 75, color = "black", show.legend = F) +
  coord_cartesian(ylim = c(0, .6), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.6, 0.05)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  geom_line(data = subset(individ_h2o2_mml, Tissue %in% "Gastrocnemius"), aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_fill_manual(breaks = c("Control", "Smoke"), values = c("white", "grey33")) +
  ggtitle("Gastrocnemius") +
  geom_point(data = subset(h2o2_mm[-c(212:216),], Tissue %in% "Gastrocnemius"), aes(x = Concentration, y = Rate, fill = Smoke), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, color = "black", show.legend = F)

h2o2_kinetics_Soleus <- ggplot(data = subset(h2o2_mm[-c(212:216),], Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, fill = Smoke)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = subset(h2o2_errors, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Smoke), inherit.aes = F, size = 1, width = 75, color = "black", show.legend = F) +
  coord_cartesian(ylim = c(0, .6), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.6, 0.05)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  geom_line(data = subset(individ_h2o2_mml, Tissue %in% "Soleus"), aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_fill_manual(breaks = c("Control", "Smoke"), values = c("white", "grey33")) +
  ggtitle("Soleus") +
  geom_point(data = subset(h2o2_mm[-c(212:216),], Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, fill = Smoke), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, color = "black", show.legend = F)

## Stats--------------
### H2O2 --------------
h2o2_long$Condition <- factor(h2o2_long$Condition, levels = unique(h2o2_long$Condition))

outliers <- data.frame()

for(y in unique(h2o2_long$Condition)){
  for(x in names(h2o2_data[c(7:14)])){
    skip_to_next <- FALSE
    
    tryCatch({
      out <- identify_outliers(data = subset(subset(h2o2_long, Condition %in% paste(y)), State %in% paste(x)), Rate)
      
      # resp_outliers <- outliers[which(isTRUE(outliers$is.extreme))]
      #
      # assign(paste("outliers", x, y, sep = "_"), resp_outliers)
      
      if(isTRUE(outliers$is.extreme)){
        print(paste(outliers[which(outliers$is.extreme == TRUE)], y, x, sep = "_"))
      }
      
      b <- data.frame(c(out[which(out$is.extreme == TRUE)]))
      
      a <- data.frame(b[1])
      
      a$State <- paste(x)
      
      a$Condition <- paste(y)
      
      outliers <- rbind(outliers, a)
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  }
}

h2o2_long$Group <- factor(h2o2_long$Group, levels = unique(h2o2_long$Group))

for(y in unique(h2o2_long$State)){
  print(subset(h2o2_long[-c(273, 279, 280, 281, 285, 287, 288),], State %in% paste(y)) |>
          group_by(Tissue) |>
          levene_test(Rate ~ Group))
  
  print(as.data.table(subset(h2o2_long[-c(273, 279, 280, 281, 285, 287, 288),], State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
                                                                                                         P.value = shapiro.test(Rate)$p.value),
                                                                                                      by = .(Tissue, Group)])
  
}

anova_states_results <- data.frame()
anova_eff_results <- data.frame()

for(y in unique(h2o2_long$State)){
  anova_states <- anova(lm(Rate ~ Group * Tissue, data = subset(na.omit(h2o2_long[-c(273, 279, 280, 281, 285, 287, 288),]), State %in% paste(y))))
  anova_states$State <- paste(y)
  
  anova_eff <- effectsize::eta_squared(anova_states)
  
  print(anova_states)
  print(effectsize::eta_squared(anova_states))
  
  anova_states_results <- rbind(anova_states_results, anova_states[c(1:3),])
  anova_eff_results <- rbind(anova_eff_results, anova_eff)
  
}

cbind(anova_states_results, anova_eff_results)


mainplot_h2o2 <- ggplot(data = subset(h2o2_long[-c(273, 279, 280, 281, 287, 288),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.1)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  geom_beeswarm(data = subset(h2o2_long[-c(273, 279, 280, 281, 287, 288),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  scale_fill_manual(values = c("white", "white")) +
  #ggtitle(expression(bold(Gastrocnemius~H['2']*O['2']*~Emission))) +
  #GM
  annotate('text', x = .51, y = 1*0.99, label = 'Smoke: p = 0.542', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 1*0.94, label = 'Tissue: p = 0.028', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 1*0.89, label = 'Interaction: p = 0.281', hjust = 0, fontface = 2) +
  #GMS
  annotate('text', x = 1.51, y = 1*0.99, label = 'Smoke: p = 0.198', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 1*0.94, label = 'Tissue: p = 0.191', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 1*0.89, label = 'Interaction: p = 0.252', hjust = 0, fontface = 2) +
  #ADP25
  annotate('text', x = 2.51, y = 1*0.99, label = 'Smoke: p = 0.168', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 1*0.94, label = 'Tissue: p = 0.046', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 1*0.89, label = 'Interaction: p = 0.196', hjust = 0, fontface = 2) +
  #ADP50
  annotate('text', x = 3.51, y = 1*0.99, label = 'Smoke: p = 0.761', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 1*0.94, label = 'Tissue: p = 0.233', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 1*0.89, label = 'Interaction: p = 0.355', hjust = 0, fontface = 2) +
  #ADP100
  annotate('text', x = 4.51, y = 1*0.99, label = 'Smoke: p = 0.423', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 1*0.94, label = 'Tissue: p = 0.080', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 1*0.89, label = 'Interaction: p = 0.153', hjust = 0, fontface = 2) +
  #ADP250
  annotate('text', x = 5.51, y = 1*0.99, label = 'Smoke: p = 0.643', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 1*0.94, label = 'Tissue: p = 0.125', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 1*0.89, label = 'Interaction: p = 0.538', hjust = 0, fontface = 2) +
  #ADP5000
  annotate('text', x = 6.51, y = 1*0.99, label = 'Smoke: p = 0.912', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 1*0.94, label = 'Tissue: p = 0.052', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 1*0.89, label = 'Interaction: p = 0.946', hjust = 0, fontface = 2)


# Import Resp data ----
resp_data <- read_excel('C:\\Users\\u1159489\\Box\\LayecLab\\CS-THNR\\O2K Analysis Files\\H2O2\\H2O2.xlsx', sheet = "Final Resp Data")

resp_data$Tissue[resp_data$Tissue == "Gastroc"] <- "Gastrocnemius"

resp_data <- na.omit(resp_data)

resp_data$Group[ resp_data$Group == "Con"] <- "Control"

resp_data$Group[ resp_data$Group == "Smo"] <- "Smoke"

resp_data <- separate( resp_data, Subject, into = "Subject", sep = " ")

resp_data$Smoke <- factor(resp_data$Group, levels = unique(resp_data$Group))

resp_data$Tissue <- factor(resp_data$Tissue, levels = unique(resp_data$Tissue))

resp_data$Sex <- factor(resp_data$Sex, levels = unique(resp_data$Sex))

resp_long <- pivot_longer( resp_data, cols = 7:14, names_to = "State", values_to = "Rate")

resp_long$State <- factor( resp_long$State, levels = unique( resp_long$State))

resp_data$RCR <- resp_data$`ADP 5000`/resp_data$GMS

resp_data$Condition <- paste(resp_data$Tissue, resp_data$Group)


# Resp kinetics ----

adp_mm <- pivot_longer(resp_data, 9:14, values_to = "Rate", names_to = "Concentration")

adp_mm <- separate(adp_mm, Concentration, into = c("ADP", "Concentration"), sep = " ")

adp_mm$Concentration[is.na(adp_mm$Concentration)] <- 0

adp_mm$Concentration <- as.numeric(adp_mm$Concentration)


## Individual ADP ----
individ_adp_kinetics <- data.frame()
individ_adp_mml <- data.frame()

for(z in unique(adp_mm$Subject)){
  for(y in unique(adp_mm$Smoke)){
    for(x in unique(adp_mm$Tissue)){
      skip_to_next <- FALSE
      
      # Note that print(b) fails since b doesn't exist
      
      tryCatch({
        adp_model <- drm(Rate ~ Concentration, data = subset(subset(subset(adp_mm, Subject %in% paste(z)), Smoke %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower","Vmax", "Km")))
        
        adp_model$Tissue <- paste(x)
        
        adp_model$Smoke <- paste(y)
        
        adp_model$Subject <- paste(z)
        
        coefs <- data.frame("Vmax" = c(adp_model$coefficients[2]), "Km" = c(adp_model$coefficients[3]))
        
        coefs$Tissue <- paste(x)
        
        coefs$Smoke <- paste(y)
        
        coefs$Subject <- paste(z)
        
        individ_adp_kinetics <- rbind(individ_adp_kinetics, coefs)
        
        mml <- data.frame(S = seq(0, max(adp_mm$Concentration), length.out = 100))
        
        mml$v <- predict(adp_model, newdata =  mml)
        
        mml$Tissue <- paste(x) 
        
        mml$Smoke <- paste(y)
        
        mml$Subject <- paste(z)
        
        individ_adp_mml <- rbind(individ_adp_mml, mml)
        
        
      }, error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { next }
    }
  }
}

individ_adp_kinetics$Condition <- paste(individ_adp_kinetics$Tissue, individ_adp_kinetics$Smoke)
individ_adp_kinetics$Smoke <- factor(individ_adp_kinetics$Smoke, levels = unique(individ_adp_kinetics$Smoke))
individ_adp_kinetics$Tissue <- factor(individ_adp_kinetics$Tissue, levels = unique(individ_adp_kinetics$Tissue))


individ_adp_mml$Condition <- paste(individ_adp_mml$Tissue, individ_adp_mml$Smoke)
adp_mm$Condition <- paste(adp_mm$Tissue, adp_mm$Smoke)


adp_errors <- adp_mm[-c(1,3,5,6,7,8,9)] |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

adp_errors$min <- adp_errors$Rate_mean - (adp_errors$Smoke == "Smoke")*adp_errors$Rate_se
adp_errors$max <- adp_errors$Rate_mean + (adp_errors$Smoke == "Control")*adp_errors$Rate_se

resp_kinetics <- ggplot(data = adp_mm, aes(x = Concentration, y = Rate, fill = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = adp_errors, aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, color = "black", show.legend = F) +
  coord_cartesian(ylim = c(0, 100), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,10)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  geom_line(data = individ_adp_mml, aes(x = S, y = v, linetype = Condition), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c(1,3,1,3)) +
  scale_fill_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("white", "grey50", "grey75", "grey25")) +
  geom_point(data = adp_mm, aes(x = Concentration, y = Rate, fill = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, color = "black", show.legend = F)

resp_long$Condition <- paste(resp_long$Tissue, resp_long$Smoke)




identify_outliers(data = subset(individ_adp_kinetics, Condition %in% "Soleus Control"), Km)
km_anova <- anova(art(data = individ_adp_kinetics, Km ~ Smoke * Tissue))
effectsize::eta_squared(km_anova)

resp_km <- ggplot(data = individ_adp_kinetics[-c(14, 17, 24),], aes(x = Tissue, y = Km, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size  = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 450, 50)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  theme(axis.title.x = element_blank()) +
  geom_beeswarm(data = individ_adp_kinetics[-c(14, 17, 24),], aes(x = Tissue, y = Km, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  scale_fill_manual(values = c("white", "white")) +
  theme(axis.title.x = element_blank())

identify_outliers(data = subset(individ_adp_kinetics, Condition %in% "Gastrocnemius Control"), Vmax)
vmax_anova <- anova(art(data = individ_adp_kinetics[-c(37, 38),], Vmax ~ Smoke * Tissue))
effectsize::eta_squared(vmax_anova)

resp_vmax <- ggplot(data = individ_adp_kinetics[-c(37, 38),], aes(x = Tissue, y = Vmax, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 140), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 140, 10)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  geom_beeswarm(data = individ_adp_kinetics[-c(37, 38),], aes(x = Tissue, y = Vmax, fill = Smoke), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2, spacing = 2) +
  scale_fill_manual(values = c("white", "white")) +
  theme(axis.title.x = element_blank()) +
  geom_vline(xintercept = 1.5, linetype = "longdash", alpha = 0.5, size = 1) +
  #stats
  annotate('text', x = .51, y = 140*0.99, label = 'Smoke: p = 0.047', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 140*0.94, label = 'Tissue: p = 0.048', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 140*0.89, label = 'Interaction: p = 0.403', hjust = 0, fontface = 2)


## Respiration --------------
resp_long$Condition <- paste(resp_long$Tissue, resp_long$Group)

resp_long$Condition <- factor(resp_long$Condition, levels = unique(resp_long$Condition))

outliers_resp <- data.frame()

for(y in unique(resp_long$Condition)){
  for(x in names(resp_data[c(9:14)])){
    skip_to_next <- FALSE
    
    tryCatch({
      out <- identify_outliers(data = subset(subset(resp_long, Condition %in% paste(y)), State %in% paste(x)), Rate)
      
      # resp_outliers <- outliers[which(isTRUE(outliers$is.extreme))]
      #
      # assign(paste("outliers", x, y, sep = "_"), resp_outliers)
      
      if(isTRUE(outliers$is.extreme)){
        print(paste(outliers[which(outliers$is.extreme == TRUE)], y, x, sep = "_"))
      }
      
      b <- data.frame(c(out[which(out$is.extreme == TRUE)]))
      
      a <- data.frame(b[1])
      
      a$State <- paste(x)
      
      a$Condition <- paste(y)
      
      outliers_resp <- rbind(outliers_resp, a)
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  }
}

# resp_descriptives <- aggregate(x = Rate ~ Smoke + Tissue + State, data = cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
# subset(resp_descriptives, State %in% "GMD 5000")

resp_long$Group <- factor(resp_long$Group, levels = unique(resp_long$Group))

for(y in c("GM", "GMS", "ADP 25", "ADP 50", "ADP 100", "ADP 250", "ADP 5000")){
  print(subset(resp_long[-c(315, 318),], State %in% paste(y)) |>
          group_by(Tissue) |>
          levene_test(Rate ~ Condition))
  
  print(as.data.table(subset(resp_long[-c(315, 318),], State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
                                                                                P.value = shapiro.test(Rate)$p.value),
                                                                             by = .(Tissue, Group)])
  
}

anova_resp_results <- data.frame()
anova_resp_eff_results <- data.frame()

for(y in c("GM", "GMS", "ADP 25", "ADP 50", "ADP 100", "ADP 250", "ADP 5000")){
  anova_states <- anova(lm(Rate ~ Group * Tissue, data = subset(resp_long[-c(315, 318),], State %in% paste(y))))
  anova_states$State <- paste(y)
  
  anova_eff <- effectsize::eta_squared(anova_states)
  
  print(anova_states)
  print(effectsize::eta_squared(anova_states))
  
  anova_resp_results <- rbind(anova_resp_results, anova_states[c(1:3),])
  anova_resp_eff_results <- rbind(anova_resp_eff_results, anova_eff)
  
}

cbind(anova_resp_results, anova_resp_eff_results)

for(y in unique(resp_long$State)){
  dunn <- subset(resp_long[-c(315, 318),], State %in% paste(y))
  
  assign(paste("resp_dunn", paste(y), sep = "_"), dunn)
  
}

resp_dunn_GM_result <- dunn.test(`resp_dunn_GM`$Rate, `resp_dunn_GM`$Condition, method = "none", list = T)
Sidak(c(resp_dunn_GM_result$P[2], resp_dunn_GM_result$P[5]))
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Smoke %in% "Control"), State %in% "GM"), formula = Rate ~ Tissue)
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Smoke %in% "Smoke"), State %in% "GM"), formula = Rate ~ Tissue)

resp_dunn_GMS_result <- dunn.test(`resp_dunn_GMS`$Rate, `resp_dunn_GMS`$Condition, method = "none", list = T)
Sidak(c(resp_dunn_GMS_result$P[2], resp_dunn_GMS_result$P[5]))
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Smoke %in% "Control"), State %in% "GMS"), formula = Rate ~ Tissue)
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Smoke %in% "Smoke"), State %in% "GMS"), formula = Rate ~ Tissue)

resp_dunn_d50_result <- dunn.test(`resp_dunn_ADP 50`$Rate, `resp_dunn_ADP 50`$Condition, method = "none", list = T)
Sidak(c(resp_dunn_d50_result$P[2], resp_dunn_d50_result$P[5]))
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Smoke %in% "Control"), State %in% "ADP 50"), formula = Rate ~ Tissue)
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Smoke %in% "Smoke"), State %in% "ADP 50"), formula = Rate ~ Tissue)

resp_dunn_d250_result <- dunn.test(`resp_dunn_ADP 250`$Rate, `resp_dunn_ADP 250`$Condition, method = "none", list = T)
Sidak(c(resp_dunn_d250_result$P[1], resp_dunn_d250_result$P[6]))
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Tissue %in% "Gastrocnemius"), State %in% "ADP 250"), formula = Rate ~ Smoke)
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Tissue %in% "Soleus"), State %in% "ADP 250"), formula = Rate ~ Smoke)

resp_dunn_d5k_result <- dunn.test(`resp_dunn_ADP 5000`$Rate, `resp_dunn_ADP 5000`$Condition, method = "none", list = T)
Sidak(c(resp_dunn_d5k_result$P[1], resp_dunn_d5k_result$P[6]))
Sidak(c(resp_dunn_d5k_result$P[2], resp_dunn_d5k_result$P[5]))
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Tissue %in% "Gastrocnemius"), State %in% "ADP 5000"), formula = Rate ~ Smoke)
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Tissue %in% "Soleus"), State %in% "ADP 5000"), formula = Rate ~ Smoke)
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Smoke %in% "Control"), State %in% "ADP 5000"), formula = Rate ~ Tissue)
rstatix::cohens_d(data = subset(subset(resp_long[-c(315, 318),], Smoke %in% "Smoke"), State %in% "ADP 5000"), formula = Rate ~ Tissue)


resp_kinetics <- ggplot(data = adp_mm, aes(x = Concentration, y = Rate, fill = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = adp_errors, aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, color = "black", show.legend = F) +
  coord_cartesian(ylim = c(0, 100), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,10)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  geom_line(data = individ_adp_mml, aes(x = S, y = v, linetype = Condition), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c(1,3,1,3)) +
  scale_fill_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("white", "grey50", "grey75", "grey25")) +
  geom_point(data = adp_mm, aes(x = Concentration, y = Rate, fill = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, color = "black", show.legend = F)


mainplot_resp <- ggplot(data = subset(resp_long[-c(315, 318),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 140), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 140, 20)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  geom_beeswarm(data = subset(resp_long[-c(315, 318),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  #GM
  annotate('text', x = .51, y = 160*0.99, label = 'Smoke: p = 0.746', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 160*0.94, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 160*0.89, label = 'Interaction: p = 0.301', hjust = 0, fontface = 2) +
  #GMS
  annotate('text', x = 1.51, y = 160*0.99, label = 'Smoke: p = 0.502', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 160*0.94, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 160*0.89, label = 'Interaction: p = 0.363', hjust = 0, fontface = 2) +
  #ADP25
  annotate('text', x = 2.51, y = 160*0.99, label = 'Smoke: p = 0.205', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 160*0.94, label = 'Tissue: p 0.042', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 160*0.89, label = 'Interaction: p = 0.620', hjust = 0, fontface = 2) +
  #ADP50
  annotate('text', x = 3.51, y = 160*0.99, label = 'Smoke: p = 0.211', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 160*0.94, label = 'Tissue: p = 0.159', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 160*0.89, label = 'Interaction: p = 0.856', hjust = 0, fontface = 2) +
  #ADP100
  annotate('text', x = 4.51, y = 160*0.99, label = 'Smoke: p = 0.076', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 160*0.94, label = 'Tissue: p = 0.096', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 160*0.89, label = 'Interaction: p = 0.602', hjust = 0, fontface = 2) +
  #ADP250
  annotate('text', x = 5.51, y = 160*0.99, label = 'Smoke: p = 0.018', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 160*0.94, label = 'Tissue: p = 0.122', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 160*0.89, label = 'Interaction: p = 0.602', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(5.75), xmax = c(6.25), y.position = 125, inherit.aes = F, label = "p = 0.032", size = 1, fontface = 2, tip.length = 0.01) +
  #ADP5000
  annotate('text', x = 6.51, y = 160*0.99, label = 'Smoke: p = 0.007', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 160*0.94, label = 'Tissue: p = 0.049', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 160*0.89, label = 'Interaction: p = 0.905', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(6.75), xmax = c(7.25), y.position = 130, inherit.aes = F, label = "p = 0.048", size = 1, fontface = 2, tip.length = 0.01)

mainplot_resp_Soleus <- ggplot(data = subset(subset(resp_long, Tissue %in% 'Soleus'), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, color = "black") +
  # geom_point(position = position_dodge(0.95), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  # ggpattern::geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", color = "black", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  # scale_fill_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("white", "grey66", "white", "grey66")) +
  scale_fill_manual(values = c("white", "grey33")) +
  coord_cartesian(ylim = c(0, 160), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 160, 20)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  new_scale_fill() +
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(subset(resp_long, Tissue %in% 'Soleus'), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Smoke), 
                dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T, cex = 1.5) +
  scale_fill_manual(values = c("white", "white")) +
  ggtitle("Soleus Respiration") +
  geom_bracket(xmin = c(5.75), xmax = c(6.25), y.position = 125, inherit.aes = F, label = "p = 0.080", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(6.75), xmax = c(7.25), y.position = 125, inherit.aes = F, label = "p = 0.038", size = 1, fontface = 2, tip.length = 0.01)

# H2O2 Per O2

h2o2_per <- data.frame(resp_data[c(1:6, 16, 18)])

h2o2_per <- cbind(h2o2_per, (h2o2_data[8:14]/resp_data[8:14]))

h2o2_per_long <- pivot_longer(h2o2_per, 9:15, names_to = "State", values_to = "Rate")

h2o2_per_long$State <- factor(h2o2_per_long$State, levels = unique(h2o2_per_long$State))

# Individual H2O2 per O2 -----------------------------------------------------------------

h2o2_per_mm <- pivot_longer(h2o2_per, 10:15, values_to = "Rate", names_to = "Concentration")

h2o2_per_mm <- separate(h2o2_per_mm, Concentration, into = c("ADP", "Concentration"), sep = " ")

h2o2_per_mm$Concentration[is.na(h2o2_per_mm$Concentration)] <- 0

h2o2_per_mm$Concentration <- as.numeric(h2o2_per_mm$Concentration)

individ_h2o2_per_kinetics <- data.frame()
individ_h2o2_per_mml <- data.frame()

for(z in unique(h2o2_per_mm$Subject)){
  for(y in unique(h2o2_per_mm$Group)){
    for(x in unique(h2o2_per_mm$Tissue)){
      skip_to_next <- FALSE
      
      # Note that print(b) fails since b doesn't exist
      
      tryCatch({
        h2o2_per_model <- drm(Rate ~ Concentration, data = subset(subset(subset(h2o2_per_mm, Subject %in% paste(z)), Smoke %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower Limit", "Upper Limit", "IC50")))
        
        h2o2_per_model$Tissue <- paste(x)
        
        h2o2_per_model$Smoke <- paste(y)
        
        h2o2_per_model$Subject <- paste(z)
        
        coefs_per <- data.frame("C" = c(h2o2_per_model$coefficients[1]), "Imax" = c(h2o2_per_model$coefficients[2]), "IC50" = c(h2o2_per_model$coefficients[3]))
        
        coefs_per$Inhibition <- -(1-coefs_per$Imax/coefs_per$C) 
        
        coefs_per$Tissue <- paste(x)
        
        coefs_per$Smoke <- paste(y)
        
        coefs_per$Subject <- paste(z)
        
        individ_h2o2_per_kinetics <- rbind(individ_h2o2_per_kinetics, coefs_per)
        
        mml_per <- data.frame(S = seq(0, max(h2o2_per_mm$Concentration), length.out = 100))
        
        mml_per$v <- predict(h2o2_per_model, newdata =  mml_per)
        
        mml_per$Tissue <- paste(x) 
        
        mml_per$Smoke <- paste(y)
        
        mml_per$Subject <- paste(z)
        
        individ_h2o2_per_mml <- rbind(individ_h2o2_per_mml, mml_per)
        
        
      }, error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { next }
    }
  }
}

individ_h2o2_per_kinetics$Smoke <- factor(individ_h2o2_per_kinetics$Smoke, levels = unique(individ_h2o2_per_kinetics$Smoke))
individ_h2o2_per_kinetics$Tissue <- factor(individ_h2o2_per_kinetics$Tissue, levels = unique(individ_h2o2_per_kinetics$Tissue))

individ_h2o2_per_kinetics$Condition <- paste(individ_h2o2_per_kinetics$Tissue, individ_h2o2_per_kinetics$Smoke)

h2o2_per_errors <- h2o2_per_mm[-c(1,3,5,6,7,8,9)] |> 
  group_by(Tissue, Group, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

h2o2_per_errors$min <- h2o2_per_errors$Rate_mean - (h2o2_per_errors$Group == "Smoke")*h2o2_per_errors$Rate_se
h2o2_per_errors$max <- h2o2_per_errors$Rate_mean + (h2o2_per_errors$Group == "Control")*h2o2_per_errors$Rate_se

identify_outliers(data = subset(individ_h2o2_per_kinetics, Condition %in% "Gastrocnemius Control"), Imax)
imax_per_anova <- anova(art(data = individ_h2o2_per_kinetics[-c(28, 29),], Imax ~ Smoke * Tissue))
effectsize::eta_squared(imax_per_anova)

aggregate(Imax ~ Smoke * Tissue, data = individ_h2o2_per_kinetics[-c(28, 29),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

h2o2_per_imax_dunn <- dunn.test(individ_h2o2_per_kinetics$Imax[-c(28, 29)], individ_h2o2_per_kinetics$Condition[-c(28, 29)], list = T, method = "none", kw = F)
Sidak(c(h2o2_per_imax_dunn$P[2], h2o2_per_imax_dunn$P[5]))
Sidak(c(h2o2_per_imax_dunn$P[1], h2o2_per_imax_dunn$P[6]))
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 29)], Tissue %in% "Gastrocnemius"), formula = Imax ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 29)], Tissue %in% "Soleus"), formula = Imax ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 29)], Smoke %in% "Control"), formula = Imax ~ Tissue)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 29)], Smoke %in% "Smoke"), formula = Imax ~ Tissue)

h2o2_per_imax <- ggplot(data = individ_h2o2_per_kinetics[-c(28, 29)], aes(x = Tissue, y = Imax, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, .006), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .006, .001)) +
  labs(y = expression(bold(bolditalic(J)*H['2']*O['2']*'/'*O['2']))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  geom_beeswarm(data = individ_h2o2_per_kinetics[-c(28, 29)], aes(x = Tissue, y = Imax, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 0.006*0.99, label = 'Smoke: p = 0.009', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.006*0.94, label = 'Tissue: p = 0.617', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.006*0.89, label = 'Interaction: p = 0.110', hjust = 0, fontface = 2) +
  geom_bracket(xmin = 0.755, xmax = 1.245, y.position = 0.0046, inherit.aes = F, label = "p = 0.005", size = 1, fontface = 2, tip.length = 0.01)


identify_outliers(data = subset(individ_h2o2_per_kinetics, Condition %in% "Gastrocnemius Control"), IC50)
km_per_anova <- anova(art(data = individ_h2o2_per_kinetics[-c(28, 31),], IC50 ~ Smoke * Tissue))
effectsize::eta_squared(km_per_anova)

aggregate(IC50 ~ Smoke * Tissue, data = individ_h2o2_per_kinetics[-c(28, 31),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

h2o2_per_km_dunn <- dunn.test(individ_h2o2_per_kinetics$IC50[-c(28, 31)], individ_h2o2_per_kinetics$Condition[-c(28, 31)], list = TRUE, method = "hs", kw = F)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 31),], Tissue %in% "Gastrocnemius"), formula = IC50 ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 31),], Tissue %in% "Soleus"), formula = IC50 ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 31),], Smoke %in% "Control"), formula = IC50 ~ Tissue)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 31),], Smoke %in% "Smoke"), formula = IC50 ~ Tissue)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 31),], Condition %in% c("Gastrocnemius Smoke", "Soleus Control")), formula = IC50 ~ Condition)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(28, 31),], Condition %in% c("Gastrocnemius Control", "Soleus Smoke")), formula = IC50 ~ Condition)

h2o2_per_km <- ggplot(data = individ_h2o2_per_kinetics[-c(28, 31),], aes(x = Tissue, y = IC50, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 80), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 80, 20)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  #ggtitle(expression(bold(IC['50']))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_h2o2_per_kinetics[-c(28, 31),], aes(x = Tissue, y = IC50, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 80*0.99, label = 'Smoke: p = 0.655', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 80*0.94, label = 'Tissue: p = 0.007', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 80*0.89, label = 'Interaction: p = 0.033', hjust = 0, fontface = 2)

identify_outliers(data = subset(individ_h2o2_per_kinetics, Condition %in% "Soleus Control"), C)
vmax_per_anova <- anova(art(data = individ_h2o2_per_kinetics[-c(4, 29),], C ~ Smoke * Tissue))
effectsize::eta_squared(vmax_per_anova)

h2o2_per_vmax <- ggplot(data = individ_h2o2_per_kinetics[-c(4, 29),], aes(x = Tissue, y = C, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 0.025), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.025, .005)) +
  labs(y = expression(bold(bolditalic(J)*H['2']*O['2']*'/'*O['2']))) +
  #ggtitle(expression(bold(Initial~Rate))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  geom_beeswarm(data = individ_h2o2_per_kinetics[-c(4, 29),], aes(x = Tissue, y = C, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 0.025*0.99, label = 'Smoke: p = 0.332', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.025*0.94, label = 'Tissue: p = 0.112', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.025*0.89, label = 'Interaction: p = 0.947', hjust = 0, fontface = 2)

identify_outliers(data = subset(individ_h2o2_per_kinetics, Condition %in% "Soleus Control"), Inhibition)
inhib_per_anova <- anova(art(data = individ_h2o2_per_kinetics, Inhibition ~ Smoke * Tissue))
effectsize::eta_squared(inhib_per_anova)

h2o2_per_inhib_dunn <- dunn.test(individ_h2o2_per_kinetics$Inhibition, individ_h2o2_per_kinetics$Condition, list = T, method = "none", kw = F)
Sidak(c(h2o2_per_imax_dunn$P[2], h2o2_per_imax_dunn$P[5]))
Sidak(c(h2o2_per_imax_dunn$P[1], h2o2_per_imax_dunn$P[6]))
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics, Tissue %in% "Gastrocnemius"), formula = Inhibition ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics, Tissue %in% "Soleus"), formula = Inhibition ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics, Smoke %in% "Control"), formula = Inhibition ~ Tissue)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics, Smoke %in% "Smoke"), formula = Inhibition ~ Tissue)

h2o2_per_inhibition <- ggplot(data = individ_h2o2_per_kinetics, aes(x = Tissue, y = Inhibition, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(-1.1, 0), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1.1, 0, 0.1), labels = scales::percent_format()) +
  labs(y=expression(bold(Inhibition~('%'~V[max])))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_discrete(position = "top") + 
  geom_beeswarm(data = individ_h2o2_per_kinetics, aes(x = Tissue, y = Inhibition, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = -1.2*0.9, label = 'Smoke: p = 0.015', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = -1.2*0.95, label = 'Tissue: p = 0.085', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = -1.2*1, label = 'Interaction: p = 0.667', hjust = 0, fontface = 2) +
  geom_bracket(xmin = 0.755, xmax = 1.245, y.position = -.98, inherit.aes = F, label = "p = 0.017", size = 1, fontface = 2, tip.length = -0.02, vjust = 1.75)

h2o2_per_errors$Condition <- paste(h2o2_per_errors$Tissue, h2o2_per_errors$Group)
individ_h2o2_per_mml$Condition <- paste(individ_h2o2_per_mml$Tissue, individ_h2o2_per_mml$Smoke)

plot_h2o2_per_kinetics <- ggplot(data = h2o2_per_mm, aes(x = Concentration, y = Rate, color = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = h2o2_per_errors, aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, show.legend = F) +
  coord_cartesian(ylim = c(0, .015), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.015, 0.003)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']*'/'*O['2']))) +
  geom_line(data = individ_h2o2_per_mml, aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  geom_point(data = h2o2_per_mm, aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)

# Stats ------------------------------------

## H2O2 per O2 --------------
h2o2_per_long$Condition <- factor(h2o2_per_long$Condition, levels = unique(h2o2_per_long$Condition))

outliers_per <- data.frame()

for(y in unique(h2o2_per_long$Condition)){
  for(x in names(h2o2_data[c(7:14)])){
    skip_to_next <- FALSE
    
    tryCatch({
      out <- identify_outliers(data = subset(subset(h2o2_per_long, Condition %in% paste(y)), State %in% paste(x)), Rate)
      
      # resp_outliers <- outliers[which(isTRUE(outliers$is.extreme))]
      #
      # assign(paste("outliers", x, y, sep = "_"), resp_outliers)
      
      if(isTRUE(outliers$is.extreme)){
        print(paste(outliers[which(outliers$is.extreme == TRUE)], y, x, sep = "_"))
      }
      
      b <- data.frame(c(out[which(out$is.extreme == TRUE)]))
      
      a <- data.frame(b[1])
      
      a$State <- paste(x)
      
      a$Condition <- paste(y)
      
      outliers_per <- rbind(outliers_per, a)
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  }
}

# resp_descriptives <- aggregate(x = Rate ~ Smoke + Tissue + State, data = cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
# subset(resp_descriptives, State %in% "GMD 5000")

h2o2_per_long$Group <- factor(h2o2_per_long$Group, levels = unique(h2o2_per_long$Group))

for(y in unique(h2o2_per_long$State)){
  print(subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], State %in% paste(y)) |>
          group_by(Tissue) |>
          levene_test(Rate ~ Group))
  
  print(as.data.table(subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
                                                                                                     P.value = shapiro.test(Rate)$p.value),
                                                                                                  by = .(Tissue, Group)])
  
}

anova_per_results <- data.frame()
anova_per_eff_results <- data.frame()

for(y in unique(h2o2_per_long$State)){
  anova_states <- anova(lm(Rate ~ Group * Tissue, data = subset(na.omit(h2o2_per_long[-c(29, 30, 246:249, 251, 252),]), State %in% paste(y))))
  anova_states$State <- paste(y)
  
  anova_eff <- effectsize::eta_squared(anova_states)
  
  print(anova_states)
  print(effectsize::eta_squared(anova_states))
  
  anova_per_results <- rbind(anova_per_results, anova_states[c(1:3),])
  anova_per_eff_results <- rbind(anova_per_eff_results, anova_eff)
  
}

cbind(anova_per_results, anova_per_eff_results)

for(y in unique(h2o2_per_long$State)){
  dunn <- subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], State %in% paste(y))
  
  assign(paste("per_dunn", paste(y), sep = "_"), dunn)
  
}

per_dunn_GMS_result <- dunn.test(`per_dunn_GMS`$Rate, `per_dunn_GMS`$Condition, method = "none", list = T)
Sidak(c(per_dunn_GMS_result$P[2], per_dunn_GMS_result$P[5]))
rstatix::cohens_d(data = subset(subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], Smoke %in% "Control"), State %in% "GMS"), formula = Rate ~ Tissue)
rstatix::cohens_d(data = subset(subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], Smoke %in% "Smoke"), State %in% "GMS"), formula = Rate ~ Tissue)


per_dunn_ADP5k <- dunn.test(`per_dunn_ADP 5000`$Rate, `per_dunn_ADP 5000`$Condition, method = "none", list = T)
Sidak(c(per_dunn_ADP5k$P[1], per_dunn_ADP5k$P[6]))
rstatix::cohens_d(data = subset(subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], Tissue %in% "Gastrocnemius"), State %in% "ADP 5000"), formula = Rate ~ Smoke)
rstatix::cohens_d(data = subset(subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], Tissue %in% "Soleus"), State %in% "ADP 5000"), formula = Rate ~ Smoke)

mainplot_per <- ggplot(data = subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  # geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  # ggpattern::geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.95), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 0.04), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.04, 0.005)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']/bolditalic(J)*O['2']))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  #GM
  annotate('text', x = 0.51, y = 0.04*0.99, label = 'Smoke: p = 0.451', hjust = 0, fontface = 2) +
  annotate('text', x = 0.51, y = 0.04*0.94, label = 'Tissue: p = 0.387', hjust = 0, fontface = 2) +
  annotate('text', x = 0.51, y = 0.04*0.89, label = 'Interaction: p = 0.297', hjust = 0, fontface = 2) +
  #GMS
  annotate('text', x = 1.51, y = 0.04*0.99, label = 'Smoke: p = 0.645', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 0.04*0.94, label = 'Tissue: p = 0.017', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 0.04*0.89, label = 'Interaction: p = 0.426', hjust = 0, fontface = 2) +
  #ADP25
  annotate('text', x = 2.51, y = 0.04*0.99, label = 'Smoke: p = 0.459', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 0.04*0.94, label = 'Tissue: p = 0.401', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 0.04*0.89, label = 'Interaction: p = 0.869', hjust = 0, fontface = 2) +
  #ADP50
  annotate('text', x = 3.51, y = 0.04*0.99, label = 'Smoke: p = 0.161', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 0.04*0.94, label = 'Tissue: p = 0.570', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 0.04*0.89, label = 'Interaction: p = 0.653', hjust = 0, fontface = 2) +
  #ADP100
  annotate('text', x = 4.51, y = 0.04*0.99, label = 'Smoke: p = 0.440', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 0.04*0.94, label = 'Tissue: p = 0.393', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 0.04*0.89, label = 'Interaction: p = 0.220', hjust = 0, fontface = 2) +
  #ADP250
  annotate('text', x = 5.51, y = 0.04*0.99, label = 'Smoke: p = 0.054', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 0.04*0.94, label = 'Tissue: p = 0.913', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 0.04*0.89, label = 'Interaction: p = 0.383', hjust = 0, fontface = 2) +
  #ADP5000
  annotate('text', x = 6.51, y = 0.04*0.99, label = 'Smoke: p = 0.025', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 0.04*0.94, label = 'Tissue: p = 0.791', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 0.04*0.89, label = 'Interaction: p = 0.904', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(6.75), xmax = c(7.25), y.position = 0.01, inherit.aes = F, label = "p = 0.057", size = 1, fontface = 2, tip.length = 0.01)


mainplot_per_Soleus <- ggplot(data = subset(subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], Tissue %in% "Soleus"), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, color = "black") +
  # geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  # ggpattern::geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.95), stat = "summary", color = "black", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  # scale_fill_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("white", "grey33")) +
  scale_fill_manual(values = c("white", "grey33")) +
  coord_cartesian(ylim = c(0, 0.04), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.04, 0.005)) +
  labs(y = expression(bold(bolditalic(J)*H['2']*O['2']*'/'*bolditalic(J)*O['2']))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  new_scale_fill() +
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(subset(h2o2_per_long[-c(29, 30, 246:249, 251, 252),], Tissue %in% "Soleus"), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Smoke), 
                dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T, cex = 1.5) +
  scale_fill_manual(values = c("white", "white")) +
  ggtitle(expression(bold(Soleus~H['2']*O['2']~Emmission~per~O['2']~Consumed)))


# Percentages--------------------------

# CAT Kinetics -----

percent_mm <- pivot_longer(h2o2_data[-c(7:14)], cols = 8:13, names_to = "Concentration", values_to = "Rate")

percent_mm <- separate(percent_mm, Concentration, into = c("ADP", "fake"), sep = "_")

percent_mm <- separate(percent_mm[-c(9)], ADP, into = c("ADP", "Concentration"), sep = " ")

percent_mm$Concentration[is.na(percent_mm$Concentration)] <- 0

percent_mm$Concentration <- as.numeric( percent_mm$Concentration)

percent_kinetic_values <- data.frame()


## Individual H2O2 Percent ----
individ_percent_kinetics <- data.frame()
individ_percent_mml <- data.frame()

for(z in unique(percent_mm$Subject)){
  for(y in unique(percent_mm$Group)){
    for(x in unique(percent_mm$Tissue)){
      skip_to_next <- FALSE
      
      # Note that print(b) fails since b doesn't exist
      
      tryCatch({
        percent_model <- drm(Rate ~ Concentration, data = subset(subset(subset(percent_mm, Subject %in% paste(z)), Smoke %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower Limit", "Upper Limit", "IC50")))
        
        percent_model$Tissue <- paste(x)
        
        percent_model$Smoke <- paste(y)
        
        percent_model$Subject <- paste(z)
        
        coefs <- data.frame("C" = c(percent_model$coefficients[1]), "Imax" = c(percent_model$coefficients[2]), "IC50" = c(percent_model$coefficients[3]))
        
        coefs$Inhibition <- -(1-coefs$Imax/coefs$C) 
        
        coefs$Tissue <- paste(x)
        
        coefs$Smoke <- paste(y)
        
        coefs$Subject <- paste(z)
        
        individ_percent_kinetics <- rbind(individ_percent_kinetics, coefs)
        
        mml <- data.frame(S = seq(0, max(percent_mm$Concentration), length.out = 100))
        
        mml$v <- predict(percent_model, newdata =  mml)
        
        mml$Tissue <- paste(x) 
        
        mml$Smoke <- paste(y)
        
        mml$Subject <- paste(z)
        
        individ_percent_mml <- rbind(individ_percent_mml, mml)
        
        
      }, error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { next }
    }
  }
}


individ_percent_kinetics$Condition <- paste(individ_percent_kinetics$Tissue, individ_percent_kinetics$Smoke)

individ_percent_mml$Condition <- paste(individ_percent_mml$Tissue, individ_percent_mml$Smoke)

individ_percent_kinetics$Smoke <- factor(individ_percent_kinetics$Smoke, levels = unique(individ_percent_kinetics$Smoke))
individ_percent_kinetics$Tissue <- factor(individ_percent_kinetics$Tissue, levels = unique(individ_percent_kinetics$Tissue))

identify_outliers(data = subset(individ_percent_kinetics, Condition %in% "Gastrocnemius Control"), Imax)
imax_percent_anova <- anova(art(data = individ_percent_kinetics, IC50 ~ Smoke * Tissue))
effectsize::eta_squared(imax_percent_anova)

percent_imax <- ggplot(data = individ_percent_kinetics, aes(x = Tissue, y = Imax, fill = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), color = "black", size = 1) +
  theme_prism() +
  scale_fill_manual(breaks = c("Control", "Smoke"), values = c("white", "grey33")) +
  geom_vline(xintercept = 1.5, linetype = "longdash", color = "grey50", size = .6) + 
  coord_cartesian(ylim = c(0, .8), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .8, .1), labels = scales::percent_format()) +
  labs(y = expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
  ggtitle(expression(bold(I[max]))) +
  theme(axis.title.x = element_blank()) + 
  new_scale_fill() +
  geom_beeswarm(data = individ_percent_kinetics, aes(x = Tissue, y = Imax, fill = Smoke), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2, spacing = 2) +
  scale_fill_manual(values = c("white", "white")) +
  theme(axis.title.x = element_blank()) +
  #stats
  annotate('text', x = .51, y = 0.8*0.99, label = 'Smoke: p = 0.357', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.8*0.94, label = 'Tissue: p = 0.288', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.8*0.89, label = 'Interaction: p = 0.357', hjust = 0, fontface = 2)

identify_outliers(data = subset(individ_percent_kinetics, Condition %in% "Gastrocnemius Control"), IC50)
km_percent_anova <- anova(art(data = individ_percent_kinetics[-c(41, 34, 36, 42),], IC50 ~ Smoke * Tissue))
effectsize::eta_squared(km_percent_anova)

h2o2_percent_IC50_dunn <- dunn.test(individ_percent_kinetics$IC50[-c(41, 34, 36, 42)], individ_percent_kinetics$Condition[-c(41, 34, 36, 42)], list = T, method = "none", kw = F)
Sidak(c(h2o2_percent_IC50_dunn$P[2], h2o2_percent_IC50_dunn$P[5]))
Sidak(c(h2o2_percent_IC50_dunn$P[1], h2o2_percent_IC50_dunn$P[6]))
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(41, 34, 36, 42),], Tissue %in% "Gastrocnemius"), formula = Inhibition ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(41, 34, 36, 42),], Tissue %in% "Soleus"), formula = Inhibition ~ Smoke)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(41, 34, 36, 42),], Smoke %in% "Control"), formula = Inhibition ~ Tissue)
rstatix::cohens_d(data = subset(individ_h2o2_per_kinetics[-c(41, 34, 36, 42),], Smoke %in% "Smoke"), formula = Inhibition ~ Tissue)

percent_km <- ggplot(data = individ_percent_kinetics[-c(41, 34, 36, 42),], aes(x = Tissue, y = IC50, fill = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), color = "black", size = 1) +
  theme_prism() +
  scale_fill_manual(breaks = c("Control", "Smoke"), values = c("white", "grey33")) +
  geom_vline(xintercept = 1.5, linetype = "longdash", color = "grey50", size = .6) + 
  coord_cartesian(ylim = c(0, 120), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  ggtitle(expression(bold(IC['50']))) +
  theme(axis.title.x = element_blank()) + 
  new_scale_fill() +
  geom_beeswarm(data = individ_percent_kinetics[-c(41, 34, 36, 42),], aes(x = Tissue, y = IC50, fill = Smoke), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2, spacing = 2) +
  scale_fill_manual(values = c("white", "white")) +
  theme(axis.title.x = element_blank()) +
  #stats
  annotate('text', x = .51, y = 120*0.99, label = 'Smoke: p = 0.048', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 120*0.94, label = 'Tissue: p = 0.021', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 120*0.89, label = 'Interaction: p = 0.195', hjust = 0, fontface = 2)


percent_errors <- percent_mm[-c(1,3,5,6,7,8,9,11)] |> 
  group_by(Tissue, Group, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

percent_errors$min_gast <- percent_errors$mean + (percent_errors$Group == "Smoke")*percent_errors$se
percent_errors$max_gast <- percent_errors$mean - (percent_errors$Group == "Control")*percent_errors$se

percent_errors$min_sol <- percent_errors$mean - (percent_errors$Group == "Smoke")*percent_errors$se
percent_errors$max_sol <- percent_errors$mean + (percent_errors$Group == "Control")*percent_errors$se

percent_errors$Condition <- paste(percent_errors$Tissue, percent_errors$Group)
percent_mm$Condition <- paste(percent_mm$Tissue, percent_mm$Smoke)

percent_kinetics <- ggplot(data = percent_mm, aes(x = Concentration, y = Rate, color = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = percent_errors, aes(x = Concentration, y = mean, ymin = min_gast, ymax = max_gast, color = Condition), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
  coord_cartesian(ylim = c(0, 1.2), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1.2, 0.2)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
  geom_line(data = individ_percent_mml, aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  geom_point(data = percent_mm, aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)

percent_kinetics_Gast <- ggplot(data = subset(percent_mm, Tissue %in% "Gastrocnemius"), aes(x = Concentration, y = Rate, fill = Smoke)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = subset(percent_errors, Tissue %in% "Gastrocnemius"), aes(x = Concentration, y = mean, ymin = min_gast, ymax = max_gast, color = Smoke), inherit.aes = F, size = 1, width = 75, color = "black", show.legend = F) +
  coord_cartesian(ylim = c(0, 1.2), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1.2, 0.2)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
  geom_line(data = subset(individ_percent_mml, Tissue %in% "Gastrocnemius"), aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_fill_manual(breaks = c("Control", "Smoke"), values = c("white", "grey33")) +
  ggtitle("Gastrocnemius (% of GMS)") +
  geom_point(data = subset(percent_mm, Tissue %in% "Gastrocnemius"), aes(x = Concentration, y = Rate, fill = Smoke), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, color = "black", show.legend = F)

percent_kinetics_Soleus <- ggplot(data = subset(percent_mm, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, fill = Smoke)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = subset(percent_errors, Tissue %in% "Soleus"), aes(x = Concentration, y = mean, ymin = min_sol, ymax = max_sol, color = Smoke), inherit.aes = F, size = 1, width = 75, color = "black", show.legend = F) +
  coord_cartesian(ylim = c(0, 1.2), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1.2, 0.2)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
  geom_line(data = subset(individ_percent_mml, Tissue %in% "Soleus"), aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_fill_manual(breaks = c("Control", "Smoke"), values = c("white", "grey33")) +
  ggtitle("Soleus (% of GMS)") +
  geom_point(data = subset(percent_mm, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, fill = Smoke), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, color = "black", show.legend = F)


## H2O2 Percent --------------
h2o2_long_percent$Condition <- paste(h2o2_long_percent$Tissue, h2o2_long_percent$Group)

h2o2_long_percent$Condition <- factor(h2o2_long_percent$Condition, levels = unique(h2o2_long_percent$Condition))

outliers_percent <- data.frame()

for(y in unique(h2o2_long_percent$Condition)){
  for(x in names(h2o2_data[c(9:14)])){
    skip_to_next <- FALSE
    
    tryCatch({
      out <- identify_outliers(data = subset(subset(h2o2_long_percent, Condition %in% paste(y)), State %in% paste(x)), Rate)
      print(out)
      # resp_outliers <- outliers[which(isTRUE(outliers$is.extreme))]
      #
      # assign(paste("outliers", x, y, sep = "_"), resp_outliers)
      
      if(isTRUE(outliers$is.extreme)){
        print(paste(outliers[which(outliers$is.extreme == TRUE)], y, x, sep = "_"))
      }
      
      b <- data.frame(c(out[which(out$is.extreme == TRUE)]))
      
      a <- data.frame(b[1])
      
      a$State <- paste(x)
      
      a$Condition <- paste(y)
      
      outliers_percent <- rbind(outliers_percent, a)
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  }
}

# resp_descriptives <- aggregate(x = Rate ~ Smoke + Tissue + State, data = cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
# subset(resp_descriptives, State %in% "GMD 5000")

h2o2_long_percent$Group <- factor(h2o2_long_percent$Group, levels = unique(h2o2_long_percent$Group))

for(y in c("ADP 25", "ADP 50", "ADP 100", "ADP 250", "ADP 5000")){
  print(subset(h2o2_long_percent[-c(29, 147),], State %in% paste(y)) |>
          group_by(Tissue) |>
          levene_test(Rate ~ Condition))
  
  print(as.data.table(subset(h2o2_long_percent[-c(200),], State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
                                                                                   P.value = shapiro.test(Rate)$p.value),
                                                                                by = .(Tissue, Group)])
  
}

anova_percent_results <- data.frame()
anova_percent_eff_results <- data.frame()

for(y in c("ADP 25", "ADP 50", "ADP 100", "ADP 250", "ADP 5000")){
  anova_states <- anova(lm(Rate ~ Group * Tissue, data = subset(h2o2_long_percent[-c(200),], State %in% paste(y))))
  anova_states$State <- paste(y)
  
  anova_eff <- effectsize::eta_squared(anova_states)
  
  print(anova_states)
  print(effectsize::eta_squared(anova_states))
  
  anova_percent_results <- rbind(anova_percent_results, anova_states[c(1:3),])
  anova_percent_eff_results <- rbind(anova_percent_eff_results, anova_eff)
  
}

cbind(anova_percent_results, anova_percent_eff_results)

for(y in unique(h2o2_long_percent$State)){
  dunn <- subset(h2o2_long_percent[-c(200),], State %in% paste(y))
  
  assign(paste("percent_dunn", paste(y), sep = "_"), dunn)
  
}


percent_dunn_ADP25 <- dunn.test(`percent_dunn_ADP 25`$Rate, `percent_dunn_ADP 25`$Condition, method = "none", list = T)
Sidak(c(percent_dunn_ADP25$P[2], percent_dunn_ADP25$P[6]))
rstatix::cohens_d(data = subset(subset(h2o2_long_percent[-c(200),], Smoke %in% "Control"), State %in% "ADP 25"), formula = Rate ~ Tissue)
rstatix::cohens_d(data = subset(subset(h2o2_long_percent[-c(200),], Smoke %in% "Smoke"), State %in% "ADP 25"), formula = Rate ~ Tissue)

mainplot_percent <- ggplot(data = subset(h2o2_long_percent[-c(200),], State %in% c('ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  # geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 1.2), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.2, 0.2), labels = scales::percent_format()) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "longdash", color = "grey50", size = .6) + 
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(h2o2_long_percent[-c(200),], State %in% c('ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  #ADP25
  annotate('text', x = 0.51, y = 1.4*0.99, label = 'Smoke: p = 0.903', hjust = 0, fontface = 2) +
  annotate('text', x = 0.51, y = 1.4*0.94, label = 'Tissue: p = 0.007', hjust = 0, fontface = 2) +
  annotate('text', x = 0.51, y = 1.4*0.89, label = 'Interaction: p = 0.703', hjust = 0, fontface = 2) +
  #ADP50
  annotate('text', x = 1.51, y = 1.4*0.99, label = 'Smoke: p = 0.249', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 1.4*0.94, label = 'Tissue: p = 0.121', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 1.4*0.89, label = 'Interaction: p = 0.585', hjust = 0, fontface = 2) +
  #ADP100
  annotate('text', x = 2.51, y = 1.4*0.99, label = 'Smoke: p = 0.373', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 1.4*0.94, label = 'Tissue: p = 0.131', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 1.4*0.89, label = 'Interaction: p = 0.558', hjust = 0, fontface = 2) +
  #ADP250
  annotate('text', x = 3.51, y = 1.4*0.99, label = 'Smoke: p = 0.242', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 1.4*0.94, label = 'Tissue: p = 0.457', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 1.4*0.89, label = 'Interaction: p = 0.514', hjust = 0, fontface = 2) +
  #ADP5000
  annotate('text', x = 4.51, y = 1.4*0.99, label = 'Smoke: p = 0.217', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 1.4*0.94, label = 'Tissue: p = 0.158', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 1.4*0.89, label = 'Interaction: p = 0.947', hjust = 0, fontface = 2)


mainplot_percent_Soleus <- ggplot(data = subset(subset(h2o2_long_percent[-c(200),], Tissue %in% "Soleus"), State %in% c('ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, color = "black") +
  theme_prism() +
  scale_fill_manual(breaks = c("Control", "Smoke"), values = c("white", "grey33")) +
  coord_cartesian(ylim = c(0, 1.4), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.4, 0.2), labels = scales::percent_format()) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "longdash", color = "grey50", size = .6) + 
  new_scale_fill() +
  geom_beeswarm(data = subset(subset(h2o2_long_percent[-c(200),], Tissue %in% "Soleus"), State %in% c('ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Smoke), 
                dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T, cex = 1.5) +
  ggtitle(expression(bold(Soleus~H['2']*O['2']*~Emission~('% of GMS')))) +
  scale_fill_manual(values = c("white", "white"))
