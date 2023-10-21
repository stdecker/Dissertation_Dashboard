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

smoke_comparisons <- list(c("Control", "Smoke"))

# Pyruvate data ----
Sol_pyr_data <- read_excel('C:\\Users\\u1159489\\Box\\LayecLab\\CS-THNR\\O2K Analysis Files\\Project 3\\Analysis File3_Pyruvate.xlsx', sheet = "Final Data")

Sol_pyr_data$Tissue[Sol_pyr_data$Tissue == "Gastroc"] <- NA

Sol_pyr_data <- na.omit(Sol_pyr_data)

names(Sol_pyr_data)[2] <- "Smoke"

Sol_pyr_data$Smoke[Sol_pyr_data$Smoke == "Con"] <- "Control"

Sol_pyr_data$Smoke[Sol_pyr_data$Smoke == "Smo"] <- "Smoke"

Sol_pyr_data <- separate(Sol_pyr_data, Subject, into = "Subject", sep = " ")

Sol_pyr_long <- pivot_longer(Sol_pyr_data, cols = 7:15, names_to = "State", values_to = "Rate")

Sol_pyr_mm <- pivot_longer(Sol_pyr_data[1:13], 8:13, values_to = "Rate", names_to = "Concentration")

Sol_pyr_mm$Concentration[Sol_pyr_mm$Concentration == "MD"] <- 0
Sol_pyr_mm$Concentration[Sol_pyr_mm$Concentration == "MDPyr1"] <- 0.1
Sol_pyr_mm$Concentration[Sol_pyr_mm$Concentration == "MDPyr2"] <- 0.25
Sol_pyr_mm$Concentration[Sol_pyr_mm$Concentration == "MDPyr3"] <- 0.5
Sol_pyr_mm$Concentration[Sol_pyr_mm$Concentration == "MDPyr4"] <- 1
Sol_pyr_mm$Concentration[Sol_pyr_mm$Concentration == "MDPyr5"] <- 5

Sol_pyr_mm$Concentration <- as.numeric(Sol_pyr_mm$Concentration)

for(y in unique(Sol_pyr_mm$Smoke)){
  model.drm <- drm(data = subset(Sol_pyr_mm, Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
  
  mml <- data.frame(S = seq(0, max(Sol_pyr_mm$Concentration), length.out = 100))
  
  mml$v <- predict(model.drm, newdata = mml)
  
  modelselect <- mselect(model.drm, fctList = list(MM.2(), MM.3()))
  fit <- modelFit(model.drm)
  
  assign(paste("Sol_pyr_summary", paste(y), sep = "_"), summary(model.drm))
  assign(paste("Sol_pyr_fit", paste(y), sep = "_"), fit)
  assign(paste("Sol_pyr_select", paste(y), sep = "_"), modelselect)
  
  assign(paste("Sol_pyr_model", paste(y), sep = "_"), model.drm)
  assign(paste("Sol_pyr_mml", paste(y), sep = "_"), mml)
}


Sol_pyr_mml_Control$Smoke <- "Control"
Sol_pyr_mml_Smoke$Smoke <- "Smoke"

Sol_mml_pyruvate <- rbind(Sol_pyr_mml_Control, Sol_pyr_mml_Smoke)

Sol_pyr_errors <- Sol_pyr_mm |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

Sol_pyr_errors$min <- Sol_pyr_errors$Rate_mean - (Sol_pyr_errors$Smoke == "Smoke")*Sol_pyr_errors$Rate_se
Sol_pyr_errors$max <- Sol_pyr_errors$Rate_mean + (Sol_pyr_errors$Smoke == "Control")*Sol_pyr_errors$Rate_se

Sol_pyr_kinetic_values <- data.frame()
Sol_pyr_predicted <- data.frame()

for(x in unique(Sol_pyr_mm$Subject)){
  for(y in unique(Sol_pyr_mm$Smoke)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      model.drm <- drm(data = subset(subset(Sol_pyr_mm, Subject %in% paste(x)), Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower","Vmax", "Km")))
      
      coefs <- data.frame("Vmax" = c(model.drm$coefficients[2]), "Km" = c(model.drm$coefficients[3]), "C" = c(model.drm$coefficients[1]))
      
      coefs$Subject <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      mml <- data.frame(S = seq(0, max(Sol_pyr_mm$Concentration), length.out = 100))
      
      mml$v <- predict(model.drm, newdata = mml)
      
      mml$Subject <- paste(x)
      
      mml$Smoke <- paste(y)
      
      Sol_pyr_predicted <- rbind(Sol_pyr_predicted, mml)
      
      coefs$Tissue <- "Soleus"
      coefs$Substrate <- "Pyruvate"
      
      Sol_pyr_kinetic_values <- rbind(Sol_pyr_kinetic_values, coefs)}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    
  }
}

Sol_pyr_kinetic_plot <- ggplot(data = Sol_pyr_mm, aes(x = Concentration, y = Rate, fill = Smoke)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1) +
  geom_line(data = Sol_pyr_predicted, aes(x = S, y = v, linetype = Smoke), size = 1.5, stat = "summary", fun.y = "mean_se()") +
  geom_errorbar(data = Sol_pyr_errors, aes(x = Concentration, y = Rate_mean, colour = Smoke, ymin = min, ymax = max), inherit.aes = F, size = 1, width = 0.1, color = "black") +
  geom_point(stat = "summary", fun.y = "mean", size = 3, pch = 21) +
  theme_prism() +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(2, 1)) +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(0, 70), xlim = c(-.1, 5.1), expand = FALSE) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 70, 10)) +
  #ggtitle("Soleus Pyruvate Kinetics") +
  ggtitle("Soleus") +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt]))), x="[Pyruvate] (mM)")# +
# geom_segment(aes(x = Sol_pyr_model_Control$coefficients[3], xend = Sol_pyr_model_Control$coefficients[3], y = 0, yend = .5*Sol_pyr_model_Control$coefficients[2] + .5*Sol_pyr_model_Control$coefficients[1]), color = "blue", linetype = "dotted") +
# geom_segment(aes(x = Sol_pyr_model_Smoke$coefficients[3], xend = Sol_pyr_model_Smoke$coefficients[3], y = 0, yend = .5*Sol_pyr_model_Smoke$coefficients[2] + .5*Sol_pyr_model_Smoke$coefficients[1]), color = "red", linetype = "dotted") +
# geom_hline(yintercept = Sol_pyr_model_Control$coefficients[2], color = "blue", linetype = "longdash") +
# geom_hline(yintercept = Sol_pyr_model_Smoke$coefficients[2], color = "red", linetype = "longdash")

Sol_pyr_vmax <- ggplot(data = Sol_pyr_kinetic_values, aes(x = Smoke, y = Vmax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 90)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(10, 90, 20)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 71, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1, tip.length = 0)


Sol_pyr_km <- ggplot(data = subset(Sol_pyr_kinetic_values, Km <0.4), aes(x = Smoke, y = Km, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 0.55)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.05, 0.55, 0.05)) +
  labs(y=expression(bold(atop('[Pyruvate]', (mM))))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .43, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1, tip.length = 0)


# Palmitate Data ----
Sol_palm_data <- read_excel('C:\\Users\\u1159489\\Box\\LayecLab\\CS-THNR\\O2K Analysis Files\\Project 3\\Analysis File3_Palmitate.xlsx', sheet = "Final Data")

Sol_palm_data$Tissue[Sol_palm_data$Tissue == "Gastroc"] <- NA

Sol_palm_data <- na.omit(Sol_palm_data)

Sol_palm_data <- separate(Sol_palm_data, Subject, into = "Subject", sep = " ")

names(Sol_palm_data)[2] <- "Smoke"

Sol_palm_data$Smoke[Sol_palm_data$Smoke == "Con"] <- "Control"

Sol_palm_data$Smoke[Sol_palm_data$Smoke == "Smo"] <- "Smoke"

## Palmitate -----
palm <- pivot_longer(Sol_palm_data, cols = 7:20, names_to = "State", values_to = "Rate")

Sol_palm_mm <- pivot_longer(Sol_palm_data[1:13], 8:13, values_to = "Rate", names_to = "Concentration")

Sol_palm_mm$Concentration[Sol_palm_mm$Concentration == "MD"] <- 0
Sol_palm_mm$Concentration[Sol_palm_mm$Concentration == "MDPalm1"] <- 0.0025
Sol_palm_mm$Concentration[Sol_palm_mm$Concentration == "MDPalm2"] <- 0.005
Sol_palm_mm$Concentration[Sol_palm_mm$Concentration == "MDPalm3"] <- 0.0125
Sol_palm_mm$Concentration[Sol_palm_mm$Concentration == "MDPalm4"] <- 0.025
Sol_palm_mm$Concentration[Sol_palm_mm$Concentration == "MDPalm5"] <- 0.04

Sol_palm_mm$Concentration <- as.numeric(Sol_palm_mm$Concentration)

for(y in unique(Sol_palm_mm$Smoke)){
  model.drm <- drm(data = subset(Sol_palm_mm, Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower","Vmax", "Km")))
  
  mml <- data.frame(S = seq(0, max(Sol_palm_mm$Concentration), length.out = 100))
  
  mml$v <- predict(model.drm, newdata = mml)
  
  modelselect <- mselect(model.drm, fctList = list(MM.2(), MM.3()))
  fit <- modelFit(model.drm)
  
  assign(paste("Sol_palm_summary", paste(y), sep = "_"), summary(model.drm))
  assign(paste("Sol_palm_fit", paste(y), sep = "_"), fit)
  assign(paste("Sol_palm_select", paste(y), sep = "_"), modelselect)
  
  assign(paste("Sol_palm_model", paste(y), sep = "_"), model.drm)
  assign(paste("Sol_palm_mml", paste(y), sep = "_"), mml)
}

Sol_palm_mml_Control$Smoke <- "Control"
Sol_palm_mml_Smoke$Smoke <- "Smoke"

Sol_mml_palm <- rbind(Sol_palm_mml_Control, Sol_palm_mml_Smoke)

Sol_palm_errors <- Sol_palm_mm |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

Sol_palm_errors$min <- Sol_palm_errors$Rate_mean + (Sol_palm_errors$Smoke == "Smoke")*Sol_palm_errors$Rate_se
Sol_palm_errors$max <- Sol_palm_errors$Rate_mean - (Sol_palm_errors$Smoke == "Control")*Sol_palm_errors$Rate_se

Sol_palm_kinetic_values <- data.frame()
Sol_palm_predicted <- data.frame()

for(x in unique(Sol_palm_mm$Subject)){
  for(y in unique(Sol_palm_mm$Smoke)){
    
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      
      model.drm <- drm(data = subset(subset(Sol_palm_mm, Subject %in% paste(x)), Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
      
      coefs <- data.frame("Vmax" = c(model.drm$coefficients[2]), "Km" = c(model.drm$coefficients[3]), "C" = c(model.drm$coefficients[1]))
      
      coefs$Subject <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      mml <- data.frame(S = seq(0, max(Sol_palm_mm$Concentration), length.out = 100))
      
      mml$v <- predict(model.drm, newdata = mml)
      
      mml$Subject <- paste(x)
      
      mml$Smoke <- paste(y)
      
      Sol_palm_predicted <- rbind(Sol_palm_predicted, mml)
      
      coefs$Tissue <- "Soleus"
      coefs$Substrate <- "PC"
      
      Sol_palm_kinetic_values <- rbind(Sol_palm_kinetic_values, coefs)}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    
  }
}

Sol_palm_kinetic_plot <- ggplot(data = Sol_palm_mm, aes(x = Concentration, y = Rate, fill = Smoke)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.001) +
  geom_line(data = Sol_palm_predicted, aes(x = S, y = v, linetype = Smoke), size = 1.5, stat = "summary", fun.y = "mean_se()") +
  geom_errorbar(data = Sol_palm_errors, aes(x = Concentration, y = Rate_mean, colour = Smoke, ymin = min, ymax = max), inherit.aes = F, size = 1, width = 0.00065, color = "black") +
  geom_point(stat = "summary", fun.y = "mean", size = 3, pch = 21) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(2, 1)) +
  theme_prism() +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(0, 25), xlim = c(-.001, 0.041), expand = FALSE)+
  ggtitle("Soleus") +
  #ggtitle("Soleus PC") +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt]))), x="[PC] (mM)")# +
# geom_segment(aes(x = Sol_palm_model_Control$coefficients[3], xend = Sol_palm_model_Control$coefficients[3], y = 0, yend = .5*Sol_palm_model_Control$coefficients[2] + .5*Sol_palm_model_Control$coefficients[1]), color = "blue", linetype = "dotted") +
# geom_segment(aes(x = Sol_palm_model_Smoke$coefficients[3], xend = Sol_palm_model_Smoke$coefficients[3], y = 0, yend = .5*Sol_palm_model_Smoke$coefficients[2] + .5*Sol_palm_model_Smoke$coefficients[1]), color = "red", linetype = "dotted") +
# geom_hline(yintercept = Sol_palm_model_Control$coefficients[2], color = "blue", linetype = "longdash") +
# geom_hline(yintercept = Sol_palm_model_Smoke$coefficients[2], color = "red", linetype = "longdash")

Sol_palm_vmax <- ggplot(data = Sol_palm_kinetic_values[-c(9),], aes(x = Smoke, y = Vmax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 35)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(5, 60, 5)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 27, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1, tip.length = 0)

Sol_palm_km <- ggplot(data = subset(Sol_palm_kinetic_values, Km < 0.065), aes(x = Smoke, y = Km, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 0.09)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.0, 1.5, 0.01)) +
  labs(y=expression(bold(atop('[PC]', (mM))))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 0.068, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1, tip.length = 0)

## Pyruvate in Palmitate ----

Sol_palm_pyr <- pivot_longer(Sol_palm_data, cols = 7:15, names_to = "State", values_to = "Rate")

Sol_palm_pyr_mm <- pivot_longer(Sol_palm_data[c(1:6, 13:18)], 7:12, values_to = "Rate", names_to = "Concentration")

Sol_palm_pyr_mm$Concentration[Sol_palm_pyr_mm$Concentration == "MDPalm5"] <- 0
Sol_palm_pyr_mm$Concentration[Sol_palm_pyr_mm$Concentration == "MDPalmPyr1"] <- 0.1
Sol_palm_pyr_mm$Concentration[Sol_palm_pyr_mm$Concentration == "MDPalmPyr2"] <- 0.25
Sol_palm_pyr_mm$Concentration[Sol_palm_pyr_mm$Concentration == "MDPalmPyr3"] <- 0.5
Sol_palm_pyr_mm$Concentration[Sol_palm_pyr_mm$Concentration == "MDPalmPyr4"] <- 1
Sol_palm_pyr_mm$Concentration[Sol_palm_pyr_mm$Concentration == "MDPalmPyr5"] <- 5

Sol_palm_pyr_mm$Concentration <- as.numeric(Sol_palm_pyr_mm$Concentration)

for(y in unique(Sol_palm_pyr_mm$Smoke)){
  model.drm <- drm(data = subset(Sol_palm_pyr_mm[-c(6, 18),], Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
  
  mml <- data.frame(S = seq(0, max(Sol_palm_pyr_mm$Concentration), length.out = 100))
  
  mml$v <- predict(model.drm, newdata = mml)
  
  modelselect <- mselect(model.drm, fctList = list(MM.2(), MM.3()))
  fit <- modelFit(model.drm)
  
  assign(paste("Sol_palm_pyr_summary", paste(y), sep = "_"), summary(model.drm))
  assign(paste("Sol_palm_pyr_fit", paste(y), sep = "_"), fit)
  assign(paste("Sol_palm_pyr_select", paste(y), sep = "_"), modelselect)
  
  assign(paste("Sol_palm_pyr_model", paste(y), sep = "_"), model.drm)
  assign(paste("Sol_palm_pyr_mml", paste(y), sep = "_"), mml)
}

Sol_palm_pyr_mml_Control$Smoke <- "Control"
Sol_palm_pyr_mml_Smoke$Smoke <- "Smoke"

mml_Sol_palm_pyruvate <- rbind(Sol_palm_pyr_mml_Control, Sol_palm_pyr_mml_Smoke)

Sol_palm_pyr_errors <- Sol_palm_pyr_mm |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

Sol_palm_pyr_errors$min <- Sol_palm_pyr_errors$Rate_mean - (Sol_palm_pyr_errors$Smoke == "Smoke")*Sol_palm_pyr_errors$Rate_se
Sol_palm_pyr_errors$max <- Sol_palm_pyr_errors$Rate_mean + (Sol_palm_pyr_errors$Smoke == "Control")*Sol_palm_pyr_errors$Rate_se
  
Sol_palm_pyr_kinetic_values <- data.frame()
Sol_palm_pyr_predicted <- data.frame()

for(x in unique(Sol_palm_pyr_mm$Subject)){
  for(y in unique(Sol_palm_pyr_mm$Smoke)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      model.drm <- drm(data = subset(subset(Sol_palm_pyr_mm[-c(6, 18, 66),], Subject %in% paste(x)), Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
      
      coefs <- data.frame("Vmax" = c(model.drm$coefficients[2]), "Km" = c(model.drm$coefficients[3]), "C" = c(model.drm$coefficients[1]))
      
      coefs$Subject <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      mml <- data.frame(S = seq(0, max(Sol_palm_pyr_mm$Concentration), length.out = 100))
      
      mml$v <- predict(model.drm, newdata = mml)
      
      mml$Subject <- paste(x)
      
      mml$Smoke <- paste(y)
      
      Sol_palm_pyr_predicted <- rbind(Sol_palm_pyr_predicted, mml)
      
      coefs$Tissue <- "Soleus"
      coefs$Substrate <- "Pyruvate in PC"
      
      Sol_palm_pyr_kinetic_values <- rbind(Sol_palm_pyr_kinetic_values, coefs)}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    
  }
}

Sol_palm_pyr_kinetic_plot <- ggplot(data = Sol_palm_pyr_mm, aes(x = Concentration, y = Rate, fill = Smoke)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1) +
  geom_line(data = Sol_palm_pyr_predicted, aes(x = S, y = v, linetype = Smoke), size = 1.5, stat = "summary", fun.y = "mean_se()") +
  geom_errorbar(data = Sol_palm_pyr_errors, aes(x = Concentration, y = Rate_mean, colour = Smoke, ymin = min, ymax = max), inherit.aes = F, size = 1, width = 0.1, color = "black") +
  geom_point(stat = "summary", fun.y = "mean", size = 3, pch = 21) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(2, 1)) +
  theme_prism() +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(0, 70), xlim = c(-.1, 5.1), expand = FALSE) +
  scale_y_continuous(breaks = seq(0,70,10)) +
  ggtitle("Soleus") +
  #ggtitle("Soleus Pyruvate in PC (0.04 mM)") +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt]))), x="[Pyuvate] (mM)")# +
# geom_segment(aes(x = Sol_palm_pyr_model_Control$coefficients[3], xend = Sol_palm_pyr_model_Control$coefficients[3], y = 0, yend = .5*Sol_palm_pyr_model_Control$coefficients[2] + .5*Sol_palm_pyr_model_Control$coefficients[1]), color = "blue", linetype = "dotted") +
# geom_segment(aes(x = Sol_palm_pyr_model_Smoke$coefficients[3], xend = Sol_palm_pyr_model_Smoke$coefficients[3], y = 0, yend = .5*Sol_palm_pyr_model_Smoke$coefficients[2] + .5*Sol_palm_pyr_model_Smoke$coefficients[1]), color = "red", linetype = "dotted") +
# geom_hline(yintercept = Sol_palm_pyr_model_Control$coefficients[2], color = "blue", linetype = "longdash") +
# geom_hline(yintercept = Sol_palm_pyr_model_Smoke$coefficients[2], color = "red", linetype = "longdash")

# for(z in names(Sol_palm_pyr_kinetic_values[1:2])){
Sol_palm_pyr_vmax <- ggplot(data = Sol_palm_pyr_kinetic_values, aes(x = Smoke, y = Vmax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 80)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(10, 80, 10)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 68, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1, tip.length = 0)

Sol_palm_pyr_km <- ggplot(data = Sol_palm_pyr_kinetic_values, aes(x = Smoke, y = Km, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 0.35)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.0, 0.35, 0.05)) +
  labs(y=expression(bold(atop('[Pyruvate]',(mM))))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .26, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1, tip.length = 0)


# ggstatsplot::ggbetweenstats(data = Sol_palm_pyr_kinetic_values, x = Smoke, y = Km, outlier.tagging = T)
# ggstatsplot::ggbetweenstats(data = Sol_palm_pyr_kinetic_values, x = Smoke, y = Vmax, outlier.tagging = T)
# 
# ggstatsplot::ggbetweenstats(data = subset(Sol_pyr_kinetic_values, Km <0.4), x = Smoke, y = Km, outlier.tagging = T)
# ggstatsplot::ggbetweenstats(data = Sol_pyr_kinetic_values, x = Smoke, y = Vmax, outlier.tagging = T)
# 
# ggstatsplot::ggbetweenstats(data = subset(Sol_palm_kinetic_values, Km < 0.069), x = Smoke, y = Km, outlier.tagging = T)
# ggstatsplot::ggbetweenstats(data = Sol_palm_kinetic_values[-c(9),], x = Smoke, y = Vmax, outlier.tagging = T)



# Bar Graphs ----
## Pyruvate data ----
Sol_pyr_vmax <- ggplot(data = Sol_pyr_kinetic_values, aes(x = Smoke, y = Vmax, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 100, 20)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 89, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0)


Sol_pyr_km <- ggplot(data = subset(Sol_pyr_kinetic_values), aes(x = Smoke, y = Km, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.7)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.7, 0.1)) +
  labs(y=expression(bold(atop('[Pyruvate]', (mM))))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .55, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0, paired = T)


## Palmitate -----
Sol_palm_vmax <- ggplot(data = Sol_palm_kinetic_values, aes(x = Smoke, y = Vmax, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 50)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 60, 10)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 43, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0)

Sol_palm_km <- ggplot(data = subset(Sol_palm_kinetic_values, Km < 0.069), aes(x = Smoke, y = Km, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.07)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.0, 1.5, 0.01)) +
  labs(y=expression(bold(atop('[PC]', (mM))))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 0.052, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0)

## Pyruvate in Palmitate ----
Sol_palm_pyr_vmax <- ggplot(data = Sol_palm_pyr_kinetic_values, aes(x = Smoke, y = Vmax, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 120)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 95, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0)

Sol_palm_pyr_km <- ggplot(data = Sol_palm_pyr_kinetic_values, aes(x = Smoke, y = Km, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.35)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.0, 0.35, 0.05)) +
  labs(y=expression(bold(atop('[Pyruvate]',(mM))))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .295, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, tip.length = 0)



pyr_ordered <- Sol_pyr_kinetic_values |> 
  arrange(Smoke)

palm_pyr_ordered <- Sol_palm_pyr_kinetic_values |> 
  arrange(Smoke)

palm_ordered <- Sol_palm_kinetic_values |> 
  arrange(Smoke)


Sol_percent <- data.frame("Km" = palm_pyr_ordered$Km/pyr_ordered$Km,
                               "Vmax" = palm_pyr_ordered$Vmax/pyr_ordered$Vmax,
                               "Subject" = pyr_ordered$Subject,
                               "Smoke" = pyr_ordered$Smoke,
                               "Tissue" = pyr_ordered$Tissue)

Sol_delta <- data.frame("Km" = palm_pyr_ordered$Km-pyr_ordered$Km,
                          "Vmax" = palm_pyr_ordered$Vmax-pyr_ordered$Vmax,
                          "Subject" = pyr_ordered$Subject,
                          "Smoke" = pyr_ordered$Smoke,
                          "Tissue" = pyr_ordered$Tissue)


Sol_contributions <- data.frame("Pyruvate" = palm_pyr_ordered$Vmax/palm_pyr_ordered$Vmax,
                          "PC" = palm_ordered$Vmax/palm_pyr_ordered$Vmax,
                          "Subject" = pyr_ordered$Subject,
                          "Smoke" = pyr_ordered$Smoke,
                          "Tissue" = pyr_ordered$Tissue)
