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
library(ggbreak)

smoke_comparisons <- list(c("Control", "Smoke"))


# Pyruvate data ----
Gast_pyr_data <- read_excel('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Project 3\\Analysis File3_Pyruvate.xlsx', sheet = "Final Data")

Gast_pyr_data$Tissue[Gast_pyr_data$Tissue == "Soleus"] <- NA

Gast_pyr_data <- na.omit(Gast_pyr_data)

names(Gast_pyr_data)[2] <- "Smoke"

Gast_pyr_data$Smoke[Gast_pyr_data$Smoke == "Con"] <- "Control"

Gast_pyr_data$Smoke[Gast_pyr_data$Smoke == "Smo"] <- "Smoke"

Gast_pyr_data <- separate(Gast_pyr_data, Subject, into = "Subject", sep = " ")

Gast_pyr_long <- pivot_longer(Gast_pyr_data, cols = 7:15, names_to = "State", values_to = "Rate")

Gast_pyr_mm <- pivot_longer(Gast_pyr_data[1:13], 8:13, values_to = "Rate", names_to = "Concentration")

Gast_pyr_mm$Concentration[Gast_pyr_mm$Concentration == "MD"] <- 0
Gast_pyr_mm$Concentration[Gast_pyr_mm$Concentration == "MDPyr1"] <- 0.1
Gast_pyr_mm$Concentration[Gast_pyr_mm$Concentration == "MDPyr2"] <- 0.25
Gast_pyr_mm$Concentration[Gast_pyr_mm$Concentration == "MDPyr3"] <- 0.5
Gast_pyr_mm$Concentration[Gast_pyr_mm$Concentration == "MDPyr4"] <- 1
Gast_pyr_mm$Concentration[Gast_pyr_mm$Concentration == "MDPyr5"] <- 5

Gast_pyr_mm$Concentration <- as.numeric(Gast_pyr_mm$Concentration)

for(y in unique(Gast_pyr_mm$Smoke)){
  model.drm <- drm(data = subset(Gast_pyr_mm, Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
  
  mml <- data.frame(S = seq(0, max(Gast_pyr_mm$Concentration), length.out = 100))
  
  mml$v <- predict(model.drm, newdata = mml)
  
  modelselect <- mselect(model.drm, fctList = list(MM.2(), MM.3()))
  fit <- modelFit(model.drm)
  
  assign(paste("Gast_pyr_summary", paste(y), sep = "_"), summary(model.drm))
  assign(paste("Gast_pyr_fit", paste(y), sep = "_"), fit)
  assign(paste("Gast_pyr_select", paste(y), sep = "_"), modelselect)
  assign(paste("Gast_pyr_model", paste(y), sep = "_"), model.drm)
  assign(paste("Gast_pyr_mml", paste(y), sep = "_"), mml)
}

Gast_pyr_stats <- psych::describeBy(Gast_pyr_data, Gast_pyr_data$Smoke)

Gast_pyr_mml_Control$Smoke <- "Control"
Gast_pyr_mml_Smoke$Smoke <- "Smoke"

Gast_mml_pyruvate <- rbind(Gast_pyr_mml_Control, Gast_pyr_mml_Smoke)

Gast_pyr_errors <- Gast_pyr_mm |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

Gast_pyr_errors$min <- Gast_pyr_errors$Rate_mean - (Gast_pyr_errors$Smoke == "Smoke")*Gast_pyr_errors$Rate_se
Gast_pyr_errors$max <- Gast_pyr_errors$Rate_mean + (Gast_pyr_errors$Smoke == "Control")*Gast_pyr_errors$Rate_se

Gast_pyr_kinetic_values <- data.frame()
Gast_pyr_predicted <- data.frame()

for(x in unique(Gast_pyr_mm$Subject)){
  for(y in unique(Gast_pyr_mm$Smoke)){
    
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      
      model.drm <- drm(data = subset(subset(Gast_pyr_mm, Subject %in% paste(x)), Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
      
      coefs <- data.frame("Vmax" = c(model.drm$coefficients[2]), "Km" = c(model.drm$coefficients[3]), "C" = c(model.drm$coefficients[1]))
      
      coefs$Subject <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      mml <- data.frame(S = seq(0, max(Gast_pyr_mm$Concentration), length.out = 100))
      
      mml$v <- predict(model.drm, newdata = mml)
      
      mml$Subject <- paste(x)
      
      mml$Smoke <- paste(y)
      
      Gast_pyr_predicted <- rbind(Gast_pyr_predicted, mml)
      
      coefs$Tissue <- "Gastrocnemius"
      coefs$Substrate <- "Pyruvate"
      
      Gast_pyr_kinetic_values <- rbind(Gast_pyr_kinetic_values, coefs)}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    
  }
}

Gast_pyr_kinetic_plot <- ggplot(data = Gast_pyr_mm, aes(x = Concentration, y = Rate, fill = Smoke)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.1) +
  geom_line(data = Gast_pyr_predicted, aes(x = S, y = v, linetype = Smoke), size = 1.5, stat = "summary", fun.y = "mean_se()") +
  geom_errorbar(data = Gast_pyr_errors, aes(x = Concentration, y = Rate_mean, colour = Smoke, ymin = min, ymax = max), inherit.aes = F, size = 1, width = 0.1, color = "black") +
  geom_point(stat = "summary", fun.y = "mean", size = 3, pch = 21) +
  theme_prism() +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(2, 1)) +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(0, 70), xlim = c(-.1, 5.1), expand = FALSE) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 70, 10)) +
  ggtitle("Gastrocnemius") +
  #ggtitle(paste0("Gastrocnemius Pyruvate")) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt]))), x="[Pyruvate] (mM)")# +
# geom_segment(aes(x = Gast_pyr_model_Control$coefficients[3], xend = Gast_pyr_model_Control$coefficients[3], y = 0, yend = .5*Gast_pyr_model_Control$coefficients[2] + .5*Gast_pyr_model_Control$coefficients[1]), color = "blue", linetype = "dotted") +
# geom_segment(aes(x = Gast_pyr_model_Smoke$coefficients[3], xend = Gast_pyr_model_Smoke$coefficients[3], y = 0, yend = .5*Gast_pyr_model_Smoke$coefficients[2] + .5*Gast_pyr_model_Smoke$coefficients[1]), color = "red", linetype = "dotted") +
# geom_hline(yintercept = Gast_pyr_model_Control$coefficients[2], color = "blue", linetype = "longdash") +
# geom_hline(yintercept = Gast_pyr_model_Smoke$coefficients[2], color = "red", linetype = "longdash")

Gast_pyr_vmax <- ggplot(data = Gast_pyr_kinetic_values, aes(x = Smoke, y = Vmax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 70)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 70, 10)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 58, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1)


Gast_pyr_km <- ggplot(data = Gast_pyr_kinetic_values, aes(x = Smoke, y = Km, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  ylab(expression(bold(atop('[Pyruvate]', '(mM)')))) +
  ggtitle(expression(bold(K[m]))) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.1, 10, 0.1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .61, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1)


# Palmitate Data ----
Gast_palm_data <- read_excel('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Project 3\\Analysis File3_Palmitate.xlsx', sheet = "Final Data")

Gast_palm_data$Tissue[Gast_palm_data$Tissue == "Soleus"] <- NA  

Gast_palm_data <- na.omit(Gast_palm_data)

Gast_palm_data <- separate(Gast_palm_data, Subject, into = "Subject", sep = " ")

names(Gast_palm_data)[2] <- "Smoke"

Gast_palm_data$Smoke[Gast_palm_data$Smoke == "Con"] <- "Control"

Gast_palm_data$Smoke[Gast_palm_data$Smoke == "Smo"] <- "Smoke"

## Palmitate -----
palm <- pivot_longer(Gast_palm_data, cols = 7:20, names_to = "State", values_to = "Rate")

Gast_palm_mm <- pivot_longer(Gast_palm_data[1:13], 8:13, values_to = "Rate", names_to = "Concentration")

Gast_palm_mm$Concentration[Gast_palm_mm$Concentration == "MD"] <- 0
Gast_palm_mm$Concentration[Gast_palm_mm$Concentration == "MDPalm1"] <- 0.0025
Gast_palm_mm$Concentration[Gast_palm_mm$Concentration == "MDPalm2"] <- 0.005
Gast_palm_mm$Concentration[Gast_palm_mm$Concentration == "MDPalm3"] <- 0.0125
Gast_palm_mm$Concentration[Gast_palm_mm$Concentration == "MDPalm4"] <- 0.025
Gast_palm_mm$Concentration[Gast_palm_mm$Concentration == "MDPalm5"] <- 0.04

Gast_palm_mm$Concentration <- as.numeric(Gast_palm_mm$Concentration)

for(y in unique(Gast_palm_mm$Smoke)){
  
  skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  
  tryCatch({
    model.drm <- drm(data = subset(Gast_palm_mm[-c(12, 24, 16, 17, 18, 36),], Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
    
    mml <- data.frame(S = seq(0, max(Gast_palm_mm$Concentration), length.out = 100))
    
    mml$v <- predict(model.drm, newdata = mml)
    
    modelselect <- mselect(model.drm, fctList = list(MM.2(), MM.3()))
    fit <- modelFit(model.drm)
    
    assign(paste("Gast_palm_summary", paste(y), sep = "_"), summary(model.drm))
    assign(paste("Gast_palm_fit", paste(y), sep = "_"), fit)
    assign(paste("Gast_palm_select", paste(y), sep = "_"), modelselect)
    
    assign(paste("Gast_palm_model", paste(y), sep = "_"), model.drm)
    assign(paste("Gast_palm_mml", paste(y), sep = "_"), mml)}, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
}

# model.drm <- drm(data = subset(Gast_palm_mm, Smoke %in% "Control"), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
# 
# mml <- data.frame(S = seq(0, max(Gast_palm_mm$Concentration), length.out = 100))
# 
# mml$v <- predict(model.drm, newdata = mml)
# 
# Gast_palm_mml_Control <- mml
# 
# 
# 
# model.drm2 <- drm(data = subset(Gast_palm_mm, Smoke %in% "Smoke"), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
# 
# mml2 <- data.frame(S = seq(0, max(Gast_palm_mm$Concentration), length.out = 100))
# 
# mml2$v <- predict(model.drm2, newdata = mml2)
# 
# Gast_palm_mml_Smoke <- mml2


Gast_palm_mml_Control$Smoke <- "Control"
Gast_palm_mml_Smoke$Smoke <- "Smoke"

Gast_mml_palm <- rbind(Gast_palm_mml_Control, Gast_palm_mml_Smoke)

Gast_palm_errors <- Gast_palm_mm |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

Gast_palm_errors$min <- Gast_palm_errors$Rate_mean - (Gast_palm_errors$Smoke == "Smoke")*Gast_palm_errors$Rate_se
Gast_palm_errors$max <- Gast_palm_errors$Rate_mean + (Gast_palm_errors$Smoke == "Control")*Gast_palm_errors$Rate_se

  
Gast_palm_kinetic_values <- data.frame()
Gast_palm_predicted <- data.frame()

for(x in unique(Gast_palm_mm$Subject)){
  for(y in unique(Gast_palm_mm$Smoke)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      
      
      model.drm <- drm(data = subset(subset(Gast_palm_mm[-c(12, 24, 16, 17, 18, 36),], Subject %in% paste(x)), Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
      
      coefs <- data.frame("Vmax" = c(model.drm$coefficients[2]), "Km" = c(model.drm$coefficients[3]), "C" = c(model.drm$coefficients[1]))
      
      coefs$Subject <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      mml <- data.frame(S = seq(0, max(Gast_palm_mm$Concentration), length.out = 100))
      
      mml$v <- predict(model.drm, newdata = mml)
      
      mml$Subject <- paste(x)
      
      mml$Smoke <- paste(y)
      
      Gast_palm_predicted <- rbind(Gast_palm_predicted, mml)
      
      coefs$Tissue <- "Gastrocnemius"
      coefs$Substrate <- "PC"
      
      Gast_palm_kinetic_values <- rbind(Gast_palm_kinetic_values, coefs)}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    
  }
}

Gast_palm_kinetic_plot <- ggplot(data = Gast_palm_mm, aes(x = Concentration, y = Rate, fill = Smoke)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.001) +
  geom_line(data = Gast_palm_predicted, aes(x = S, y = v, linetype = Smoke), size = 1.5, stat = "summary", fun.y = "mean_se()") +
  geom_errorbar(data = Gast_palm_errors, aes(x = Concentration, y = Rate_mean, colour = Smoke, ymin = min, ymax = max), inherit.aes = F, size = 1, width = 0.00065, color = "black") +
  geom_point(stat = "summary", fun.y = "mean", size = 3, pch = 21) +
  theme_prism() +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(2, 1)) +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(0, 25), xlim = c(-.001, 0.041), expand = FALSE)+
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 25, 5)) +
  ggtitle("Gastrocnemius") +
  #ggtitle(paste0("Gastrocnemius PC")) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt]))), x="[PC] (mM)")# +
# geom_segment(aes(x = Gast_palm_model_Control$coefficients[3], xend = Gast_palm_model_Control$coefficients[3], y = 0, yend = .5*Gast_palm_model_Control$coefficients[2] + .5*Gast_palm_model_Control$coefficients[1]), color = "blue", linetype = "dotted") +
# geom_segment(aes(x = Gast_palm_model_Smoke$coefficients[3], xend = Gast_palm_model_Smoke$coefficients[3], y = 0, yend = .5*Gast_palm_model_Smoke$coefficients[2] + .5*Gast_palm_model_Smoke$coefficients[1]), color = "red", linetype = "dotted") +
# geom_hline(yintercept = Gast_palm_model_Control$coefficients[2], color = "blue", linetype = "longdash") +
# geom_hline(yintercept = Gast_palm_model_Smoke$coefficients[2], color = "red", linetype = "longdash")


Gast_palm_vmax <- ggplot(data = subset(Gast_palm_kinetic_values[-c(8),]), aes(x = Smoke, y = Vmax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  labs(y=expression(bold(V[max]~Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 20)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 20, 5)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 17, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1)


Gast_palm_km <- ggplot(data = subset(Gast_palm_kinetic_values, Km <0.07), aes(x = Smoke, y = Km, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  ylab(expression(bold(K[m]~(mM)))) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 0.04)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.01, 1, 0.01)) +
  labs(y=expression(bold(atop('[PC]', (mM))))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .028, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1)


## Pyruvate in Palmitate ----

Gast_palm_pyr <- pivot_longer(Gast_palm_data, cols = 7:15, names_to = "State", values_to = "Rate")

Gast_palm_pyr_mm <- pivot_longer(Gast_palm_data[c(1:6, 13:18)], 7:12, values_to = "Rate", names_to = "Concentration")

Gast_palm_pyr_mm$Concentration[Gast_palm_pyr_mm$Concentration == "MDPalm5"] <- 0
Gast_palm_pyr_mm$Concentration[Gast_palm_pyr_mm$Concentration == "MDPalmPyr1"] <- 0.1
Gast_palm_pyr_mm$Concentration[Gast_palm_pyr_mm$Concentration == "MDPalmPyr2"] <- 0.25
Gast_palm_pyr_mm$Concentration[Gast_palm_pyr_mm$Concentration == "MDPalmPyr3"] <- 0.5
Gast_palm_pyr_mm$Concentration[Gast_palm_pyr_mm$Concentration == "MDPalmPyr4"] <- 1
Gast_palm_pyr_mm$Concentration[Gast_palm_pyr_mm$Concentration == "MDPalmPyr5"] <- 5

Gast_palm_pyr_mm$Concentration <- as.numeric(Gast_palm_pyr_mm$Concentration)

for(y in unique(Gast_palm_pyr_mm$Smoke)){
  model.drm <- drm(data = subset(Gast_palm_pyr_mm, Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
  
  mml <- data.frame(S = seq(0, max(Gast_palm_pyr_mm$Concentration), length.out = 100))
  
  mml$v <- predict(model.drm, newdata = mml)
  
  modelselect <- mselect(model.drm, fctList = list(MM.2(), MM.3()))
  fit <- modelFit(model.drm)
  
  assign(paste("Gast_palm_pyr_summary", paste(y), sep = "_"), summary(model.drm))
  assign(paste("Gast_palm_pyr_fit", paste(y), sep = "_"), fit)
  assign(paste("Gast_palm_pyr_select", paste(y), sep = "_"), modelselect)
  
  assign(paste("Gast_palm_pyr_model", paste(y), sep = "_"), model.drm)
  assign(paste("Gast_palm_pyr_mml", paste(y), sep = "_"), mml)
}

Gast_palm_pyr_mml_Control$Smoke <- "Control"
Gast_palm_pyr_mml_Smoke$Smoke <- "Smoke"

summary(model.drm)

mml_Gast_palm_pyr <- rbind(Gast_palm_pyr_mml_Control, Gast_palm_pyr_mml_Smoke)

Gast_palm_pyr_errors <- Gast_palm_pyr_mm |>
  group_by(Tissue, Smoke, Concentration) |>
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

Gast_palm_pyr_errors$min <- Gast_palm_pyr_errors$Rate_mean - (Gast_palm_pyr_errors$Smoke == "Smoke")*Gast_palm_pyr_errors$Rate_se
Gast_palm_pyr_errors$max <- Gast_palm_pyr_errors$Rate_mean + (Gast_palm_pyr_errors$Smoke == "Control")*Gast_palm_pyr_errors$Rate_se

Gast_palm_pyr_kinetic_values <- data.frame()
Gast_palm_pyr_predicted <- data.frame()

for(x in unique(Gast_palm_pyr_mm$Subject)){
  for(y in unique(Gast_palm_pyr_mm$Smoke)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      model.drm <- drm(data = subset(subset(Gast_palm_pyr_mm[-c(12, 17, 18, 66),], Subject %in% paste(x)), Smoke %in% paste(y)), formula = Rate ~ Concentration, fct = MM.3(names = c("Lower", "Vmax", "Km")))
      
      coefs <- data.frame("Vmax" = c(model.drm$coefficients[2]), "Km" = c(model.drm$coefficients[3]), "C" = c(model.drm$coefficients[1]))
      
      coefs$Subject <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      mml <- data.frame(S = seq(0, max(Gast_palm_pyr_mm$Concentration), length.out = 100))
      
      mml$v <- predict(model.drm, newdata = mml)
      
      mml$Subject <- paste(x)
      
      mml$Smoke <- paste(y)
      
      Gast_palm_pyr_predicted <- rbind(Gast_palm_pyr_predicted, mml)
      
      coefs$Tissue <- "Gastrocnemius"
      coefs$Substrate <- "Pyruvate in PC"
      
      Gast_palm_pyr_kinetic_values <- rbind(Gast_palm_pyr_kinetic_values, coefs)}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    
  }
}

Gast_palm_pyr_kinetic_plot <- ggplot(data = Gast_palm_pyr_mm, aes(x = Concentration, y = Rate, fill = Smoke)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.1) +
  geom_line(data = Gast_palm_pyr_predicted, aes(x = S, y = v, linetype = Smoke), size = 1.5, stat = "summary", fun.y = "mean_se()") +
  geom_errorbar(data = Gast_palm_pyr_errors, aes(x = Concentration, y = Rate_mean, colour = Smoke, ymin = min, ymax = max), inherit.aes = F, size = 1, width = 0.1, color = "black") +
  geom_point(stat = "summary", fun.y = "mean", size = 3, pch = 21) +
  theme_prism() +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(2, 1)) +
  scale_fill_manual(values = c("white", "grey25")) +
  coord_cartesian(ylim = c(0, 70), xlim = c(-.1, 5.1), expand = FALSE) +
  scale_y_continuous(breaks = seq(0, 70, 10)) +
  ggtitle("Gastrocnemius") +
  #ggtitle(paste0("Gastrocnemius Pyruvate in PC (0.04 mM)")) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt]))), x= "[Pyruvate] (mM)")# +
# geom_segment(aes(x = Gast_palm_pyr_model_Control$coefficients[3], xend = Gast_palm_pyr_model_Control$coefficients[3], y = 0, yend = .5*Gast_palm_pyr_model_Control$coefficients[2] + .5*Gast_palm_pyr_model_Control$coefficients[1]), color = "blue", linetype = "dotted") +
# geom_segment(aes(x = Gast_palm_pyr_model_Smoke$coefficients[3], xend = Gast_palm_pyr_model_Smoke$coefficients[3], y = 0, yend = .5*Gast_palm_pyr_model_Smoke$coefficients[2] + .5*Gast_palm_pyr_model_Smoke$coefficients[1]), color = "red", linetype = "dotted") +
# geom_hline(yintercept = Gast_palm_pyr_model_Control$coefficients[2], color = "blue", linetype = "longdash") +
# geom_hline(yintercept = Gast_palm_pyr_model_Smoke$coefficients[2], color = "red", linetype = "longdash")


# for(z in names(Gast_palm_Gast_pyr_kinetic_values[1:2])){
Gast_palm_pyr_vmax <- ggplot(data = Gast_palm_pyr_kinetic_values, aes(x = Smoke, y = Vmax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 60)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(10, 60, 10)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 40, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1)


Gast_palm_pyr_km <- ggplot(data = subset(Gast_palm_pyr_kinetic_values), aes(x = Smoke, y = Km, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  ylab(expression(bold(K[m]~(mM)))) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.02, 0.2, 0.02)) +
  labs(y=expression(bold(atop('[Pyruvate]','(mM)')))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .16, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1)

 
# ggstatsplot::ggbetweenstats(data = subset(Gast_palm_pyr_kinetic_values), x = Smoke, y = Km, outlier.tagging = T)
# ggstatsplot::ggbetweenstats(data = subset(Gast_palm_pyr_kinetic_values), x = Smoke, y = Vmax, outlier.tagging = T)
# 
# ggstatsplot::ggbetweenstats(data = Gast_pyr_kinetic_values, x = Smoke, y = Km, outlier.tagging = T)
# ggstatsplot::ggbetweenstats(data = Gast_pyr_kinetic_values, x = Smoke, y = Vmax, outlier.tagging = T)
# 
# ggstatsplot::ggbetweenstats(data = subset(Gast_palm_kinetic_values), x = Smoke, y = Km, outlier.tagging = T)
# ggstatsplot::ggbetweenstats(data = subset(Gast_palm_kinetic_values), x = Smoke, y = Vmax, outlier.tagging = T)

# Bar Graphs ----
## Pyruvate data ----
Gast_pyr_vmax <- ggplot(data = Gast_pyr_kinetic_values, aes(x = Smoke, y = Vmax, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 90)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 90, 10)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 72, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0)


Gast_pyr_km <- ggplot(data = Gast_pyr_kinetic_values, aes(x = Smoke, y = Km, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  ylab(expression(bold(atop('[Pyruvate]', '(mM)')))) +
  ggtitle(expression(bold(K[m]))) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.9)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.9, 0.1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .75, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, tip.length = 0)

## Palmitate -----
Gast_palm_vmax <- ggplot(data = subset(Gast_palm_kinetic_values), aes(x = Smoke, y = Vmax, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  labs(y=expression(bold(V[max]~Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 25)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 25, 5)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 22, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, paired = T, tip.length = 0)


Gast_palm_km <- ggplot(data = subset(Gast_palm_kinetic_values, Km <0.07), aes(x = Smoke, y = Km, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  ylab(expression(bold(K[m]~(mM)))) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.007)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.00, 0.007, 0.001)) +
  labs(y=expression(bold(atop('[PC]', (mM))))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .0055, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0)


## Pyruvate in Palmitate ----
Gast_palm_pyr_vmax <- ggplot(data = subset(Gast_palm_pyr_kinetic_values, Vmax < 48), aes(x = Smoke, y = Vmax, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 60)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 60, 10)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = 51, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0)


Gast_palm_pyr_km <- ggplot(data = subset(Gast_palm_pyr_kinetic_values, Km < 0.26), aes(x = Smoke, y = Km, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  ylab(expression(bold(K[m]~(mM)))) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 0.25)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.0, 0.25, 0.05)) +
  labs(y=expression(bold(atop('[Pyruvate]','(mM)')))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(label.y = .19, aes(label = paste0("p = ", ..p.format..)), comparisons = smoke_comparisons, bracket.size = 1, method = "t.test", tip.length = 0)

pyr_ordered <- Gast_pyr_kinetic_values |> 
  arrange(Smoke)

palm_ordered <- Gast_palm_kinetic_values |> 
  arrange(Smoke)

palm_pyr_ordered <- Gast_palm_pyr_kinetic_values |> 
  arrange(Smoke)


Gast_percent <- data.frame("Km" = palm_pyr_ordered$Km/pyr_ordered$Km,
                               "Vmax" = palm_pyr_ordered$Vmax/pyr_ordered$Vmax,
                               "Subject" = pyr_ordered$Subject,
                               "Smoke" = pyr_ordered$Smoke,
                               "Tissue" = pyr_ordered$Tissue)

Gast_delta <- data.frame("Km" = palm_pyr_ordered$Km-pyr_ordered$Km,
                               "Vmax" = palm_pyr_ordered$Vmax-pyr_ordered$Vmax,
                               "Subject" = pyr_ordered$Subject,
                               "Smoke" = pyr_ordered$Smoke,
                               "Tissue" = pyr_ordered$Tissue)

Gast_contributions <- data.frame("Pyruvate" = palm_pyr_ordered$Vmax/palm_pyr_ordered$Vmax,
                                "PC" = palm_ordered$Vmax/palm_pyr_ordered$Vmax,
                                "Subject" = pyr_ordered$Subject,
                                "Smoke" = pyr_ordered$Smoke,
                                "Tissue" = pyr_ordered$Tissue)
