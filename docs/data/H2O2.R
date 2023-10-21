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

# inhib_kinetics <- nls(Rate ~ -(Vmax *
#                                 Concentration / (Km + Concentration*(1+Concentration/Ks))),
#                       data= subset(h2o2_mm, Group %in% "Smoke"), start=list(Km=0.025,
#                                               Vmax=11, Ks=1))
# summary(inhib_kinetics)

h2o2_kinetic_values <- data.frame()

for(y in unique(h2o2_mm$Group)){
  for(x in unique(h2o2_mm$Tissue)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      h2o2_model <- drm(Rate ~ Concentration, data = subset(subset(h2o2_mm[-c(40, 53, 54),], Group %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower Limit", "Upper Limit", "IC50")))
      
      h2o2_model$Tissue <- paste(x)
      
      h2o2_model$Smoke <- paste(y)
      
      coefs <- data.frame("Imax" = c(h2o2_model$coefficients[1]), "Vmax" = c(h2o2_model$coefficients[2]), "IC50" = c(h2o2_model$coefficients[3]))
      
      coefs$Tissue <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      h2o2_kinetic_values <- rbind(h2o2_kinetic_values, coefs)
      
      mml <- data.frame(S = seq(0, max(h2o2_mm$Concentration), length.out = 100))
      
      mml$v <- predict(h2o2_model, newdata =  mml)
      
      # modelselect <- mselect(h2o2_model, fctList = list(MM.2(), MM.3(), LL.3(), LL.4(), W2.3(), W2.4(), W1.3(), W1.4()))
      # fit <- modelFit(h2o2_model)
      # 
      # assign(paste(paste(x), "h2o2_summary", paste(y), sep = "_"), summary(h2o2_model))
      # assign(paste(paste(x), "h2o2_fit", paste(y), sep = "_"), fit)
      # assign(paste(paste(x), "h2o2_select", paste(y), sep = "_"), modelselect)
      assign(paste(paste(x), "h2o2_residuals", sep = "_"), residuals(h2o2_model))
      assign(paste(paste(x), "h2o2_model", paste(y), sep = "_"), h2o2_model)
      assign(paste(paste(x), "h2o2_mml", paste(y), sep = "_"), mml)
      
      
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  }
}


# Gastrocnemius_h2o2_mml_Control$Smoke <- "Control"
# 
# Gastrocnemius_h2o2_mml_Smoke$Smoke <- "Smoke"
# 
# Gastrocnemius_h2o2_mml_Control$Tissue <- "Gastrocnemius"
# Gastrocnemius_h2o2_mml_Smoke$Tissue <- "Gastrocnemius"
# 
# Soleus_h2o2_mml_Control$Smoke <- "Control"
# 
# Soleus_h2o2_mml_Smoke$Smoke <- "Smoke"
# 
# Soleus_h2o2_mml_Control$Tissue <- "Soleus"
# Soleus_h2o2_mml_Smoke$Tissue <- "Soleus"

# h2o2_mml <- rbind(Gastrocnemius_h2o2_mml_Control, Gastrocnemius_h2o2_mml_Smoke, Soleus_h2o2_mml_Smoke)#Soleus_h2o2_mml_Control, 



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

# individ_h2o2_kinetics$Sex <- h2o2_data$Sex

individ_h2o2_kinetics$Condition <- paste(individ_h2o2_kinetics$Tissue, individ_h2o2_kinetics$Smoke)
individ_h2o2_kinetics$Smoke <- factor(individ_h2o2_kinetics$Smoke, levels = unique(individ_h2o2_kinetics$Smoke))
individ_h2o2_kinetics$Smoke <- factor(individ_h2o2_kinetics$Smoke, levels = unique(individ_h2o2_kinetics$Smoke))

individ_h2o2_mml$Condition <- paste(individ_h2o2_mml$Tissue, individ_h2o2_mml$Smoke)
# h2o2_mml$Condition <- paste(h2o2_mml$Tissue, h2o2_mml$Smoke)
# h2o2_mm$Condition <- paste(h2o2_mm$Tissue, h2o2_mm$Smoke)

individ_h2o2_kinetics$Smoke <- factor(individ_h2o2_kinetics$Smoke, levels = unique(individ_h2o2_kinetics$Smoke))
individ_h2o2_kinetics$Tissue <- factor(individ_h2o2_kinetics$Tissue, levels = unique(individ_h2o2_kinetics$Tissue))

# h2o2_imax <- ggplot(data = individ_h2o2_kinetics, aes(x = Condition, y = Imax, fill = Condition)) +
#   # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   stat_summary(geom = "bar", position = position_dodge(0.9), fill = "white", size = 1) +
#   theme_prism() +
#   scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   # geom_bar_pattern(aes(x = Condition, y = Imax, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
#   # theme_prism() +
#   # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
#   # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   coord_cartesian(ylim = c(0, 0.25)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.25, 0.05)) +
#   labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
#   #ggtitle(expression(bold(I[max]))) +
#   theme(axis.title.x = element_blank()) + 
#   new_scale_fill() +
#   geom_beeswarm(data = individ_h2o2_kinetics, aes(x = Condition, y = Imax, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
#   scale_fill_manual(values = c("white", "white", "white", "white")) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank())
# 
# anova(art(Imax ~ Smoke * Tissue, data = individ_h2o2_kinetics))
# 
# h2o2_km <- ggplot(data = individ_h2o2_kinetics, aes(x = Condition, y = IC50, fill = Condition)) +
#   # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   # stat_summary(geom = "bar", position = position_dodge(0.9), fill = "white", size = 1) +
#   # theme_prism() +
#   # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   geom_bar_pattern(aes(x = Condition, y = IC50, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
#   theme_prism() +
#   scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
#   scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   coord_cartesian(ylim = c(0, 120)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
#   labs(y= expression(bold('[ADP] (μM)'))) +
#   #ggtitle(expression(bold(IC['50']))) +
#   theme(axis.title.x = element_blank()) + 
#   new_scale_fill() +
#   # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   geom_beeswarm(data = individ_h2o2_kinetics, aes(x = Condition, y = IC50, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
#   scale_fill_manual(values = c("white", "white", "white", "white")) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank())
# 
# anova(art(IC50 ~ Smoke * Tissue, data = individ_h2o2_kinetics))
# 
# h2o2_vmax <- ggplot(data = individ_h2o2_kinetics, aes(x = Condition, y = C, fill = Condition)) +
#   # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   # stat_summary(geom = "bar", position = position_dodge(0.9), fill = "white", size = 1) +
#   # theme_prism() +
#   # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   geom_bar_pattern(aes(x = Condition, y = C, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
#   theme_prism() +
#   scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
#   scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.2)) +
#   labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
#   #ggtitle(expression(bold(Initial~Rate))) +
#   theme(axis.title.x = element_blank()) + 
#   new_scale_fill() +
#   # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   geom_beeswarm(data = individ_h2o2_kinetics, aes(x = Condition, y = C, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
#   scale_fill_manual(values = c("white", "white", "white", "white")) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank())
# 
# anova(art(C ~ Smoke * Tissue, data = individ_h2o2_kinetics))
# 
# h2o2_percent <- ggplot(data = individ_h2o2_kinetics, aes(x = Condition, y = Inhibition, fill = Condition)) +
#   # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
#   # stat_summary(geom = "bar", position = position_dodge(0.9), fill = "white", size = 1) +
#   # theme_prism() +
#   geom_bar_pattern(aes(x = Condition, y = Inhibition, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
#   theme_prism() +
#   scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
#   scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   coord_cartesian(ylim = c(-1, 0)) +
#   scale_y_continuous(expand = c(0,0), breaks = seq(-1, 0, 0.1), labels = scales::percent_format()) +
#   labs(y=expression(bold(Inhibition~('%'~V[max])))) +
#   #ggtitle(expression(bold(H['2']*O['2']~Inhibition))) +
#   theme(axis.title.x = element_blank()) +
#   scale_x_discrete(position = "top") + 
#   new_scale_fill() +
#   # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   geom_beeswarm(data = individ_h2o2_kinetics, aes(x = Condition, y = Inhibition, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
#   scale_fill_manual(values = c("white", "white", "white", "white")) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank())
# 
# anova(art(Inhibition ~ Smoke * Tissue, data = individ_h2o2_kinetics))

# h2o2_mml$Condition <- paste(h2o2_mml$Tissue, h2o2_mml$Smoke)

h2o2_errors <- h2o2_mm[-c(1,3,5,6,7,8,9)] |> 
  group_by(Tissue, Group, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

h2o2_errors$min <- h2o2_errors$Rate_mean - (h2o2_errors$Group == "Smoke")*h2o2_errors$Rate_se
h2o2_errors$max <- h2o2_errors$Rate_mean + (h2o2_errors$Group == "Control")*h2o2_errors$Rate_se

# h2o2_kinetics <- ggplot(data = h2o2_mm, aes(x = Concentration, y = Rate, fill = Condition)) +
#   theme_prism() +
#   #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
#   geom_errorbar(data = h2o2_errors, aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
#   coord_cartesian(ylim = c(0, .8), xlim = c(-100, 5500), expand = F) +
#   scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.8, 0.1)) +
#   labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
#   geom_line(data = individ_h2o2_mml, aes(x = S, y = v, linetype = Condition), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
#   scale_linetype_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c(1,3,1,3)) +
#   scale_fill_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("white", "grey50", "grey75", "grey25")) +
#   geom_point(data = h2o2_mm, aes(x = Concentration, y = Rate, fill = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)

h2o2_long$Condition <- paste(h2o2_long$Tissue, h2o2_long$Group)

# mainplot_h2o2 <- ggplot(data = subset(h2o2_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Condition)) +
#   # stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
#   #geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   # theme_prism() +
#   # scale_fill_manual(values = c("white", "grey33", "grey66", "black")) +
#   geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
#   theme_prism() +
#   scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
#   scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   coord_cartesian(ylim = c(0, 1), clip = "off") +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.1)) +
#   labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
#   theme(axis.title.x = element_blank(),
#         legend.position = "bottom") +
#   geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
#   new_scale_fill() +
#   # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
#   #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
#   #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   geom_beeswarm(data = subset(h2o2_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Condition), 
#                 dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
#   scale_fill_manual(values = c("white", "white", "white", "white"))


# Import Resp data ----
resp_data <- read_excel('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\H2O2\\H2O2.xlsx', sheet = "Final Resp Data")

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

adp_kinetic_values <- data.frame()

for(y in unique(adp_mm$Smoke)){
  for(x in unique( adp_mm$Tissue)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      adp_model <- drm(Rate ~ Concentration, data = subset(subset(adp_mm, Smoke %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower","Vmax", "Km")))
      
      adp_model$Tissue <- paste(x)
      
      adp_model$Smoke <- paste(y)
      
      coefs <- data.frame("Vmax" = c(adp_model$coefficients[2]), "Km" = c(adp_model$coefficients[3]))
      
      coefs$Tissue <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      adp_kinetic_values <- rbind(adp_kinetic_values, coefs)
      
      mml <- data.frame(S = seq(0, max(adp_mm$Concentration), length.out = 100))
      
      mml$v <- predict(adp_model, newdata =  mml)
      
      modelselect <- mselect(adp_model, fctList = list(MM.2(), MM.3()))
      fit <- modelFit(adp_model)
      
      assign(paste(paste(x), "adp_residuals", paste(y), sep = "_"), residuals(adp_model))
      assign(paste(paste(x), "adp_summary", paste(y), sep = "_"), summary(adp_model))
      assign(paste(paste(x), "adp_fit", paste(y), sep = "_"), fit)
      assign(paste(paste(x), "adp_select", paste(y), sep = "_"), modelselect)
      
      assign(paste(paste(x), "adp_model", paste(y), sep = "_"), adp_model)
      assign(paste(paste(x), "adp_mml", paste(y), sep = "_"), mml)
      
      
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  }
}

Gastrocnemius_adp_mml_Control$Smoke <- "Control"

Gastrocnemius_adp_mml_Smoke$Smoke <- "Smoke"

Gastrocnemius_adp_mml_Control$Tissue <- "Gastrocnemius"
Gastrocnemius_adp_mml_Smoke$Tissue <- "Gastrocnemius"

Soleus_adp_mml_Control$Smoke <- "Control"

Soleus_adp_mml_Smoke$Smoke <- "Smoke"

Soleus_adp_mml_Control$Tissue <- "Soleus"
Soleus_adp_mml_Smoke$Tissue <- "Soleus"

adp_mml <- rbind(Gastrocnemius_adp_mml_Control, Gastrocnemius_adp_mml_Smoke, Soleus_adp_mml_Smoke, Soleus_adp_mml_Control)


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

# individ_adp_kinetics$Sex <- resp_data$Sex
individ_adp_kinetics$Condition <- paste(individ_adp_kinetics$Tissue, individ_adp_kinetics$Smoke)
individ_adp_kinetics$Smoke <- factor(individ_adp_kinetics$Smoke, levels = unique(individ_adp_kinetics$Smoke))
individ_adp_kinetics$Tissue <- factor(individ_adp_kinetics$Tissue, levels = unique(individ_adp_kinetics$Tissue))


individ_adp_mml$Condition <- paste(individ_adp_mml$Tissue, individ_adp_mml$Smoke)
adp_mml$Condition <- paste(adp_mml$Tissue, adp_mml$Smoke)
adp_mm$Condition <- paste(adp_mm$Tissue, adp_mm$Smoke)

resp_km <- ggplot(data = individ_adp_kinetics[-c(14,23, 17, 18, 13),], aes(x = Condition, y = Km, fill = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.99), width = 0.3, size = 1, na.rm = TRUE) +
  # stat_summary(geom = "bar", position = position_dodge(0.9), fill = "white", size  = 1) +
  # theme_prism() +
  # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  geom_bar_pattern(aes(x = Condition, y = Km, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  theme_prism() +
  scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 450, 50)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  #ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  new_scale_fill() +
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_adp_kinetics[-c(14,23, 17, 18, 13),], aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

resp_vmax <- ggplot(data = individ_adp_kinetics[], aes(x = Condition, y = Vmax, fill = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.99), width = 0.3, size = 1, na.rm = TRUE) +
  # stat_summary(geom = "bar", position = position_dodge(0.9), fill = "white", size = 1) +
  # theme_prism() +
  # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  ggpattern::geom_bar_pattern(aes(x = Condition, y = Vmax, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  theme_prism() +
  scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 140)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 140, 10)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  #ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  new_scale_fill() +
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_adp_kinetics, aes(x = Condition, y = Vmax, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  scale_fill_manual(values = c("white", "white", "white", "white")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

adp_errors <- adp_mm[-c(1,3,5,6,7,8,9)] |> 
  group_by(Tissue, Smoke, Concentration) |> 
  summarise_each(funs(mean, se=sd(.)/sqrt(n())))

adp_errors$min <- adp_errors$Rate_mean - (adp_errors$Smoke == "Smoke")*adp_errors$Rate_se
adp_errors$max <- adp_errors$Rate_mean + (adp_errors$Smoke == "Control")*adp_errors$Rate_se

resp_kinetics <- ggplot(data = adp_mm, aes(x = Concentration, y = Rate, fill = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = adp_errors, aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
  coord_cartesian(ylim = c(0, 100), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,10)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  geom_line(data = individ_adp_mml, aes(x = S, y = v, linetype = Condition), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c(1,3,1,3)) +
  scale_fill_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("white", "grey50", "grey75", "grey25")) +
  geom_point(data = adp_mm, aes(x = Concentration, y = Rate, fill = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)

resp_long$Condition <- paste(resp_long$Tissue, resp_long$Smoke)

mainplot_resp <- ggplot(data = subset(resp_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  #geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  # theme_prism() +
  # ggpattern::geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  # scale_fill_manual(values = c("white", "grey33", "grey66", "black")) +
  coord_cartesian(ylim = c(0, 140), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 140, 20)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(resp_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T)


# H2O2 Per O2

h2o2_per <- data.frame(resp_data[c(1:6, 16, 18)])

h2o2_per <- cbind(h2o2_per, (h2o2_data[8:14]/resp_data[8:14]))

h2o2_per_long <- pivot_longer(h2o2_per, 9:15, names_to = "State", values_to = "Rate")

h2o2_per_long$State <- factor(h2o2_per_long$State, levels = unique(h2o2_per_long$State))

mainplot_per <- ggplot(data = subset(h2o2_per_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  #geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  # theme_prism() +
  # ggpattern::geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  # scale_fill_manual(values = c("white", "grey33", "grey66", "black")) +
  coord_cartesian(ylim = c(0, 0.04), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.04, 0.01)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']/bolditalic(J)*O['2']))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(h2o2_per_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T)


# ggplot(data = subset(resp_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Condition)) +
#   stat_summary(geom = "crossbar", position = position_dodge(0.99), fun = "mean", fill = "white", show.legend = FALSE, fatten = 2.5) +
#   #geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   theme_prism() +
#   # scale_color_manual(values = c("white", "grey33", "grey67", "black")) +
#   coord_cartesian(ylim = c(0, 100), clip = "off") +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 100, 10)) +
#   labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
#   theme(axis.title.x = element_blank(),
#         legend.position = "bottom") +
#   geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
#   new_scale_fill() +
#   # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
#   #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
#   #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   geom_beeswarm(data = subset(resp_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, fill = Condition), 
#                 dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, groupOnX = T) +
#   scale_fill_manual(values = c("white", "grey33", "grey67", "black"))



## Plots -------------------
identify_outliers(data = subset(individ_h2o2_kinetics, Condition %in% "Soleus Control"), Imax)
imax_anova <- anova(art(data = individ_h2o2_kinetics[-c(35, 36),], Imax ~ Smoke * Tissue))
effectsize::eta_squared(imax_anova)

h2o2_imax <- ggplot(data = individ_h2o2_kinetics[-c(35, 36),], aes(x = Tissue, y = Imax, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  # geom_bar_pattern(aes(x = Condition, y = Imax, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 0.25), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.25, 0.05)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  #ggtitle(expression(bold(I[max]))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +  
  geom_beeswarm(data = individ_h2o2_kinetics[-c(35, 36),], aes(x = Tissue, y = Imax, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 0.249, label = 'Smoke: p = 0.772', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.239, label = 'Tissue: p = 0.029', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.229, label = 'Interaction: p = 0.665', hjust = 0, fontface = 2)

identify_outliers(data = subset(individ_h2o2_kinetics, Condition %in% "Gastrocnemius Control"), IC50)
IC50_anova <- anova(art(data = individ_h2o2_kinetics[-c(34, 36),], IC50 ~ Smoke * Tissue))
effectsize::eta_squared(IC50_anova)

h2o2_km <- ggplot(data = individ_h2o2_kinetics[-c(34, 36)], aes(x = Tissue, y = IC50, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  # geom_bar_pattern(aes(x = Condition, y = IC50, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 120), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  #ggtitle(expression(bold(IC['50']))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +  
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_h2o2_kinetics[-c(34, 36)], aes(x = Tissue, y = IC50, color = Condition), dodge.width = .95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 119, label = 'Smoke: p = 0.075', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 114, label = 'Tissue: p = 0.014', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 109, label = 'Interaction: p = 0.144', hjust = 0, fontface = 2)


identify_outliers(data = subset(individ_h2o2_kinetics, Condition %in% "Soleus Control"), C)
C_anova <- anova(art(data = individ_h2o2_kinetics, C ~ Smoke * Tissue))
effectsize::eta_squared(C_anova)

h2o2_vmax <- ggplot(data = individ_h2o2_kinetics, aes(x = Tissue, y = C, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  # geom_bar_pattern(aes(x = Condition, y = C, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.2)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  #ggtitle(expression(bold(Initial~Rate))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +  
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_h2o2_kinetics, aes(x = Tissue, y = C, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 0.99, label = 'Smoke: p = 0.489', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.94, label = 'Tissue: p = 0.100', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.88, label = 'Interaction: p = 0.487', hjust = 0, fontface = 2)


identify_outliers(data = subset(individ_h2o2_kinetics, Condition %in% "Soleus Control"), Inhibition)
Inhibition_anova <- anova(art(data = individ_h2o2_kinetics, Inhibition ~ Smoke * Tissue))
effectsize::eta_squared(Inhibition_anova)

h2o2_percent <- ggplot(data = individ_h2o2_kinetics, aes(x = Tissue, y = Inhibition, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  # geom_bar_pattern(aes(x = Condition, y = Inhibition, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(-1, 0), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1, 0, 0.1), labels = scales::percent_format()) +
  labs(y=expression(bold(Inhibition~('%'~V[max])))) +
  #ggtitle(expression(bold(H['2']*O['2']~Inhibition))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  scale_x_discrete(position = "top") + 
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_h2o2_kinetics, aes(x = Tissue, y = Inhibition, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = -0.86, label = 'Smoke: p = 0.433', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = -0.92, label = 'Tissue: p = 0.257', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = -0.98, label = 'Interaction: p = 0.302', hjust = 0, fontface = 2)

h2o2_errors$Condition <- paste(h2o2_errors$Tissue, h2o2_errors$Group)

h2o2_kinetics <- ggplot(data = h2o2_mm, aes(x = Concentration, y = Rate, color = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = h2o2_errors, aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
  coord_cartesian(ylim = c(0, .8), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.8, 0.1)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  geom_line(data = individ_h2o2_mml, aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  geom_point(data = h2o2_mm, aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)

h2o2_kinetics_Soleus <- ggplot(data = subset(h2o2_mm, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, color = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = subset(h2o2_errors, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Smoke), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
  coord_cartesian(ylim = c(0, .8), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.8, 0.1)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  geom_line(data = subset(individ_h2o2_mml, Tissue %in% "Soleus"), aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  #ggtitle("Soleus") +
  geom_point(data = subset(h2o2_mm, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)


mainplot_h2o2 <- ggplot(data = subset(h2o2_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  # geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.1)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(h2o2_long[-c(273, 269:272, 281, 285:288),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  scale_fill_manual(values = c("white", "white")) +
  #ggtitle(expression(bold(Gastrocnemius~H['2']*O['2']*~Emission))) +
  #GM
  annotate('text', x = .51, y = 0.99, label = 'Smoke: p = 0.544', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.94, label = 'Tissue: p = 0.014', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.89, label = 'Interaction: p = 0.209', hjust = 0, fontface = 2) +
  #GMS
  annotate('text', x = 1.51, y = 0.99, label = 'Smoke: p = 0.202', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 0.94, label = 'Tissue: p = 0.100', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 0.89, label = 'Interaction: p = 0.251', hjust = 0, fontface = 2) +
  #ADP25
  annotate('text', x = 2.51, y = 0.99, label = 'Smoke: p = 0.146', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 0.94, label = 'Tissue: p = 0.020', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 0.89, label = 'Interaction: p = 0.159', hjust = 0, fontface = 2) +
  #ADP50
  annotate('text', x = 3.51, y = 0.99, label = 'Smoke: p = 0.771', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 0.94, label = 'Tissue: p = 0.136', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 0.89, label = 'Interaction: p = 0.305', hjust = 0, fontface = 2) +
  #ADP100
  annotate('text', x = 4.51, y = 0.99, label = 'Smoke: p = 0.757', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 0.94, label = 'Tissue: p = 0.070', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 0.89, label = 'Interaction: p = 0.269', hjust = 0, fontface = 2) +
  #ADP250
  annotate('text', x = 5.51, y = 0.99, label = 'Smoke: p = 0.926', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 0.94, label = 'Tissue: p = 0.204', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 0.89, label = 'Interaction: p = 0.246', hjust = 0, fontface = 2) +
  #ADP5000
  annotate('text', x = 6.51, y = 0.99, label = 'Smoke: p = 0.694', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 0.94, label = 'Tissue: p = 0.132', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 0.89, label = 'Interaction: p = 0.568', hjust = 0, fontface = 2)

mainplot_h2o2_Soleus <- ggplot(data = subset(subset(h2o2_long[-c(273, 269:272, 281, 285:288),], Tissue %in% "Soleus"), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  # geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.1)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~(pmol[H['2']*O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  new_scale_fill() +
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(subset(h2o2_long, Tissue %in% "Soleus"), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  scale_fill_manual(values = c("white", "white"))
  #ggtitle(expression(bold(Soleus~H['2']*O['2']*~Emission)))

km_anova <- anova(art(data = individ_adp_kinetics, Km ~ Smoke * Tissue))
effectsize::eta_squared(km_anova)

resp_km <- ggplot(data = individ_adp_kinetics[-c(14,23, 17, 18, 13),], aes(x = Tissue, y = Km, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.99), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size  = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  # geom_bar_pattern(aes(x = Condition, y = Km, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  coord_cartesian(ylim = c(0, 450)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 450, 50)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  #ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank()) +
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_adp_kinetics[-c(14,23, 17, 18, 13),], aes(x = Tissue, y = Km, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  scale_fill_manual(values = c("white", "white")) +
  theme(axis.title.x = element_blank())

vmax_anova <- anova(art(data = individ_adp_kinetics, Vmax ~ Smoke * Tissue))
effectsize::eta_squared(vmax_anova)

resp_vmax <- ggplot(data = individ_adp_kinetics[], aes(x = Tissue, y = Vmax, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.99), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.95), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  # ggpattern::geom_bar_pattern(aes(x = Condition, y = Vmax, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 140), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 140, 10)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  #ggtitle(expression(bold(V[max]))) +
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_adp_kinetics, aes(x = Tissue, y = Vmax, color = Condition), dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  theme(axis.title.x = element_blank()) +
  #stats
  annotate('text', x = .51, y = 139, label = 'Smoke: p = 0.150', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 134, label = 'Tissue: p = 0.096', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 129, label = 'Interaction: p = 0.795', hjust = 0, fontface = 2)


resp_kinetics <- ggplot(data = adp_mm, aes(x = Concentration, y = Rate, fill = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = adp_errors, aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Condition), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
  coord_cartesian(ylim = c(0, 100), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,10)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  geom_line(data = individ_adp_mml, aes(x = S, y = v, linetype = Condition), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c(1,3,1,3)) +
  scale_fill_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("white", "grey50", "grey75", "grey25")) +
  geom_point(data = adp_mm, aes(x = Concentration, y = Rate, fill = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)


mainplot_resp <- ggplot(data = subset(resp_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  # geom_point(position = position_dodge(0.95), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  # ggpattern::geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 140), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 140, 20)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(resp_long, State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  #ggtitle("Gastrocnemius Respiration") +
  #GM
  annotate('text', x = .51, y = 135, label = 'Smoke: p = 0.946', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 125, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 115, label = 'Interaction: p = 0.104', hjust = 0, fontface = 2) +
  #GMS
  annotate('text', x = 1.51, y = 135, label = 'Smoke: p = 0.679', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 125, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 115, label = 'Interaction: p = 0.147', hjust = 0, fontface = 2) +
  #ADP25
  annotate('text', x = 2.51, y = 135, label = 'Smoke: p = 0.153', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 125, label = 'Tissue: p = 0.032', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 115, label = 'Interaction: p = 0.772', hjust = 0, fontface = 2) +
  #ADP50
  annotate('text', x = 3.51, y = 135, label = 'Smoke: p = 0.203', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 125, label = 'Tissue: p = 0.105', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 115, label = 'Interaction: p = 0.784', hjust = 0, fontface = 2) +
  #ADP100
  annotate('text', x = 4.51, y = 135, label = 'Smoke: p = 0.099', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 125, label = 'Tissue: p = 0.082', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 115, label = 'Interaction: p = 0.969', hjust = 0, fontface = 2) +
  #ADP250
  annotate('text', x = 5.51, y = 135, label = 'Smoke: p = 0.032', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 125, label = 'Tissue: p = 0.093', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 115, label = 'Interaction: p = 0.956', hjust = 0, fontface = 2) +
  #ADP5000
  annotate('text', x = 6.51, y = 135, label = 'Smoke: p = 0.018', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 125, label = 'Tissue: p = 0.051', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 115, label = 'Interaction: p = 0.746', hjust = 0, fontface = 2)

mainplot_resp_Soleus <- ggplot(data = subset(subset(resp_long, Tissue %in% 'Soleus'), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.95), size = 1, fill = "white") +
  # geom_point(position = position_dodge(0.95), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  # ggpattern::geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  # scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  scale_fill_manual(values = c("white", "grey33")) +
  coord_cartesian(ylim = c(0, 140), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 140, 20)) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  new_scale_fill() +
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(subset(resp_long, Tissue %in% 'Soleus'), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  scale_fill_manual(values = c("white", "white"))

h2o2_per_long$Condition <- factor(h2o2_per_long$Condition, levels = unique(h2o2_per_long$Condition))

mainplot_per <- ggplot(data = subset(h2o2_per_long[-c(29, 30, 176, 246:249),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
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
  geom_beeswarm(data = subset(h2o2_per_long[-c(29, 30, 176, 246:249),], State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  #GM
  annotate('text', x = .51, y = 0.039, label = 'Smoke: p = 0.913', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.036, label = 'Tissue: p = 0.103', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.033, label = 'Interaction: p = 0.587', hjust = 0, fontface = 2) +
  #GMS
  annotate('text', x = 1.51, y = 0.039, label = 'Smoke: p = 0.569', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 0.036, label = 'Tissue: p = 0.028', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 0.033, label = 'Interaction: p = 0.303', hjust = 0, fontface = 2) +
  #ADP25
  annotate('text', x = 2.51, y = 0.039, label = 'Smoke: p = 0.413', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 0.036, label = 'Tissue: p = 0.249', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 0.033, label = 'Interaction: p = 0.797', hjust = 0, fontface = 2) +
  #ADP50
  annotate('text', x = 3.51, y = 0.039, label = 'Smoke: p = 0.140', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 0.036, label = 'Tissue: p = 0.346', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 0.033, label = 'Interaction: p = 0.601', hjust = 0, fontface = 2) +
  #ADP100
  annotate('text', x = 4.51, y = 0.039, label = 'Smoke: p = 0.460', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 0.036, label = 'Tissue: p = 0.180', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 0.033, label = 'Interaction: p = 0.203', hjust = 0, fontface = 2) +
  #ADP250
  annotate('text', x = 5.51, y = 0.039, label = 'Smoke: p = 0.257', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 0.036, label = 'Tissue: p = 0.354', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 0.033, label = 'Interaction: p = 0.206', hjust = 0, fontface = 2) +
  #ADP5000
  annotate('text', x = 6.51, y = 0.039, label = 'Smoke: p = 0.230', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 0.036, label = 'Tissue: p = 0.152', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 0.033, label = 'Interaction: p = 0.480', hjust = 0, fontface = 2)

mainplot_per_Soleus <- ggplot(data = subset(subset(h2o2_per_long[-c(29, 30, 176, 246:249),], Tissue %in% "Soleus"), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
  # geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  # ggpattern::geom_bar_pattern(aes(x = State, y = Rate, fill = Condition, pattern = Tissue), position = position_dodge(0.95), stat = "summary", fill = "white", size = 1, pattern_fill = "black", pattern_angle = 135, show.legend = F) +
  # theme_prism() +
  # scale_pattern_manual(breaks = c("Gastrocnemius", "Soleus"), values = c('none', "stripe")) +
  # scale_fill_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("white", "grey33")) +
  scale_fill_manual(values = c("white", "grey33")) +
  coord_cartesian(ylim = c(0, 0.04), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.04, 0.005)) +
  labs(y=expression(bold(bolditalic(J)*H['2']*O['2']/bolditalic(J)*O['2']))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), linetype = "longdash", color = "grey50", size = .6) + 
  new_scale_fill() +
  # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
  #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
  #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = subset(subset(h2o2_per_long[-c(29, 30, 176, 246:249),], Tissue %in% "Soleus"), State %in% c("GM", "GMS", 'ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.95, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  scale_fill_manual(values = c("white", "white"))

# ggsave("Main_Per.png", mainplot_per, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# ggsave("Main_Resp.png", mainplot_resp, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# ggsave("Main_H2O2.png", mainplot_h2o2, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# 
# ggsave("Main_Per_Gastroc.png", mainplot_per_Gastroc, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# ggsave("Main_Resp_Gastroc.png", mainplot_resp_Gastroc, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# ggsave("Main_H2O2_Gastroc.png", mainplot_h2o2_Gastroc, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# 
# ggsave("Main_Per_Soleus.png", mainplot_per_Soleus, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# ggsave("Main_Resp_Soleus.png", mainplot_resp_Soleus, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# ggsave("Main_H2O2_Soleus.png", mainplot_h2o2_Soleus, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# 
# H2O2_plot <- ggarrange(h2o2_vmax, h2o2_km, h2o2_imax, h2o2_percent, common.legend = T, nrow = 2, ncol = 2, legend = "bottom", labels = "AUTO")
# ggsave("H2O2_values.png", H2O2_plot, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 12, height = 8, dpi = 300)
# 
# adp_plot <- ggarrange(resp_km, resp_vmax, common.legend = T, nrow = 1, legend = "bottom", labels = "AUTO")
# ggsave("ADP_values.png", adp_plot, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 12, height = 5, dpi = 300)
# 
# ggsave("Gastroc_kinetics.png", h2o2_kinetics_Gast, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 10, height = 5, dpi = 300)
# ggsave("Soleus_kinetics.png", h2o2_kinetics_Soleus, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 10, height = 5, dpi = 300)
# 
# h2o2_kinetics_plot <- ggarrange(h2o2_kinetics_Gast, h2o2_kinetics_Soleus, common.legend = T, nrow = 2, legend = "bottom", labels = "AUTO")
# ggsave("h2o2_kinetics_plot.png", h2o2_kinetics_plot, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 10, height = 8, dpi = 300)
# 
# all_resp <- ggarrange(mainplot_resp_Gastroc, mainplot_resp_Soleus, nrow = 2, legend = "bottom", labels = "AUTO", common.legend = TRUE)
# ggsave('all_resp.png', all_resp, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 7.5, dpi = 300)
# 
# all_h2o2 <- ggarrange(mainplot_h2o2_Gastroc, mainplot_h2o2_Soleus, nrow = 2, legend = "bottom", labels = "AUTO", common.legend = TRUE)
# ggsave('all_h2o2.png', all_h2o2, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 7.5, dpi = 300)
# 
# all_h2o2_per <- ggarrange(mainplot_per_Gastroc, mainplot_per_Soleus, nrow = 2, legend = "bottom", labels = "AUTO", common.legend = TRUE)
# ggsave('all_h2o2_per.png', all_h2o2_per, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 7.5, dpi = 300)

# Percentages--------------------------

# CAT Kinetics -----

percent_mm <- pivot_longer(h2o2_data[-c(7:14)], cols = 8:13, names_to = "Concentration", values_to = "Rate")

percent_mm <- separate(percent_mm, Concentration, into = c("ADP", "fake"), sep = "_")

percent_mm <- separate(percent_mm[-c(9)], ADP, into = c("ADP", "Concentration"), sep = " ")

percent_mm$Concentration[is.na(percent_mm$Concentration)] <- 0

percent_mm$Concentration <- as.numeric( percent_mm$Concentration)

# inhib_kinetics <- nls(Rate ~ -(Vmax *
#                                 Concentration / (Km + Concentration*(1+Concentration/Ks))),
#                       data= subset(h2o2_mm, Group %in% "Smoke"), start=list(Km=0.025,
#                                               Vmax=11, Ks=1))
# summary(inhib_kinetics)

percent_kinetic_values <- data.frame()

for(y in unique(percent_mm$Group)){
  for(x in unique(percent_mm$Tissue)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      percent_model <- drm(Rate ~ Concentration, data = subset(subset(percent_mm[-c(40, 53, 54),], Group %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower Limit", "Upper Limit", "IC50")))
      
      percent_model$Tissue <- paste(x)
      
      percent_model$Smoke <- paste(y)
      
      coefs <- data.frame("Imax" = c(percent_model$coefficients[1]), "Vmax" = c(percent_model$coefficients[2]), "IC50" = c(percent_model$coefficients[3]))
      
      coefs$Tissue <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      percent_kinetic_values <- rbind(h2o2_kinetic_values, coefs)
      
      mml <- data.frame(S = seq(0, max(percent_mm$Concentration), length.out = 100))
      
      mml$v <- predict(percent_model, newdata =  mml)
      
      # modelselect <- mselect(h2o2_model, fctList = list(MM.2(), MM.3(), LL.3(), LL.4(), W2.3(), W2.4(), W1.3(), W1.4()))
      # fit <- modelFit(h2o2_model)
      # 
      # assign(paste(paste(x), "h2o2_summary", paste(y), sep = "_"), summary(h2o2_model))
      # assign(paste(paste(x), "h2o2_fit", paste(y), sep = "_"), fit)
      # assign(paste(paste(x), "h2o2_select", paste(y), sep = "_"), modelselect)
      assign(paste(paste(x), "percent_residuals", sep = "_"), residuals(percent_model))
      assign(paste(paste(x), "percent_model", paste(y), sep = "_"), percent_model)
      assign(paste(paste(x), "percent_mml", paste(y), sep = "_"), mml)
      
      
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  }
}


# Gastrocnemius_h2o2_mml_Control$Smoke <- "Control"
# 
# Gastrocnemius_h2o2_mml_Smoke$Smoke <- "Smoke"
# 
# Gastrocnemius_h2o2_mml_Control$Tissue <- "Gastrocnemius"
# Gastrocnemius_h2o2_mml_Smoke$Tissue <- "Gastrocnemius"
# 
# Soleus_h2o2_mml_Control$Smoke <- "Control"
# 
# Soleus_h2o2_mml_Smoke$Smoke <- "Smoke"
# 
# Soleus_h2o2_mml_Control$Tissue <- "Soleus"
# Soleus_h2o2_mml_Smoke$Tissue <- "Soleus"

# h2o2_mml <- rbind(Gastrocnemius_h2o2_mml_Control, Gastrocnemius_h2o2_mml_Smoke, Soleus_h2o2_mml_Smoke)#Soleus_h2o2_mml_Control, 



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

# individ_h2o2_kinetics$Sex <- h2o2_data$Sex

individ_percent_kinetics$Condition <- paste(individ_percent_kinetics$Tissue, individ_percent_kinetics$Smoke)

individ_percent_mml$Condition <- paste(individ_percent_mml$Tissue, individ_percent_mml$Smoke)
# percent_mml$Condition <- paste(percent_mml$Tissue, percent_mml$Smoke)
# percent_mm$Condition <- paste(percent_mm$Tissue, percent_mm$Smoke)

individ_percent_kinetics$Smoke <- factor(individ_percent_kinetics$Smoke, levels = unique(individ_percent_kinetics$Smoke))
individ_percent_kinetics$Tissue <- factor(individ_percent_kinetics$Tissue, levels = unique(individ_percent_kinetics$Tissue))

identify_outliers(data = subset(individ_percent_kinetics, Condition %in% "Gastrocnemius Smoke"), Imax)
imax_percent_anova <- anova(art(data = individ_percent_kinetics, IC50 ~ Smoke * Tissue))
effectsize::eta_squared(imax_percent_anova)

percent_imax <- ggplot(data = individ_percent_kinetics, aes(x = Tissue, y = Imax, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, .8), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .8, .1), labels = scales::percent_format()) +
  labs(y = expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +  
  geom_beeswarm(data = individ_percent_kinetics, aes(x = Tissue, y = Imax, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 0.79, label = 'Smoke: p = 0.503', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.75, label = 'Tissue: p = 0.256', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.71, label = 'Interaction: p = 0.334', hjust = 0, fontface = 2)

identify_outliers(data = subset(individ_percent_kinetics, Condition %in% "Gastrocnemius Smoke"), IC50)
km_percent_anova <- anova(art(data = individ_percent_kinetics[-c(34, 36),], IC50 ~ Smoke * Tissue))
effectsize::eta_squared(km_percent_anova)

percent_km <- ggplot(data = individ_percent_kinetics[-c(34, 36),], aes(x = Tissue, y = IC50, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 120), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_percent_kinetics[-c(34, 36),], aes(x = Tissue, y = IC50, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 119, label = 'Smoke: p = 0.064', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 114, label = 'Tissue: p = 0.014', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 109, label = 'Interaction: p = 0.104', hjust = 0, fontface = 2)



# percent_mml$Condition <- paste(percent_mml$Tissue, percent_mml$Smoke)

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

percent_kinetics_Soleus <- ggplot(data = subset(percent_mm, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, color = Condition)) +
  theme_prism() +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
  geom_errorbar(data = subset(percent_errors, Tissue %in% "Soleus"), aes(x = Concentration, y = mean, ymin = min_sol, ymax = max_sol, color = Smoke), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
  coord_cartesian(ylim = c(0, 1.2), xlim = c(-100, 5500), expand = F) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1.2, 0.2)) +
  labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
  geom_line(data = subset(individ_percent_mml, Tissue %in% "Soleus"), aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
  scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  geom_point(data = subset(percent_mm, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)

h2o2_long_percent$Condition <- paste(h2o2_long_percent$Tissue, h2o2_long_percent$Group)

h2o2_long_percent$Condition <- factor(h2o2_long_percent$Condition, levels = unique(h2o2_long_percent$Condition))

mainplot_percent <- ggplot(data = subset(h2o2_long_percent[-c(29, 147),], State %in% c('ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
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
  geom_beeswarm(data = subset(h2o2_long_percent[-c(29, 147),], State %in% c('ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
                dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
  #ADP25
  annotate('text', x = 0.51, y = 1.19, label = 'Smoke: p = 0.434', hjust = 0, fontface = 2) +
  annotate('text', x = 0.51, y = 1.13, label = 'Tissue: p = 0.024', hjust = 0, fontface = 2) +
  annotate('text', x = 0.51, y = 1.07, label = 'Interaction: p = 0.393', hjust = 0, fontface = 2) +
  #ADP50
  annotate('text', x = 1.51, y = 1.19, label = 'Smoke: p = 0.657', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 1.13, label = 'Tissue: p = 0.124', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 1.07, label = 'Interaction: p = 0.172', hjust = 0, fontface = 2) +
  #ADP100
  annotate('text', x = 2.51, y = 1.19, label = 'Smoke: p = 0.855', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 1.13, label = 'Tissue: p = 0.127', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 1.07, label = 'Interaction: p = 0.250', hjust = 0, fontface = 2) +
  #ADP250
  annotate('text', x = 3.51, y = 1.19, label = 'Smoke: p = 0.831', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 1.13, label = 'Tissue: p = 0.243', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 1.07, label = 'Interaction: p = 0.068', hjust = 0, fontface = 2) +
  #ADP5000
  annotate('text', x = 4.51, y = 1.19, label = 'Smoke: p = 0.551', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 1.13, label = 'Tissue: p = 0.119', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 1.07, label = 'Interaction: p = 0.476', hjust = 0, fontface = 2)


# mainplot_percent_Soleus <- ggplot(data = subset(subset(h2o2_long_percent[-c(29, 147),], Tissue %in% "Soleus"), State %in% c('ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition)) +
#   stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, fill = "white") +
#   # geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   theme_prism() +
#   scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   coord_cartesian(ylim = c(0, 1.2), clip = "off") +
#   scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.2, 0.2), labels = scales::percent_format()) +
#   labs(y=expression(bold(bolditalic(J)*H['2']*O['2']~('%'~of~GMS)))) +
#   theme(axis.title.x = element_blank(),
#         legend.position = "bottom") +
#   geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "longdash", color = "grey50", size = .6) + 
#   new_scale_fill() +
#   # geom_point(data = subset(cat_long, State %in% c('GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
#   #                                                 , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Condition), 
#   #            position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   geom_beeswarm(data = subset(subset(h2o2_long_percent[-c(29, 147),], Tissue %in% "Soleus"), State %in% c('ADP 25', 'ADP 50', 'ADP 100', 'ADP 250', 'ADP 5000')), aes(x = State, y = Rate, color = Condition), 
#                 dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, groupOnX = T) +
#   #ggtitle(expression(bold(Soleus~H['2']*O['2']*~Emission~('% of GMS')))) +
#   scale_fill_manual(values = c("white", "white"))

# ggsave("Percent_Gastroc_kinetics.png", percent_kinetics_Gast, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 10, height = 5, dpi = 300)
# ggsave("Percent_Soleus_kinetics.png", percent_kinetics_Soleus, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 10, height = 5, dpi = 300)
# 
# ggsave("Percent_Main_H2O2_Gastroc.png", mainplot_percent_Gastroc, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# ggsave("Percent_Main_H2O2_Soleus.png", mainplot_percent_Soleus, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 5, dpi = 300)
# 
# all_h2o2_percent <- ggarrange(mainplot_percent_Gastroc, mainplot_percent_Soleus, nrow = 2, legend = "bottom", labels = "AUTO", common.legend = TRUE)
# ggsave('all_h2o2_percent.png', all_h2o2_percent, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 15, height = 7.5, dpi = 300)
# 
# percent_plot <- ggarrange(percent_km, percent_imax, common.legend = T, nrow = 1, legend = "bottom", labels = "AUTO")
# ggsave("percent_values.png", percent_plot, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 12, height = 5, dpi = 300)
# 
# h2o2_percent_kinetics_plot <- ggarrange(percent_kinetics_Gast, percent_kinetics_Soleus, common.legend = T, nrow = 2, legend = "bottom", labels = "AUTO")
# ggsave("h2o2_percent_kinetics_plot.png", h2o2_percent_kinetics_plot, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 10, height = 8, dpi = 300)

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
imax_per_anova <- anova(art(data = individ_h2o2_per_kinetics[-c(29, 30),], Imax ~ Smoke * Tissue))
effectsize::eta_squared(imax_per_anova)

h2o2_per_imax <- ggplot(data = individ_h2o2_per_kinetics[-c(29, 30),], aes(x = Tissue, y = Imax, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, .006), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, .006, .001)) +
  labs(y = expression(bold(bolditalic(J)*H['2']*O['2']*'/'*O['2']))) +
  #ggtitle(expression(bold(I[max]))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  geom_beeswarm(data = individ_h2o2_per_kinetics[-c(29, 30),], aes(x = Tissue, y = Imax, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 0.0059, label = 'Smoke: p = 0.006', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.0055, label = 'Tissue: p = 0.970', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.0051, label = 'Interaction: p = 0.177', hjust = 0, fontface = 2)


identify_outliers(data = subset(individ_h2o2_per_kinetics, Condition %in% "Soleus Smoke"), IC50)
km_per_anova <- anova(art(data = individ_h2o2_per_kinetics[-c(29, 30),], IC50 ~ Smoke * Tissue))
effectsize::eta_squared(km_per_anova)

h2o2_per_km <- ggplot(data = individ_h2o2_per_kinetics[-c(29),], aes(x = Tissue, y = IC50, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(0, 120), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
  labs(y= expression(bold('[ADP] (μM)'))) +
  #ggtitle(expression(bold(IC['50']))) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  # geom_point(data = individ_adp_kinetics, aes(x = Condition, y = Km, fill = Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(data = individ_h2o2_per_kinetics[-c(29),], aes(x = Tissue, y = IC50, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 119, label = 'Smoke: p = 0.369', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 114, label = 'Tissue: p = 0.004', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 109, label = 'Interaction: p = 0.230', hjust = 0, fontface = 2)

identify_outliers(data = subset(individ_h2o2_per_kinetics, Condition %in% "Soleus Control"), C)
vmax_per_anova <- anova(art(data = individ_h2o2_per_kinetics[-c(4, 30),], C ~ Smoke * Tissue))
effectsize::eta_squared(vmax_per_anova)

h2o2_per_vmax <- ggplot(data = individ_h2o2_per_kinetics[-c(4, 30),], aes(x = Tissue, y = C, color = Condition)) +
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
  geom_beeswarm(data = individ_h2o2_per_kinetics[-c(4, 30),], aes(x = Tissue, y = C, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = 0.0245, label = 'Smoke: p = 0.333', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.023, label = 'Tissue: p = 0.426', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = 0.0215, label = 'Interaction: p = 0.278', hjust = 0, fontface = 2)

identify_outliers(data = subset(individ_h2o2_per_kinetics, Condition %in% "Soleus Control"), Inhibition)
inhib_per_anova <- anova(art(data = individ_h2o2_per_kinetics[-c(29),], Inhibition ~ Smoke * Tissue))
effectsize::eta_squared(inhib_per_anova)

h2o2_per_inhibition <- ggplot(data = individ_h2o2_per_kinetics[-c(29),], aes(x = Tissue, y = Inhibition, color = Condition)) +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.99), fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
  coord_cartesian(ylim = c(-1.1, 0), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1.1, 0, 0.1), labels = scales::percent_format()) +
  labs(y=expression(bold(Inhibition~('%'~V[max])))) +
  #ggtitle(expression(bold(H['2']*O['2']~Inhibition))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_discrete(position = "top") + 
  geom_beeswarm(data = individ_h2o2_per_kinetics[-c(29),], aes(x = Tissue, y = Inhibition, color = Condition), dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 2) +
  #stats
  annotate('text', x = .51, y = -0.95, label = 'Smoke: p = 0.088', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = -1.02, label = 'Tissue: p = 0.190', hjust = 0, fontface = 2) +
  annotate('text', x = .51, y = -1.09, label = 'Interaction: p = 0.100', hjust = 0, fontface = 2)

all_H2O2_per <- ggarrange(h2o2_per_vmax, h2o2_per_km, h2o2_per_imax, h2o2_per_inhibition, common.legend = T, nrow = 2, ncol = 2, legend = "bottom", labels = "AUTO")
# ggsave("H2O2_per_values.png", all_H2O2_per, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 12, height = 8, dpi = 300)

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

# h2o2_per_kinetics_Soleus <- ggplot(data = subset(h2o2_per_mm, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, color = Condition)) +
#   theme_prism() +
#   #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 100) +
#   geom_errorbar(data = subset(h2o2_per_errors, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate_mean, ymin = min, ymax = max, color = Smoke), inherit.aes = F, size = 1, width = 75, fill = "white", show.legend = F) +
#   coord_cartesian(ylim = c(0, .015), xlim = c(-100, 5500), expand = F) +
#   scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.015, 0.003)) +
#   labs(x = "[ADP] (µM)", y = expression(bold(bolditalic(J)*H['2']*O['2']*'/'*O['2']))) +
#   geom_line(data = subset(individ_h2o2_per_mml, Tissue %in% "Soleus"), aes(x = S, y = v, linetype = Smoke), size = 1, fun.y = "mean", stat = "summary", show.legend = F) +
#   scale_linetype_manual(breaks = c("Control", "Smoke"), values = c(1,3)) +
#   scale_color_manual(breaks = c("Gastrocnemius Control", "Gastrocnemius Smoke","Soleus Control", "Soleus Smoke"), values = c("red", "darkred", "blue", "darkblue")) +
#   #ggtitle("Soleus") +
#   geom_point(data = subset(h2o2_per_mm, Tissue %in% "Soleus"), aes(x = Concentration, y = Rate, color = Condition), stat = "summary", fun.y = "mean", size = 3.5, pch = 21, fill = "white", show.legend = F)

# h2o2_per_kinetics_plot <- ggarrange(h2o2_per_kinetics_Gast, h2o2_per_kinetics_Soleus, common.legend = T, nrow = 2, legend = "bottom", labels = "AUTO")
# ggsave("h2o2_per_kinetics_plot.png", h2o2_per_kinetics_plot, path = "C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\Plots\\H2O2", width = 10, height = 8, dpi = 300)

# Stats ------------------------------------

## H2O2 --------------
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

# resp_descriptives <- aggregate(x = Rate ~ Smoke + Tissue + State, data = cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
# subset(resp_descriptives, State %in% "GMD 5000")

h2o2_long$Group <- factor(h2o2_long$Group, levels = unique(h2o2_long$Group))

for(y in unique(h2o2_long$State)){
  print(subset(h2o2_long[-c(273, 269:272, 281, 285:288),], State %in% paste(y)) |>
    group_by(Tissue) |>
    levene_test(Rate ~ Group))

  print(as.data.table(subset(h2o2_long[-c(273, 269:272, 281, 285:288),], State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
                                         P.value = shapiro.test(Rate)$p.value),
                                      by = .(Tissue, Group)])

}

anova_states_results <- data.frame()
anova_eff_results <- data.frame()

for(y in unique(h2o2_long$State)){
  anova_states <- anova(lm(Rate ~ Group * Tissue, data = subset(na.omit(h2o2_long[-c(273, 269:272, 281, 285:288),]), State %in% paste(y))))
  anova_states$State <- paste(y)

  anova_eff <- effectsize::eta_squared(anova_states)

  print(anova_states)
  print(effectsize::eta_squared(anova_states))

  anova_states_results <- rbind(anova_states_results, anova_states[c(1:3),])
  anova_eff_results <- rbind(anova_eff_results, anova_eff)

}

cbind(anova_states_results, anova_eff_results)

## H2O2 per O2 --------------

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
  print(subset(h2o2_per_long[-c(29, 30, 176, 246:249),], State %in% paste(y)) |>
          group_by(Tissue) |>
          levene_test(Rate ~ Group))
  
  print(as.data.table(subset(h2o2_per_long[-c(29, 30, 176, 246:249),], State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
                                                                                                  P.value = shapiro.test(Rate)$p.value),
                                                                                               by = .(Tissue, Group)])
  
}

anova_per_results <- data.frame()
anova_per_eff_results <- data.frame()

for(y in unique(h2o2_per_long$State)){
  anova_states <- anova(lm(Rate ~ Group * Tissue, data = subset(na.omit(h2o2_per_long[-c(29, 30, 176, 246:249),]), State %in% paste(y))))
  anova_states$State <- paste(y)
  
  anova_eff <- effectsize::eta_squared(anova_states)
  
  print(anova_states)
  print(effectsize::eta_squared(anova_states))
  
  anova_per_results <- rbind(anova_per_results, anova_states[c(1:3),])
  anova_per_eff_results <- rbind(anova_per_eff_results, anova_eff)
  
}

cbind(anova_per_results, anova_per_eff_results)


## H2O2 Percent --------------


outliers_percent <- data.frame()

for(y in unique(h2o2_long_percent$Condition)){
  for(x in names(h2o2_data[c(9:14)])){
    skip_to_next <- FALSE
    
    tryCatch({
      out <- identify_outliers(data = subset(subset(h2o2_long_percent, Condition %in% paste(y)), State %in% paste(x)), Rate)
      
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
  
  print(as.data.table(subset(h2o2_long_percent[-c(29, 147),], State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
                                                                                                  P.value = shapiro.test(Rate)$p.value),
                                                                                               by = .(Tissue, Group)])
  
}

anova_percent_results <- data.frame()
anova_percent_eff_results <- data.frame()

for(y in c("ADP 25", "ADP 50", "ADP 100", "ADP 250", "ADP 5000")){
  anova_states <- anova(lm(Rate ~ Group * Tissue, data = subset(h2o2_long_percent[-c(29, 147),], State %in% paste(y))))
  anova_states$State <- paste(y)
  
  anova_eff <- effectsize::eta_squared(anova_states)
  
  print(anova_states)
  print(effectsize::eta_squared(anova_states))
  
  anova_percent_results <- rbind(anova_percent_results, anova_states[c(1:3),])
  anova_percent_eff_results <- rbind(anova_percent_eff_results, anova_eff)
  
}

cbind(anova_percent_results, anova_percent_eff_results)

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
  print(subset(resp_long, State %in% paste(y)) |>
          group_by(Tissue) |>
          levene_test(Rate ~ Condition))
  
  print(as.data.table(subset(resp_long, State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
                                                                                       P.value = shapiro.test(Rate)$p.value),
                                                                                    by = .(Tissue, Group)])
  
}

anova_percent_results <- data.frame()
anova_percent_eff_results <- data.frame()

for(y in c("GM", "GMS", "ADP 25", "ADP 50", "ADP 100", "ADP 250", "ADP 5000")){
  anova_states <- anova(lm(Rate ~ Group * Tissue, data = subset(resp_long, State %in% paste(y))))
  anova_states$State <- paste(y)
  
  anova_eff <- effectsize::eta_squared(anova_states)
  
  print(anova_states)
  print(effectsize::eta_squared(anova_states))
  
  anova_percent_results <- rbind(anova_percent_results, anova_states[c(1:3),])
  anova_percent_eff_results <- rbind(anova_percent_eff_results, anova_eff)
  
}

cbind(anova_percent_results, anova_percent_eff_results)

