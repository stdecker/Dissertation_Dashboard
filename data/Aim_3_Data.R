library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(drc)
library(EnvStats)
library(ggprism)
library(flextable)
library(dunn.test)
library(ARTool)
library(psych)
library(plotrix)
library(psych)
library(scales)
library(ggnewscale)
library(ggbeeswarm)
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

# Import data ----
cat_data <- read_excel('Analysis File4.xlsx', sheet = "Final Data")

cat_data$Tissue[cat_data$Tissue == "Gastroc"] <- "Gastrocnemius"

# cat_data <- na.omit(cat_data)
cat_data <- cat_data[-c(3,46:299),]

cat_data$Smoke[ cat_data$Smoke == "Con"] <- "Control"

cat_data$Smoke[ cat_data$Smoke == "Smo"] <- "Smoke"

cat_data <- separate( cat_data, Subject, into = "Subject", sep = " ")

cat_data$Smoke <- factor(cat_data$Smoke, levels = unique(cat_data$Smoke))

cat_data$Tissue <- factor(cat_data$Tissue, levels = unique(cat_data$Tissue))

cat_data$Sex <- factor(cat_data$Sex, levels = unique(cat_data$Sex))

cat_long <- pivot_longer( cat_data, cols = 7:23, names_to = "State", values_to = "Rate")

cat_long$State <- factor( cat_long$State, levels = unique(cat_long$State))

cat_data$RCR <- cat_data$GMDS/cat_data$GM

cat_data$Thermodynamic_Coupling <- sqrt(1-(cat_data$`CAT 5.0`/cat_data$`FCCP Peak`))

cat_data$SpecificRate <- cat_data$GM/cat_data$Rot

cat_data$CIIPercent <- cat_data$Rot/cat_data$`FCCP Peak`

cat_data$FCCP_per_CAT <- cat_data$`FCCP Peak`/cat_data$`CAT 5.0`

cat_data$GMDS_per_CAT <- cat_data$GMDS/cat_data$`CAT 5.0`

cat_data$GMDS_per_FCCP <- cat_data$GMDS/cat_data$`FCCP Peak`

cat_data$CAT_cont <- cat_data$GMDS - cat_data$`CAT 5.0`

cat_data$Condition <- paste(cat_data$Tissue, cat_data$Smoke)

# cat_data$CAT_inhib_percent <- -(1-(cat_data$`CAT 5.0`/cat_data$GMDS))
# 
# aggregate(CAT_inhib_percent ~ Condition, data = cat_data, mean)
# 
# anova(art(CAT_inhib_percent ~ Smoke * Tissue, data = cat_data))
# 
# dunn.test::dunn.test(cat_data$CAT_inhib_percent, cat_data$Condition, list = T, method = "hs")

# ADP kinetics ----

adp_mm <- pivot_longer(cat_data, 8:13, values_to = "Rate", names_to = "Concentration")

adp_mm <- na.omit(adp_mm[-c(7:25)])

adp_mm <- separate(adp_mm, Concentration, into = c("ADP", "Concentration"), sep = " ")

adp_mm$Concentration[is.na(adp_mm$Concentration)] <- 0

adp_mm$Concentration <- as.numeric( adp_mm$Concentration)

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

# CAT Kinetics -----
cat_mm <- pivot_longer(cat_data, 15:20, values_to = "Rate", names_to = "Concentration")

cat_mm <- na.omit(cat_mm[-c(7:25)])

cat_mm <- separate(cat_mm, Concentration, into = c("cat", "Concentration"), sep = " ")

cat_mm$Concentration[is.na(cat_mm$Concentration)] <- 0

cat_mm$Concentration <- as.numeric( cat_mm$Concentration)

# inhib_kinetics <- nls(Rate ~ (Vmax *
#                                 Concentration / (Km + Concentration*(1+Concentration/Ks))),
#                       data= cat_mm, start=list(Km=0.025,
#                                               Vmax=11, Ks=1))
# summary(inhib_kinetics)

cat_kinetic_values <- data.frame()

for(y in unique(cat_mm$Smoke)){
  for(x in unique(cat_mm$Tissue)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      cat_model <- drm(Rate ~ Concentration, data = subset(subset(cat_mm, Smoke %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower Limit", "Upper Limit", "IC50")))
      
      cat_model$Tissue <- paste(x)
      
      cat_model$Smoke <- paste(y)
      
      coefs <- data.frame("Imax" = c(cat_model$coefficients[1]), "Vmax" = c(cat_model$coefficients[2]), "IC50" = c(cat_model$coefficients[3]))
      
      coefs$Tissue <- paste(x)
      
      coefs$Smoke <- paste(y)
      
      cat_kinetic_values <- rbind(cat_kinetic_values, coefs)
      
      mml <- data.frame(S = seq(0, max(cat_mm$Concentration), length.out = 100))
      
      mml$v <- predict(cat_model, newdata =  mml)
      
      # modelselect <- mselect(cat_model, fctList = list(MM.2(), MM.3(), LL.3(), LL.4(), W2.3(), W2.4(), W1.3(), W1.4()))
      # fit <- modelFit(cat_model)
      # 
      # assign(paste(paste(x), "cat_summary", paste(y), sep = "_"), summary(cat_model))
      # assign(paste(paste(x), "cat_fit", paste(y), sep = "_"), fit)
      # assign(paste(paste(x), "cat_select", paste(y), sep = "_"), modelselect)
      assign(paste(paste(x), "cat_residuals", sep = "_"), residuals(cat_model))
      assign(paste(paste(x), "cat_model", paste(y), sep = "_"), cat_model)
      assign(paste(paste(x), "cat_mml", paste(y), sep = "_"), mml)
      
      
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  }
}


Gastrocnemius_cat_mml_Control$Smoke <- "Control"

Gastrocnemius_cat_mml_Smoke$Smoke <- "Smoke"

Gastrocnemius_cat_mml_Control$Tissue <- "Gastrocnemius"
Gastrocnemius_cat_mml_Smoke$Tissue <- "Gastrocnemius"

# Soleus_cat_mml_Control$Smoke <- "Control"

# Soleus_cat_mml_Smoke$Smoke <- "Smoke"

Soleus_cat_mml_Control$Tissue <- "Soleus"
# Soleus_cat_mml_Smoke$Tissue <- "Soleus"

# cat_mml <- rbind(Gastrocnemius_cat_mml_Control, Gastrocnemius_cat_mml_Smoke, Soleus_cat_mml_Control, Soleus_cat_mml_Smoke)

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

# individ_adp_kinetics$Sex <- cat_data$Sex
individ_adp_mml$Condition <- paste(individ_adp_mml$Tissue, individ_adp_mml$Smoke)


adp_km <- ggplot(data = individ_adp_kinetics, aes(x = Tissue, y = Km, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 500)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 500, 50)) +
  labs(y= expression(bold(atop('[ADP]', '(μM)')))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank())

adp_vmax <- ggplot(data = individ_adp_kinetics, aes(x = Tissue, y = Vmax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 90)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(10, 90, 10)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank()) +
  stat_compare_means(label.y = .43, label = "p.format", comparisons = smoke_comparisons, bracket.size = 1)


## Individual CAT ----
individ_cat_kinetics <- data.frame()
individ_cat_mml <- data.frame()

for(z in unique(cat_mm$Subject)){
  for(y in unique(cat_mm$Smoke)){
    for(x in unique(cat_mm$Tissue)){
      skip_to_next <- FALSE
      
      # Note that print(b) fails since b doesn't exist
      
      tryCatch({
        cat_model <- drm(Rate ~ Concentration, data = subset(subset(subset(cat_mm, Subject %in% paste(z)), Smoke %in% paste(y)), Tissue %in% paste(x)), fct = MM.3(names = c("Lower Limit", "Upper Limit", "IC50")))
        
        cat_model$Tissue <- paste(x)
        
        cat_model$Smoke <- paste(y)
        
        cat_model$Subject <- paste(z)
        
        coefs <- data.frame("C" = c(cat_model$coefficients[1]), "Imax" = c(cat_model$coefficients[2]), "IC50" = c(cat_model$coefficients[3]))
        
        coefs$Inhibition <- -(1-coefs$Imax/coefs$C) 
        
        coefs$Tissue <- paste(x)
        
        coefs$Smoke <- paste(y)
        
        coefs$Subject <- paste(z)
        
        individ_cat_kinetics <- rbind(individ_cat_kinetics, coefs)
        
        mml <- data.frame(S = seq(0, max(cat_mm$Concentration), length.out = 100))
        
        mml$v <- predict(cat_model, newdata =  mml)
        
        mml$Tissue <- paste(x) 
        
        mml$Smoke <- paste(y)
        
        mml$Subject <- paste(z)
        
        individ_cat_mml <- rbind(individ_cat_mml, mml)
        
        
      }, error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { next }
    }
  }
}

# individ_cat_kinetics$Sex <- cat_data$Sex
individ_cat_mml$Condition <- paste(individ_cat_mml$Tissue, individ_cat_mml$Smoke)


cat_imax <- ggplot(data = individ_cat_kinetics, aes(x = Tissue, y = Imax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 70)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(10, 70, 10)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(I[max]))) +
  theme(axis.title.x = element_blank())


cat_km <- ggplot(data = individ_cat_kinetics, aes(x = Tissue, y = IC50, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 3)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 10, 0.5)) +
  labs(y= expression(bold(atop('[CAT]', '(μM)')))) +
  ggtitle(expression(bold(IC[50]))) +
  theme(axis.title.x = element_blank())

cat_vmax <- ggplot(data = individ_cat_kinetics, aes(x = Tissue, y = Vmax, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 150)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(10, 150, 20)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank())

cat_percent <- ggplot(data = individ_cat_kinetics, aes(x = Tissue, y = Inhibition, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(-1, 0)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1, 0, 0.1), labels = scales::percent_format()) +
  labs(y=expression(bold(Inhibition~('%'~V[max])))) +
  ggtitle(expression(bold(CAT~Inhibition))) +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(position = "top")



## Main Plot ----
mainplot_gast <- ggplot(data = subset(subset(cat_long, Tissue %in% "Gastrocnemius"), State %in% c("Basal", 'GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
                                                                                                  , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 120)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(CAT~Inhibition))) +
  theme(axis.title.x = element_blank())

mainplot_sol <- ggplot(data = subset(subset(cat_long, Tissue %in% "Soleus"), State %in% c("Basal", 'GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
                                                                                          , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, fill = Smoke)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width=0.9), width = 0.3, size = 1, na.rm = TRUE) +
  stat_summary(geom = "bar", position = position_dodge(0.9), color = "black", size = 1) +
  theme_prism() +
  scale_fill_prism() +
  coord_cartesian(ylim = c(0, 120)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 20)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(CAT~Inhibition))) +
  theme(axis.title.x = element_blank())


#stats ----

individ_cat_kinetics$Tissue <- factor(individ_cat_kinetics$Tissue)
individ_cat_kinetics$Smoke <- factor(individ_cat_kinetics$Smoke)

individ_adp_kinetics$Tissue <- factor(individ_adp_kinetics$Tissue)
individ_adp_kinetics$Smoke <- factor(individ_adp_kinetics$Smoke)

individ_adp_kinetics$Km[individ_adp_kinetics$Km == individ_adp_kinetics[19,2]] <- NA
individ_adp_kinetics$Km[individ_adp_kinetics$Km == individ_adp_kinetics[15,2]] <- NA
individ_adp_kinetics$Km[individ_adp_kinetics$Km == individ_adp_kinetics[4,2]] <- NA

# for(x in unique(individ_adp_kinetics$Smoke)){
#   for(y in unique(individ_adp_kinetics$Tissue)){
#       print(
#         identify_outliers(subset(subset(individ_adp_kinetics, Smoke %in% paste(x)), Tissue %in% paste(y)), Vmax)
# 
#       )
# 
# 
# 
#   }
# }
# 
# for(x in unique(individ_cat_kinetics$Smoke)){
#   for(y in unique(individ_cat_kinetics$Tissue)){
#     print(
#       identify_outliers(subset(subset(individ_cat_kinetics, Smoke %in% paste(x)), Tissue %in% paste(y)), Inhibition)
# 
#     )
# 
# 
# 
#   }
# }
resp_anova <- anova(art(Rate ~ State * Smoke * Tissue, data = na.omit(cat_long)))

# summary(aov(Inhibition ~ Smoke * Tissue, data = na.omit(individ_cat_kinetics)))
# 
# TukeyHSD(aov(Inhibition ~ Smoke * Tissue, data = na.omit(individ_cat_kinetics)))
# 
# ## Normality & homogeneity ----
# individ_adp_kinetics |>
#   group_by(Tissue) |>
#   levene_test(Vmax ~ Smoke)
# 
# as.data.table(individ_adp_kinetics)[,.(Statistic = shapiro.test(Vmax)$statistic,
#                                         P.value = shapiro.test(Vmax)$p.value),
#                                      by = .(Tissue, Smoke)]
# 
# individ_adp_kinetics |>
#   group_by(Tissue) |>
#   levene_test(Km ~ Smoke)
# 
# as.data.table(individ_adp_kinetics)[,.(Statistic = shapiro.test(Km)$statistic,
#                                        P.value = shapiro.test(Km)$p.value),
#                                     by = .(Tissue, Smoke)]
# 
# individ_cat_kinetics |>
#   group_by(Tissue) |>
#   levene_test(Imax ~ Smoke)
# 
# as.data.table(individ_cat_kinetics)[,.(Statistic = shapiro.test(Imax)$statistic,
#                                        P.value = shapiro.test(Imax)$p.value),
#                                     by = .(Tissue, Smoke)]
# 
# individ_cat_kinetics |>
#   group_by(Tissue) |>
#   levene_test(IC50 ~ Smoke)
# 
# as.data.table(individ_cat_kinetics)[,.(Statistic = shapiro.test(IC50)$statistic,
#                                        P.value = shapiro.test(IC50)$p.value),
#                                     by = .(Tissue, Smoke)]
# 
# individ_cat_kinetics |>
#   group_by(Tissue) |>
#   levene_test(Inhibition ~ Smoke)
# 
# as.data.table(individ_cat_kinetics)[,.(Statistic = shapiro.test(Inhibition)$statistic,
#                                        P.value = shapiro.test(Inhibition)$p.value),
#                                     by = .(Tissue, Smoke)]


# main_results <- cat_long |> 
#   subset(State %in% c("Basal", 'GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
#                            , 'Rot', 'AmA & Omy')) |> 
#   group_by(State, Tissue) |> 
#   t_test(Rate ~ Smoke)# |> 
# # adjust_pvalue(method = "BH") |> 
# # add_significance("p.adj")
# 
# Gastroc_main <- subset(main_results, Tissue %in% "Gastrocnemius") |> 
#   add_xy_position(x = "State", dodge = 0.9)
# 
# Gastroc_main$xmin <- seq(0.8, 7.8, 1)
# Gastroc_main$xmax <- seq(1.2, 8.2, 1)
# 
# Soleus_main <- subset(main_results, Tissue %in% "Soleus")|> 
#   add_xy_position(x = "State", dodge = 0.9)
# 
# Soleus_main$xmin <- seq(0.8, 7.8, 1)
# Soleus_main$xmax <- seq(1.2, 8.2, 1)
# 
# 
# adp_km_results <- individ_adp_kinetics |> 
#   group_by(Tissue) |> 
#   t_test(Km ~ Smoke) |> 
#   add_xy_position(x = "Tissue", dodge = 0.9)
# 
# adp_vmax_results <- individ_adp_kinetics |> 
#   group_by(Tissue) |> 
#   t_test(Vmax ~ Smoke) |> 
#   add_xy_position(x = "Tissue", dodge = 0.9)
# 
# cat_imax_results <- individ_cat_kinetics |> 
#   group_by(Tissue) |> 
#   t_test(Imax ~ Smoke) |> 
#   add_xy_position(x = "Tissue", dodge = 0.9)
# 
# cat_IC50_results <- individ_cat_kinetics[-c(25),] |> 
#   group_by(Tissue) |> 
#   t_test(IC50 ~ Smoke) |> 
#   add_xy_position(x = "Tissue", dodge = 0.9)
# 
# cat_inhib_results <- individ_cat_kinetics |> 
#   group_by(Tissue) |> 
#   t_test(Inhibition ~ Smoke) |> 
#   add_xy_position(x = "Tissue", dodge = 0.9)

# pairwise.t.test(individ_adp_kinetics$Vmax, c(individ_adp_kinetics$Tissue, individ_adp_kinetics$Smoke))
# 
# individ_adp_kinetics$Condition <- paste(individ_adp_kinetics$Tissue, individ_adp_kinetics$Smoke)
# individ_cat_kinetics$Condition <- paste(individ_cat_kinetics$Tissue, individ_cat_kinetics$Smoke)

# Graphs ----
## ADP ----

individ_adp_kinetics$Condition <- paste(individ_adp_kinetics$Tissue, individ_adp_kinetics$Smoke)

Km_anova <- anova(art(Km ~ Smoke * Tissue, data = na.omit(individ_adp_kinetics)))
effectsize::eta_squared(Km_anova)
aggregate(x = Km ~ Smoke + Tissue, data = individ_adp_kinetics, FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

adp_km_dunn <- dunn.test(individ_adp_kinetics$Km, individ_adp_kinetics$Condition, list = TRUE, method = "hs", kw = F)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Tissue %in% "Gastrocnemius"), formula = Km ~ Smoke)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Tissue %in% "Soleus"), formula = Km ~ Smoke)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Smoke %in% "Control"), formula = Km ~ Tissue)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Smoke %in% "Smoke"), formula = Km ~ Tissue)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Condition %in% c("Gastrocnemius Smoke", "Soleus Control")), formula = Km ~ Condition)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Condition %in% c("Gastrocnemius Control", "Soleus Smoke")), formula = Km ~ Condition)

adp_km <- ggplot(data = individ_adp_kinetics, aes(x = Condition, y = Km, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(values = c("blue", "red", "blue4", "red4")) +
  coord_cartesian(ylim = c(0, 900)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 900, 100)) +
  labs(y= expression(bold(atop('[ADP]', '(μM)')))) +
  ggtitle(expression(bold(K[m]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 880, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 840, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 800, label = 'Interaction: p < 0.001', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 790, inherit.aes = F, label = "p = 0.042", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Soleus Control", y.position = 660, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Soleus Control", y.position = 730, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.01)

Vmax_anova <- anova(art(Vmax ~ Smoke * Tissue, data = na.omit(individ_adp_kinetics)))
effectsize::eta_squared(Vmax_anova)
aggregate(x = Vmax ~ Smoke + Tissue, data = individ_adp_kinetics, FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

adp_vmax_dunn <- dunn.test(individ_adp_kinetics$Vmax, individ_adp_kinetics$Condition, list = TRUE, method = "none", kw = F)
Sidak(c(adp_vmax_dunn$P[1], adp_vmax_dunn$P[6]))
Sidak(c(adp_vmax_dunn$P[2], adp_vmax_dunn$P[5]))
rstatix::cohens_d(data = subset(individ_adp_kinetics, Tissue %in% "Gastrocnemius"), formula = Vmax ~ Smoke)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Tissue %in% "Soleus"), formula = Vmax ~ Smoke)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Smoke %in% "Control"), formula = Vmax ~ Tissue)
rstatix::cohens_d(data = subset(individ_adp_kinetics, Smoke %in% "Smoke"), formula = Vmax ~ Tissue)

adp_vmax <- ggplot(data = individ_adp_kinetics, aes(x = Condition, y = Vmax, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(values = c("blue", "red", "blue4", "red4")) +
  coord_cartesian(ylim = c(0, 160)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 160, 20)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 157, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 150, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 143, label = 'Interaction: p = 0.813', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 145, inherit.aes = F, label = "p = 0.013", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Soleus Control", y.position = 120, inherit.aes = F, label = "p = 0.028", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Soleus Smoke", y.position = 133, inherit.aes = F, label = "p = 0.032", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Gastrocnemius Control", y.position = 66, inherit.aes = F, label = "p = 0.014", size = 1, fontface = 2, tip.length = 0.02)

aggregate(x = Vmax ~ Smoke + Tissue, data = individ_adp_kinetics, FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
# identify_outliers(data = subset(individ_adp_kinetics, Condition %in% "Soleus Control"), Vmax)

## CAT ----
individ_cat_kinetics$Condition <- paste(individ_cat_kinetics$Tissue, individ_cat_kinetics$Smoke)

Imax_anova <- anova(art(Imax ~ Smoke * Tissue, data = na.omit(individ_cat_kinetics[-c(28,36),])))
effectsize::eta_squared(Imax_anova)

aggregate(x = Imax ~ Smoke + Tissue, data = individ_cat_kinetics[-c(28,36),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

cat_imax <- individ_cat_kinetics[-c(28,36),]
Imax_dunn <- dunn.test(cat_imax$Imax, cat_imax$Condition, list = T, method = "none", kw = F)
Sidak(c(Imax_dunn$P[2], Imax_dunn$P[5]))
Sidak(c(Imax_dunn$P[1], Imax_dunn$P[6]))

rstatix::cohens_d(data = subset(cat_imax, Smoke %in% "Control"), formula = Imax ~ Tissue)
rstatix::cohens_d(data = subset(cat_imax, Smoke %in% "Smoke"), formula = Imax ~ Tissue)

cat_imax <- ggplot(data = individ_cat_kinetics, aes(x = Condition, y = Imax, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(values = c("blue", "red", "blue4", "red4")) +
  coord_cartesian(ylim = c(0, 180)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 180, 20)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(I[max]))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 176, label = 'Smoke: p = 0.059', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 168, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 160, label = 'Interaction: p = 0.506', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Soleus Control", y.position = 132, inherit.aes = F, label = "p = 0.002", size = 1, fontface = 2, tip.length = 0.02) +
  geom_bracket(xmin = "Gastrocnemius Smoke", xmax = "Soleus Smoke", y.position = 147, inherit.aes = F, label = "p = 0.007", size = 1, fontface = 2, tip.length = 0.02) 


IC50_anova <- anova(art(IC50 ~ Smoke * Tissue, data = na.omit(individ_cat_kinetics[-c(36, 7),])))
effectsize::eta_squared(IC50_anova)

aggregate(x = IC50 ~ Smoke + Tissue, data = individ_cat_kinetics[-c(36, 7),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

cat_km <- ggplot(data = individ_cat_kinetics[-c(36, 7),], aes(x = Condition, y = IC50, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(values = c("blue", "red", "blue4", "red4")) +
  coord_cartesian(ylim = c(0, 10)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 10, 1)) +
  labs(y= expression(bold('[CAT] (μM)'))) +
  ggtitle(expression(bold(IC['50']))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  annotate('text', x = 0.55, y = 9.8, label = 'Smoke: p = 0.172', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 9.3, label = 'Tissue: p = 0.715', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = 8.8, label = 'Interaction: p = 0.290', hjust = 0, fontface = 2)


cat_vmax_anova <- anova(art(C ~ Smoke * Tissue, data = na.omit(individ_cat_kinetics)))
effectsize::eta_squared(cat_vmax_anova)
aggregate(x = C ~ Smoke + Tissue, data = individ_cat_kinetics[-c(36, 7),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

cat_Vmax_dunn <- dunn.test(individ_cat_kinetics$C, individ_cat_kinetics$Condition, list = T, method = "none", kw = F)
Sidak(c(cat_Vmax_dunn$P[2], cat_Vmax_dunn$P[5]))
Sidak(c(cat_Vmax_dunn$P[1], cat_Vmax_dunn$P[6]))
rstatix::cohens_d(data = subset(individ_cat_kinetics, Tissue %in% "Gastrocnemius"), formula = C ~ Smoke)
rstatix::cohens_d(data = subset(individ_cat_kinetics, Tissue %in% "Soleus"), formula = C ~ Smoke)
rstatix::cohens_d(data = subset(individ_cat_kinetics, Smoke %in% "Control"), formula = C ~ Tissue)
rstatix::cohens_d(data = subset(individ_cat_kinetics, Smoke %in% "Smoke"), formula = C ~ Tissue)

cat_vmax <- ggplot(data = individ_cat_kinetics, aes(x = Condition, y = C, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(values = c("blue", "red", "blue4", "red4")) +
  coord_cartesian(ylim = c(0, 150)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 150, 25)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(V[max]))) +
  theme(axis.title.x = element_blank())


Inhibition_anova <- anova(art(Inhibition ~ Smoke * Tissue, data = na.omit(individ_cat_kinetics[-c(7, 24),])))
effectsize::eta_squared(Inhibition_anova)

aggregate(x = Inhibition ~ Smoke + Tissue, data = individ_cat_kinetics[-c(7, 24),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

cat_inhib <- individ_cat_kinetics[-c(7, 24),]

cat_Inhibition_dunn <- dunn.test(cat_inhib$Inhibition, cat_inhib$Condition, list = T, method = "hs", kw = F)
rstatix::cohens_d(data = subset(individ_cat_kinetics[-c(7, 24),], Tissue %in% "Gastrocnemius"), formula = Inhibition ~ Smoke)
rstatix::cohens_d(data = subset(individ_cat_kinetics[-c(7, 24),], Tissue %in% "Soleus"), formula = Inhibition ~ Smoke)
rstatix::cohens_d(data = subset(individ_cat_kinetics[-c(7, 24),], Smoke %in% "Control"), formula = Inhibition ~ Tissue)
rstatix::cohens_d(data = subset(individ_cat_kinetics[-c(7, 24),], Smoke %in% "Smoke"), formula = Inhibition ~ Tissue)
rstatix::cohens_d(data = subset(individ_cat_kinetics[-c(7, 24),], Condition %in% c("Gastrocnemius Smoke", "Soleus Control")), formula = Inhibition ~ Condition)
rstatix::cohens_d(data = subset(individ_cat_kinetics[-c(7, 24),], Condition %in% c("Gastrocnemius Control", "Soleus Smoke")), formula = Inhibition ~ Condition)


cat_percent <- ggplot(data = individ_cat_kinetics[-c(7, 24),], aes(x = Condition, y = Inhibition, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(values = c("blue", "red", "blue4", "red4")) +
  coord_cartesian(ylim = c(-1.3, 0)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-1.3, 0, 0.1), labels = scales::percent_format()) +
  labs(y=expression(bold(Inhibition~('%'~V[max])))) +
  ggtitle(expression(bold(CAT~Inhibition))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_x_discrete(position = "top") +
  annotate('text', x = 0.55, y = -1.14, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = -1.2, label = 'Tissue: p = 0.067', hjust = 0, fontface = 2) +
  annotate('text', x = 0.55, y = -1.26, label = 'Interaction: p = 0.007', hjust = 0, fontface = 2) +
  geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = -1, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = -0.02, vjust = 1.8) +
  geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = -.73, inherit.aes = F, label = "p = 0.003", size = 1, fontface = 2, tip.length = -0.02, vjust = 1.8) 


## Main Plot ----
mainplot_gast <- ggplot(data = subset(subset(cat_long, Tissue %in% "Gastrocnemius"), State %in% c("Basal", 'GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
                                                                                                  , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, 10)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(Gastrocnemius))) +
  theme(axis.title.x = element_blank())

mainplot_sol <- ggplot(data = subset(subset(cat_long, Tissue %in% "Soleus"), State %in% c("Basal", 'GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
                                                                                          , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, color = Smoke)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_prism() +
  coord_cartesian(ylim = c(0, 175)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 175, 25)) +
  labs(y=expression(bold(atop(Respiration~Rate,(pmol[O[2]]/sec/mg[wt]))))) +
  ggtitle(expression(bold(Soleus))) +
  theme(axis.title.x = element_blank())

cat_long$Condition <- paste(cat_long$Tissue, cat_long$Smoke)

# main_results_all <- cat_long |> 
#   subset(State %in% c("Basal", 'GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
#                       , 'Rot', 'AmA & Omy')) |> 
#   group_by(State) |> 
#   t_test(Rate ~ Condition)


mainplot_all <- ggplot(data = subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], State %in% c("Basal", 'GM', 'GMD 5000', 'GMDS', "CAT 5.0", 'FCCP Peak'
                                                                                                       , 'Rot', 'AmA & Omy')), aes(x = State, y = Rate, color = Condition)) +
  stat_summary(geom = "bar", position = position_dodge(0.9), fill = NA, size = 1) +
  geom_point(position = position_dodge(0.9), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  theme_prism() +
  scale_color_manual(values = c("blue", "red", "blue4", "red4")) +
  coord_cartesian(ylim = c(0, 200)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 200, 25)) +
  labs(y=expression(bold(Respiration~Rate~(pmol[O[2]]/sec/mg[wt])))) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5), linetype = "longdash", color = "grey50", size = .6) +
  #Basal
  annotate('text', x = 0.45, y = 198, label = 'Smoke: p = 0.408', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 191, label = 'Tissue: p = 0.283', hjust = 0, fontface = 2) +
  annotate('text', x = 0.45, y = 184, label = 'Interaction: p = 0.829', hjust = 0, fontface = 2) +
  #GM
  annotate('text', x = 1.51, y = 198, label = 'Smoke: p = 0.068', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 1.51, y = 184, label = 'Interaction: p = 0.406', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(1.665), xmax = c(2.115), y.position = 55, inherit.aes = F, label = "p = 0.048", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(1.885), xmax = c(2.345), y.position = 70, inherit.aes = F, label = "p = 0.019", size = 1, fontface = 2, tip.length = 0.01) +
  #GMD
  annotate('text', x = 2.51, y = 198, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 191, label = 'Tissue: p = 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 2.51, y = 184, label = 'Interaction: p = 0.604', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(2.665), xmax = c(2.885), y.position = 70, inherit.aes = F, label = "p = 0.014", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(3.115), xmax = c(3.345), y.position = 135, inherit.aes = F, label = "p = 0.038", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(2.665), xmax = c(3.115), y.position = 125, inherit.aes = F, label = "p = 0.048", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(2.885), xmax = c(3.345), y.position = 85, inherit.aes = F, label = "p = 0.019", size = 1, fontface = 2, tip.length = 0.01) +
  #GMDS
  annotate('text', x = 3.51, y = 198, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 3.51, y = 184, label = 'Interaction: p = 0.604', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(3.665), xmax = c(3.885), y.position = 75, inherit.aes = F, label = "p = 0.052", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(4.115), xmax = c(4.345), y.position = 160, inherit.aes = F, label = "p = 0.080", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(3.665), xmax = c(4.115), y.position = 150, inherit.aes = F, label = "p = 0.012", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(3.885), xmax = c(4.345), y.position = 125, inherit.aes = F, label = "p = 0.008", size = 1, fontface = 2, tip.length = 0.01) +
  #CAT
  annotate('text', x = 4.51, y = 198, label = 'Smoke: p = 1.000', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 4.51, y = 184, label = 'Interaction: p = 0.835', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(4.665), xmax = c(5.115), y.position = 130, inherit.aes = F, label = "p < 0.001", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(4.885), xmax = c(5.345), y.position = 140, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.01) +
  #FCCP
  annotate('text', x = 5.51, y = 198, label = 'Smoke: p = 0.004', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 5.51, y = 184, label = 'Interaction: p = 0.857', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(5.665), xmax = c(6.115), y.position = 162, inherit.aes = F, label = "p = 0.008", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(5.885), xmax = c(6.345), y.position = 172, inherit.aes = F, label = "p = 0.015", size = 1, fontface = 2, tip.length = 0.01) +
  #Rot
  annotate('text', x = 6.51, y = 198, label = 'Smoke: p = 0.022', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 6.51, y = 184, label = 'Interaction: p = 0.468', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(6.665), xmax = c(7.115), y.position = 100, inherit.aes = F, label = "p = 0.002", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(6.885), xmax = c(7.345), y.position = 110, inherit.aes = F, label = "p = 0.015", size = 1, fontface = 2, tip.length = 0.01) +
  #Ama & Omy
  annotate('text', x = 7.51, y = 198, label = 'Smoke: p = 0.086', hjust = 0, fontface = 2) +
  annotate('text', x = 7.51, y = 191, label = 'Tissue: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 7.51, y = 184, label = 'Interaction: p = 0.750', hjust = 0, fontface = 2) +
  geom_bracket(xmin = c(7.665), xmax = c(8.115), y.position = 60, inherit.aes = F, label = "p = 0.008", size = 1, fontface = 2, tip.length = 0.01) +
  geom_bracket(xmin = c(7.885), xmax = c(8.345), y.position = 70, inherit.aes = F, label = "p = 0.010", size = 1, fontface = 2, tip.length = 0.01)

# for(y in unique(cat_long$State)){
#   print(subset(cat_long, State %in% paste(y)) |>
#     group_by(Tissue) |>
#     levene_test(Rate ~ Smoke))
# 
#   print(as.data.table(subset(cat_long, State %in% paste(y)))[,.(Statistic = shapiro.test(Rate)$statistic,
#                                          P.value = shapiro.test(Rate)$p.value),
#                                       by = .(Tissue, Smoke)])
# 
# }
# 
# anova_states_results <- data.frame()
# anova_eff_results <- data.frame()
# 
# for(y in unique(cat_long$State)){
#   anova_states <- anova(art(Rate ~ Smoke * Tissue, data = subset(na.omit(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),]), State %in% paste(y))))
#   anova_states$State <- paste(y)
# 
#   anova_eff <- effectsize::eta_squared(anova_states)
# 
#   print(anova_states)
#   print(effectsize::eta_squared(anova_states))
# 
#   anova_states_results <- rbind(anova_states_results, anova_states)
#   anova_eff_results <- rbind(anova_eff_results, anova_eff)
# 
# }
# 
# cbind(anova_states_results, anova_eff_results)
# 
# for(y in unique(cat_long$State)){
#   dunn <- subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], State %in% paste(y))
# 
#   assign(paste("dunn", paste(y), sep = "_"), dunn)
# 
# }
# 
# 
# GM_dunn <- dunn.test(dunn_GM$Rate, dunn_GM$Condition, method = "none", list = T)
# Sidak(c(GM_dunn$P[2], GM_dunn$P[5]))
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Control"), State %in% "GM"), formula = Rate ~ Tissue)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Smoke"), State %in% "GM"), formula = Rate ~ Tissue)
# 
# GMD_dunn <- dunn.test(`dunn_GMD 5000`$Rate, `dunn_GMD 5000`$Condition, method = "none", list = T)
# Sidak(c(GMD_dunn$P[2], GMD_dunn$P[5]))
# Sidak(c(GMD_dunn$P[1], GMD_dunn$P[6]))
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Tissue %in% "Gastrocnemius"), State %in% "GMD 5000"), formula = Rate ~ Smoke)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Tissue %in% "Soleus"), State %in% "GMD 5000"), formula = Rate ~ Smoke)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Control"), State %in% "GMD 5000"), formula = Rate ~ Tissue)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Smoke"), State %in% "GMD 5000"), formula = Rate ~ Tissue)
# 
# GMDS_dunn <- dunn.test(dunn_GMDS$Rate, dunn_GMDS$Condition, method = "none", list = T)
# Sidak(c(GMDS_dunn$P[2], GMDS_dunn$P[5]))
# Sidak(c(GMDS_dunn$P[1], GMDS_dunn$P[6]))
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Tissue %in% "Gastrocnemius"), State %in% "GMDS"), formula = Rate ~ Smoke)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Tissue %in% "Soleus"), State %in% "GMDS"), formula = Rate ~ Smoke)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Control"), State %in% "GMDS"), formula = Rate ~ Tissue)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Smoke"), State %in% "GMDS"), formula = Rate ~ Tissue)
# 
# CAT_dunn <- dunn.test(`dunn_CAT 5.0`$Rate, `dunn_CAT 5.0`$Condition, method = "none", list = T)
# Sidak(c(CAT_dunn$P[2], CAT_dunn$P[5]))
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Control"), State %in% "CAT 5.0"), formula = Rate ~ Tissue)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Smoke"), State %in% "CAT 5.0"), formula = Rate ~ Tissue)
# 
# FCCP_dunn <- dunn.test(`dunn_FCCP Peak`$Rate, `dunn_FCCP Peak`$Condition, method = "none", list = T)
# Sidak(c(FCCP_dunn$P[2], FCCP_dunn$P[5]))
# Sidak(c(FCCP_dunn$P[1], FCCP_dunn$P[6]))
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Tissue %in% "Gastrocnemius"), State %in% "FCCP Peak"), formula = Rate ~ Smoke)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Tissue %in% "Soleus"), State %in% "FCCP Peak"), formula = Rate ~ Smoke)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Control"), State %in% "FCCP Peak"), formula = Rate ~ Tissue)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Smoke"), State %in% "FCCP Peak"), formula = Rate ~ Tissue)
# 
# Rot_dunn <- dunn.test(dunn_Rot$Rate, dunn_Rot$Condition, method = "none", list = T)
# Sidak(c(Rot_dunn$P[2], Rot_dunn$P[5]))
# Sidak(c(Rot_dunn$P[1], Rot_dunn$P[6]))
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Tissue %in% "Gastrocnemius"), State %in% "Rot"), formula = Rate ~ Smoke)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Tissue %in% "Soleus"), State %in% "Rot"), formula = Rate ~ Smoke)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Control"), State %in% "Rot"), formula = Rate ~ Tissue)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Smoke"), State %in% "Rot"), formula = Rate ~ Tissue)
# 
# 
# AA_dunn <- dunn.test(`dunn_AmA & Omy`$Rate, `dunn_AmA & Omy`$Condition, method = "none", list = T)
# Sidak(c(AA_dunn$P[2], AA_dunn$P[5]))
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Control"), State %in% "AmA & Omy"), formula = Rate ~ Tissue)
# rstatix::cohens_d(data = subset(subset(cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], Smoke %in% "Smoke"), State %in% "AmA & Omy"), formula = Rate ~ Tissue)
# 
# 
# resp_descriptives <- aggregate(x = Rate ~ Smoke + Tissue + State, data = cat_long[-c(69, 85, 578, 595, 597, 612, 614, 629),], FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
# subset(resp_descriptives, State %in% "GMD 5000")
# 
# outliers <- data.frame()
# 
# for(y in unique(cat_long$Condition)){
#   for(x in names(cat_data[c(7,8,13,14,20:23)])){
#     skip_to_next <- FALSE
#     
#     # Note that print(b) fails since b doesn't exist
#     
#     tryCatch({
#     out <- identify_outliers(data = subset(subset(cat_long, Condition %in% paste(y)), State %in% paste(x)), Rate)
#     
#     # resp_outliers <- outliers[which(isTRUE(outliers$is.extreme))]
#     # 
#     # assign(paste("outliers", x, y, sep = "_"), resp_outliers)
#     
#     if(isTRUE(outliers$is.extreme)){ 
#        print(paste(outliers[which(outliers$is.extreme == TRUE)], y, x, sep = "_"))
#     }
#     
#     b <- data.frame(c(out[which(out$is.extreme == TRUE)]))
#     
#     a <- data.frame(b[1])
#     
#     a$State <- paste(x)
#     
#     a$Condition <- paste(y)
#     
#     outliers <- rbind(outliers, a)
#     }, error = function(e) { skip_to_next <<- TRUE})
#     
#     if(skip_to_next) { next }
#   }
# }


## New ETC analysis

etc_analysis <- pivot_wider(na.omit(cat_long), id_cols = Subject, values_from = Rate, names_from = c(Condition, State))

etc_analysis$Gastrocnemius_difference <- (etc_analysis$`Gastrocnemius Smoke_FCCP Peak` - etc_analysis$`Gastrocnemius Control_FCCP Peak`)/etc_analysis$`Gastrocnemius Control_FCCP Peak`

etc_analysis$Soleus_difference <- (etc_analysis$`Soleus Smoke_FCCP Peak` - etc_analysis$`Soleus Control_FCCP Peak`)/etc_analysis$`Soleus Control_FCCP Peak`

etc_analysis$Gastrocnemius_percentage <- (((etc_analysis$`Gastrocnemius Smoke_GMDS` - etc_analysis$`Gastrocnemius Smoke_CAT 5.0`) - (etc_analysis$`Gastrocnemius Control_GMDS` - etc_analysis$`Gastrocnemius Control_CAT 5.0`))/etc_analysis$`Gastrocnemius Control_GMDS`)

etc_analysis$Soleus_percentage <- (((etc_analysis$`Soleus Smoke_GMDS` - etc_analysis$`Soleus Smoke_CAT 5.0`) - (etc_analysis$`Soleus Control_GMDS` - etc_analysis$`Soleus Control_CAT 5.0`))/etc_analysis$`Soleus Control_GMDS`)

etc_analysis_long <- pivot_longer(etc_analysis, cols = 2:73, names_to = c("Condition", "State"), names_sep = "_")

etc_analysis_long <- na.omit(etc_analysis_long)

FCCP_diff <- ggplot(data = subset(etc_analysis_long, State %in% "difference"), aes(x = Condition, y = value, color = Condition)) +
  # stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, color = "black") +
  geom_bar(aes(x = Condition, y = value, fill = Condition), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(values = c("grey50", "black")) +
  coord_cartesian(ylim = c(-.75, .5), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-.75, .5, .25), labels = label_percent()
  ) +
  labs(y=expression(bold(Delta~FCCP[Peak]))) +
  #ggtitle(expression(bold(Delta~FCCP[Peak]~"(Smoke"-"Control)"))) +
  # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = CIIPercent, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
  geom_beeswarm(dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3, fill = "white") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  geom_hline(yintercept = 0, size = 1) # +
# annotate('text', x = 3.55, y = 3.45, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
# annotate('text', x = 3.55, y = 3.32, label = 'Tissue: p = 0.006', hjust = 0, fontface = 2) +
# annotate('text', x = 3.55, y = 3.19, label = 'Interaction: p = 0.006', hjust = 0, fontface = 2) +
# geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = 3.37, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.02) +
# geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 2.5, inherit.aes = F, label = "p = 0.075", size = 1, fontface = 2, tip.length = 0.02)

identify_outliers(subset(subset(etc_analysis_long, Condition %in% "Soleus"), State %in% "difference"), value)

wilcox_test(subset(etc_analysis_long, State %in% "difference"), value ~ Condition)

cohens_d(subset(etc_analysis_long, State %in% "difference"), value ~ Condition)

aggregate(value ~ Condition, data = subset(etc_analysis_long, State %in% "difference"), FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))

CAT_diff <- ggplot(data = subset(etc_analysis_long[-c(466),], State %in% "percentage"), aes(x = Condition, y = value, color = Condition)) +
  geom_bar(position = position_dodge(0.99), stat = "summary", fill = "white", size = 1) +
  theme_prism() +
  scale_color_manual(values = c("grey50", "black")) +
  coord_cartesian(ylim = c(-.7, .1), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-.7, .1, .1), labels = label_percent()
  ) +
  labs(y=expression(bold('CSC-Induced'~ANT~Inhibition))) +
  geom_beeswarm(dodge.width = 0.99, size = 2, pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3, fill = "white") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  geom_hline(yintercept = 0, size = 1) +
  geom_bracket(xmin = "Gastrocnemius", xmax = "Soleus", y.position = -.66, inherit.aes = F, label = "p = 0.009", size = 1, fontface = 2, tip.length = -0.02, vjust = 2.5)

identify_outliers(subset(subset(etc_analysis_long, Condition %in% "Soleus"), State %in% "percentage"), value)

wilcox_test(subset(etc_analysis_long[-c(466),], State %in% "percentage"), value ~ Condition)

cohens_d(subset(etc_analysis_long, State %in% "percentage"), value ~ Condition)

aggregate(value ~ Condition, data = subset(etc_analysis_long[-c(466),], State %in% "percentage"), FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))



###### GGIRAPH ------
## New ETC analysis

# etc_analysis <- pivot_wider(na.omit(cat_long), id_cols = Subject, values_from = Rate, names_from = c(Condition, State))
# 
# etc_analysis$Gastrocnemius_difference <- (etc_analysis$`Gastrocnemius Smoke_FCCP Peak` - etc_analysis$`Gastrocnemius Control_FCCP Peak`)/etc_analysis$`Gastrocnemius Control_FCCP Peak`
# 
# etc_analysis$Soleus_difference <- (etc_analysis$`Soleus Smoke_FCCP Peak` - etc_analysis$`Soleus Control_FCCP Peak`)/etc_analysis$`Soleus Control_FCCP Peak`
# 
# etc_analysis$Gastrocnemius_percentage <- (((etc_analysis$`Gastrocnemius Smoke_GMDS` - etc_analysis$`Gastrocnemius Smoke_CAT 5.0`) - (etc_analysis$`Gastrocnemius Control_GMDS` - etc_analysis$`Gastrocnemius Control_CAT 5.0`))/etc_analysis$`Gastrocnemius Control_GMDS`)
# 
# etc_analysis$Soleus_percentage <- (((etc_analysis$`Soleus Smoke_GMDS` - etc_analysis$`Soleus Smoke_CAT 5.0`) - (etc_analysis$`Soleus Control_GMDS` - etc_analysis$`Soleus Control_CAT 5.0`))/etc_analysis$`Soleus Control_GMDS`)
# 
# etc_analysis_long <- pivot_longer(etc_analysis, cols = 2:73, names_to = c("Condition", "State"), names_sep = "_")
# 
# etc_analysis_long <- na.omit(etc_analysis_long)
# 
# ggiraph::set_girafe_defaults(opts_hover = ggiraph::opts_hover(css = ggiraph::girafe_css_bicolor(primary = "darkred", secondary = "black")))
# 
# FCCP_diff <- ggplot(data = subset(etc_analysis_long, State %in% "difference"), aes(x = Condition, y = value, color = Condition)) +
#   # stat_summary(geom = "bar", position = position_dodge(0.99), size = 1, color = "black") +
#   ggiraph::geom_bar_interactive(aes(x = Condition, y = value, fill = Condition, tooltip = value), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1) +
#   theme_prism() +
#   scale_color_manual(values = c("grey50", "black")) +
#   coord_cartesian(ylim = c(-.75, .5), clip = "off") +
#   scale_y_continuous(expand = c(0,0), breaks = seq(-.75, .5, .25), labels = label_percent()
#   ) +
#   labs(y=expression(bold(Delta~FCCP[Peak]))) +
#   #ggtitle(expression(bold(Delta~FCCP[Peak]~"(Smoke"-"Control)"))) +
#   # geom_point(data = individ_cat_kinetics, aes(x = Condition, y = CIIPercent, fill =  Condition), position = position_dodge(0.99), size = 2, pch = 21, stroke = 1.5, show.legend = FALSE) +
#   ggiraph::geom_point_interactive(pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3, fill = "white", position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0),
#                                   aes(tooltip = paste0("\U0394 FCCP<sub>Peak</sub> = ", scales::percent(round(value,3)), "%; Subject = ", Subject), data_id = Subject)) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank(),
#         axis.line.x = element_blank()) +
#   geom_hline(yintercept = 0, size = 1) # +
# # annotate('text', x = 3.55, y = 3.45, label = 'Smoke: p < 0.001', hjust = 0, fontface = 2) +
# # annotate('text', x = 3.55, y = 3.32, label = 'Tissue: p = 0.006', hjust = 0, fontface = 2) +
# # annotate('text', x = 3.55, y = 3.19, label = 'Interaction: p = 0.006', hjust = 0, fontface = 2) +
# # geom_bracket(xmin = "Gastrocnemius Control", xmax = "Gastrocnemius Smoke", y.position = 3.37, inherit.aes = F, label = "p = 0.001", size = 1, fontface = 2, tip.length = 0.02) +
# # geom_bracket(xmin = "Soleus Control", xmax = "Soleus Smoke", y.position = 2.5, inherit.aes = F, label = "p = 0.075", size = 1, fontface = 2, tip.length = 0.02)
# 
# ggiraph::girafe(ggobj = FCCP_diff)
# 
# identify_outliers(subset(subset(etc_analysis_long, Condition %in% "Soleus"), State %in% "difference"), value)
# 
# wilcox_test(subset(etc_analysis_long, State %in% "difference"), value ~ Condition)
# 
# cohens_d(subset(etc_analysis_long, State %in% "difference"), value ~ Condition)
# 
# aggregate(value ~ Condition, data = subset(etc_analysis_long, State %in% "difference"), FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
# 
# CAT_diff <- ggplot(data = subset(etc_analysis_long[-c(466),], State %in% "percentage"), aes(x = Condition, y = value, color = Condition)) +
#   ggiraph::geom_bar_interactive(aes(x = Condition, y = value, fill = Condition, tooltip = value), position = position_dodge(0.99), stat = "summary", fill = "white", size = 1) +
#   theme_prism() +
#   scale_color_manual(values = c("grey50", "black")) +
#   coord_cartesian(ylim = c(-.7, .1), clip = "off") +
#   scale_y_continuous(expand = c(0,0), breaks = seq(-.7, .1, .1), labels = label_percent()
#   ) +
#   labs(y=expression(bold('CSC-Induced'~ANT~Inhibition))) +
#   ggiraph::geom_point_interactive(pch = 21, stroke = 1.5, show.legend = FALSE, cex = 3, fill = "white", position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0),
#                                   aes(tooltip = paste0("CSC-induced ANT Inhibition = ", scales::percent(round(value,3)), "%; Subject = ", Subject), data_id = Subject)) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank(),
#         axis.line.x = element_blank()) +
#   geom_hline(yintercept = 0, size = 1) +
#   geom_bracket(xmin = "Gastrocnemius", xmax = "Soleus", y.position = -.66, inherit.aes = F, label = "p = 0.009", size = 1, fontface = 2, tip.length = -0.02, vjust = 2.5)
# 
# ggiraph::girafe(ggobj = CAT_diff)
# 
# identify_outliers(subset(subset(etc_analysis_long, Condition %in% "Soleus"), State %in% "percentage"), value)
# 
# wilcox_test(subset(etc_analysis_long[-c(466),], State %in% "percentage"), value ~ Condition)
# 
# cohens_d(subset(etc_analysis_long, State %in% "percentage"), value ~ Condition)
# 
# aggregate(value ~ Condition, data = subset(etc_analysis_long[-c(466),], State %in% "percentage"), FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x)))
# 


