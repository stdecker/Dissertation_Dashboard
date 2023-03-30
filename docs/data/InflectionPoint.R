library(tidyr)
library(ggpubr)
library(rstatix)
library(ggnewscale)
library(scales)
library(ggprism)
library(flextable)
library(officer)
library(officedown)
library(segmented)
library(emmeans)
library(ggnewscale)
library(ggbeeswarm)
source("data/Aim_1_Data.R")

line_log <- data.frame(line$Muscle)

line_log <- line

line_log$Rate <- line$Rate

line_log$Concentration <- log10(line$Concentration)

line_log$Concentration[is.infinite(line_log$Concentration)] <- NA

line_log <- na.omit(line_log)


plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

line <- pivot_longer(data, c(7:19), names_to = "Smoke", values_to = "Rate")

line <- line[-c(2, 5, 6, 8:43)]

line <- na.omit(separate(line, Smoke, into = 'Concentration', sep = " ", remove = FALSE))

line$Concentration[line$Concentration == "GMDS"] <- 0

line$Concentration <- as.numeric(line$Concentration)

line$Rate_CS <- (line$Mass*line$Rate)/line$CS_activity

line$Muscle <- factor(line$Muscle, levels = c("Heart", "Soleus", "Gastrocnemius", "Aorta"))

line <- na.omit(line)

# Change lines 

abs_line <- ggplot(line, aes(x = Concentration, y = Rate, pch = Muscle, group = Muscle, color = Muscle)) +
  stat_summary(geom = "errorbar", width = 0.15, na.rm = TRUE, size = 1) +
  geom_line(aes(linetype = Muscle), stat = "summary", fun = "mean", size = .75) +
  geom_point(stat = "summary", fun = "mean", size = 2.5, stroke = .75, fill = "white") +
  scale_color_manual(values = c("green", "red", "blue", "purple")) +
  scale_shape_manual(values = c(23,21,24,22)) +
  xlab("") +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "none",
        strip.background = element_blank(),
        axis.title.y = element_text(size = 12)
  ) +
  theme_prism() +
  scale_x_continuous(breaks = c(0, 0.004, 0.04, 0.4, 4, 40, 400, 4000),
               labels = label_comma(drop0trailing = TRUE),
               trans = pseudo_log_trans(sigma = 10^-4, base = 10),
               ) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O[2]]/sec/mg[wt]))), x = "CSC Concentration (μg/mL)") +
  annotate('text', x = 25, y = 343, label = 'CSC Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 25, y = 327, label = 'Tissue Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 25, y = 312, label = 'Interaction Effect: p < 0.001', hjust = 0, fontface = 2) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 350, 50)) +
  coord_cartesian(ylim = c(0, 350))


for(y in unique(line$Muscle)){
  lineplots <- ggplot(subset(line, Muscle %in% paste(y)), aes(x = Concentration, y = Rate)) +
    geom_point(color = "black", alpha = 0.5, fill = "white") +
    geom_line(stat = "summary", fun.y = "mean", size = 1) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, color = ifelse(y == "Gastrocnemius", "blue",
                                                                                   ifelse(y == "Soleus", "red",
                                                                                          ifelse(y == "Heart", "green", "purple")))) +
    geom_point(stat = "summary", fun.y = "mean", stroke = 0.5, fill = "white", size = ifelse(y == "Heart", 4,3),
               shape = ifelse(y == "Gastrocnemius", 21,
                                  ifelse(y == "Soleus", 22,
                                      ifelse(y == "Heart", 24, 23))), color = ifelse(y == "Gastrocnemius", "blue",
                                                                                     ifelse(y == "Soleus", "red",
                                                                                            ifelse(y == "Heart", "green", "purple")))) +
    # geom_smooth(color = "black") +
    xlab("") +
    scale_x_continuous(breaks = c(0, 0.004, 0.04, 0.4, 4, 40, 400, 4000),
                       labels = label_comma(drop0trailing = TRUE),
                       trans = pseudo_log_trans(sigma = 10^-4, base = 10),
    ) +
    # scale_x_log10(breaks = c(0, 0.004, 0.04, 0.4, 4, 40, 400, 4000),
    #               labels = label_comma(drop0trailing = TRUE),
    # ) +
    theme_prism() +
    theme(axis.text.x = element_text(size = 10),
          legend.position = "none",
          strip.background = element_blank(),
          axis.title.y = element_text(size = 13)
    ) +
    labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/mg[wt]))), x = "CSC Concentration (μg/mL)")
  
  assign(paste(y, "preline", sep = "_"), lineplots)
  
}



line_ag <- aggregate(line_log$Rate, list(line_log$Concentration, line_log$Muscle), FUN = mean)

for(y in unique(line_ag$Group.2)){
  
  my.lm2 <- lm(data = subset(line_ag, Group.2 %in% paste(y)), x ~ Group.1)
  
  bp2 <- segmented(my.lm2,
                   seg.Z =  ~Group.1)
  
  assign(paste(y, "lm2", sep = "_"), my.lm2)
  
  assign(paste(y, "breaks2", sep = "_"), bp2)
  
  fit <- fitted(bp2)
  
  assign(paste(y, "fitted2", sep = "_"), fit)
  
  b0 <- coef(bp2)[[1]]
  b1 <- coef(bp2)[[2]]
  
  c1 <- coef(bp2)[[2]] + coef(bp2)[[3]]
  break1 <- bp2$psi[[2]]
  
  c0 <- b0 + b1 * break1 - c1 * break1
  
  lm_slope <- slope(bp2)
  
  assign(paste(y, "b0", sep = "_"), b0)
  assign(paste(y, "b1", sep = "_"), b1)
  
  assign(paste(y, "c0", sep = "_"), c0)
  assign(paste(y, "c1", sep = "_"), c1)
  
  assign(paste(y, "slopes", sep = '_'), lm_slope)
}




individual <- data.frame(matrix(ncol = 9, nrow = 0))

line_log$Subject <- as.character(line_log$Subject)
line_log$Subject <- factor(line_log$Subject, levels = unique(line_log$Subject))

line_log$Muscle <- factor(line_log$Muscle, levels = unique(line_log$Muscle))


for(z in unique(line_log$Subject)){
  for(y in unique(line_log$Muscle)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
  my.lm3 <- lm(data = subset(subset(line_log, Muscle %in% paste(y)), Subject %in% paste(z)), Rate ~ Concentration)
  
  summary <- summary(my.lm3)
  
  bp3 <- segmented(my.lm3,
                   seg.Z =  ~Concentration)
  
  fit <- fitted(bp3)
  
  b0 <- coef(bp3)[[1]]
  b1 <- coef(bp3)[[2]]
  
  c1 <- coef(bp3)[[2]] + coef(bp3)[[3]]
  break1 <- bp3$psi[[2]]
  
  c0 <- b0 + b1 * break1 - c1 * break1
  
  lm_slope <- slope(bp3)
  
  adjrsquared <- summary$adj.r.squared
  
  rsquared <- summary$r.squared
  
  
  a <- cbind(z, y, bp3$psi[2], b0, b1, c0, c1, rsquared, adjrsquared)
  
  individual <- rbind(individual, a)}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
  
  }
}

colnames(individual) <- c("Subject", "Tissue", 'Break', "Y1", 'Slope1', "Y2", "Slope2", "R-Squared", "Adj R-Squared")

individual$Break <- as.numeric(individual$Break)

individual$Slope2 <- as.numeric(individual$Slope2)

individual$Break <- 10^individual$Break

individual$Tissue <- factor(individual$Tissue, levels = c("Heart", "Soleus", "Gastrocnemius", "Aorta"))

individual$`Adj R-Squared` <- as.numeric(individual$`Adj R-Squared`)

individual$`R-Squared` <- as.numeric(individual$`R-Squared`)

my_comparisons2 <- list(c("Gastrocnemius", "Aorta"),
                        c("Soleus", "Aorta"),
                        c("Heart", "Aorta"))


breakpoint_pvals <- compare_means(Break~Tissue, data = individual, p.adjust.method = "holm")

slope_pvals <- compare_means(Slope2~Tissue, data = individual, p.adjust.method = "holm")

m_break <- lm(Break~Tissue, data = individual)

es_break <- eff_size(emmeans(m_break, ~Tissue), sigma = sigma(m_break), edf = df.residual(m_break))

m_slope <- lm(Slope2~Tissue, data = individual)

es_slope <- eff_size(emmeans(m_slope, ~Tissue), sigma = sigma(m_slope), edf = df.residual(m_slope))

individual_inflection <- ggplot(data = individual, aes(x = Tissue, y = Break, color = Tissue)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.3) +
  stat_summary(geom = "bar", fun = "mean", size = .75, fill = 'white') +
  theme_prism() +
  scale_color_manual(values = c("green", "red", "blue", "purple")) +
  coord_cartesian(ylim = c(0, 1000), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1000, 100)) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("[CSC] at Break (μg/mL)") +
  geom_beeswarm(cex = 3, data = individual, aes(x = Tissue, y = Break, pch = Tissue, color = Tissue),
                position = position_dodge(0.9), size = 2.5, stroke = 0.75, show.legend = FALSE, fill = "white") +
  scale_shape_manual(values = c(24, 22, 21, 23)) +
  annotate(geom = "text", x = 0.5, y = 1050, label = "Main Effect: p = 0.001", fontface = 2, hjust = 0) +
  #stat_compare_means(label.y = 870, label.x = .75, aes(label = paste0("Main Effect: p = ", ..p.format..)), fontface =2) +
  geom_bracket(xmin = "Aorta", xmax = "Gastrocnemius", y.position = 960, label = "bold(paste('p = 0.005;', ~bolditalic(d),' = 2.3'))", type = 'expression', size = 1, fontface = 2, tip.length = 0, inherit.aes = F) +
  #geom_bracket(xmin = "Aorta", xmax = "Soleus", y.position = 960, label = "bold(paste('p = 0.081;', ~bolditalic(d),' = 1.8'))", type = 'expression', size = 1, fontface = 2, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Aorta", xmax = "Heart", y.position = 880, label = "bold(paste('p = 0.001;', ~bolditalic(d),' = 2.6'))", type = 'expression', size = 1, fontface = 2, tip.length = 0, inherit.aes = F)



individual_slope2 <- ggplot(data = individual, aes(x = Tissue, y = Slope2, color = Tissue)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.3) +
  stat_summary(geom = "bar", fun = "mean", size = .75, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("green", "red", "blue", "purple")) +
  coord_cartesian(ylim = c(-300, 0), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-300, 0, 25)) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  scale_x_discrete(position = "top") +
  ylab("Slope After Break") +
  geom_beeswarm(cex = 3, data = individual, aes(x = Tissue, y = Slope2, pch = Tissue, color = Tissue),
                position = position_dodge(0.9), size = 2.5, stroke = 0.75, show.legend = FALSE, fill = 'white') +
  scale_shape_manual(values = c(24, 22, 21, 23)) +
  #stat_compare_means(label.y = -290, label.x = .75, aes(label = paste0("Main Effect: p = ", ..p.format..)), fontface =2) +
  annotate(geom = "text", x = 0.5, y = -315, label = "Main Effect: p < 0.001", fontface = 2, hjust = 0) +
  geom_bracket(xmin = "Aorta", xmax = "Soleus", y.position = -160, label = "bold(paste('p = 0.002;', ~bolditalic(d),' = 1.5'))", type = 'expression', size = 1, fontface = 2, vjust = 1.5, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Aorta", xmax = "Heart", y.position = -290, label = "bold(paste('p < 0.001;', ~bolditalic(d),' = 5.0'))", type = 'expression', size = 1, fontface = 2, vjust = 1.5, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Gastrocnemius", xmax = "Soleus", y.position = -120, label = "bold(paste('p = 0.018;', ~bolditalic(d),' = 1.3'))", type = 'expression', size = 1, fontface = 2, vjust = 1.5, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Gastrocnemius", xmax = "Heart", y.position = -260, label = "bold(paste('p < 0.001;', ~bolditalic(d),' = 4.8'))", type = 'expression', size = 1, fontface = 2, vjust = 1.5, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Soleus", xmax = "Heart", y.position = -230, label = "bold(paste('p = 0.047;', ~bolditalic(d),' = 3.5'))", type = 'expression', size = 1, fontface = 2, vjust = 1.5, tip.length = 0, inherit.aes = F)
  



individual_describe <- describeBy(individual, individual$Tissue)


#Graphs ---------------

# Gastroc ----
Gastrocnemius_line <- Gastrocnemius_preline +
  scale_y_continuous(breaks = seq(0,150,10), expand = c(0,0)) +
#  coord_cartesian(ylim = c(0, 120)) +
  coord_cartesian(ylim = c(0, 90)) +
  geom_vline(xintercept = (10^Gastrocnemius_breaks2$psi[2]), linetype = "dotted", size = .75) +
  geom_abline(intercept = Gastrocnemius_b0+2, slope = Gastrocnemius_b1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_abline(intercept = Gastrocnemius_c0*1.95, slope = Gastrocnemius_c1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_bracket(xmin = 1200, xmax = 2400, y.position = 50, label = "*", size = 1.25, fontface = 2, label.size = 8, tip.length = 0) +
  annotate(x = 0, y = 10, geom = "text", label = substitute(bold(paste('Adjusted ', R^'2', ' = ', x)), list(x = as.character(round(individual_describe$Gastrocnemius[9,3], 3)))), fontface = 2, hjust = 0)


#Soleus ----
Soleus_line <- Soleus_preline +
  scale_y_continuous(breaks = seq(0,225,25), expand = c(0,0)) +
#  coord_cartesian(ylim = c(0, 220)) +
  coord_cartesian(ylim = c(0, 225)) +
  geom_vline(xintercept = (10^Soleus_breaks2$psi[2]), linetype = "dotted", size = .75) +
  geom_abline(intercept = Soleus_b0 + 11, slope = Soleus_b1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_abline(intercept = Soleus_c0*2.025, slope = Soleus_c1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_bracket(xmin = 400, xmax = 2400, y.position = 100, label = "*", size = 1.25, fontface = 2, label.size = 8, tip.length = 0) +
  annotate(x = 0, y = 20, geom = "text", label = substitute(bold(paste('Adjusted ', R^'2', ' = ', x)), list(x = as.character(round(individual_describe$Soleus[9,3], 3)))), fontface = 2, hjust = 0)


#Aorta ----
Aorta_line <- Aorta_preline +
  scale_y_continuous(breaks = seq(0,50,5), expand = c(0,0)) +
#  coord_cartesian(ylim = c(0, 45)) +
  coord_cartesian(ylim = c(0, 45)) +
  geom_vline(xintercept = (10^Aorta_breaks2$psi[2]), linetype = "dotted", size = .75) +
  geom_abline(intercept = Aorta_b0+2.5, slope = Aorta_b1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_abline(intercept = Aorta_c0*1.975, slope = Aorta_c1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_bracket(xmin = 2000, xmax = 2400, y.position = 30, label = "*", size = 1.25, fontface = 2, label.size = 8, tip.length = 0) +
  annotate(x = 0, y = 5, geom = "text", label = substitute(bold(paste('Adjusted ', R^'2', ' = ', x)), list(x = as.character(round(individual_describe$Aorta[9,3], 3)))), fontface = 2, hjust = 0)


#Heart ----
Heart_line <- Heart_preline +
  scale_y_continuous(breaks = seq(0,450,50), expand = c(0,0)) +
#  coord_cartesian(ylim = c(0, 450)) +
  coord_cartesian(ylim = c(0, 450)) +
  geom_vline(xintercept = (10^Heart_breaks2$psi[2]), linetype = "dotted", size = .75) +
  geom_abline(intercept = Heart_b0+20, slope = Heart_b1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_abline(intercept = Heart_c0*2.025, slope = Heart_c1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_bracket(xmin = 1600, xmax = 2400, y.position = 320, label = "*", size = 1.25, fontface = 2, label.size = 8, tip.length = 0)  +
  annotate(x = 0, y = 50, geom = "text", label = substitute(bold(paste('Adjusted ', R^'2', ' = ', x)), list(x = as.character(round(individual_describe$Heart[9,3], 3)))), fontface = 2, hjust = 0)



## Per CS


CS_line <- ggplot(line, aes(x = Concentration, y = Rate_CS, pch = Muscle, group = Muscle, color = Muscle)) +
  stat_summary(geom = "errorbar", width = 0.15, na.rm = TRUE, size = 1) +
  geom_line(aes(linetype = Muscle), stat = "summary", fun = "mean", size = 1) +
  geom_point(stat = "summary", fun = "mean", size = 2, stroke = 0.5, fill = "white") +
  scale_color_manual(values = c("green", "red", "blue", "purple")) +
  scale_shape_manual(values = c(23,21,24,22)) +
  xlab("") +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "none",
        strip.background = element_blank(),
        axis.title.y = element_text(size = 12)
  ) +
  theme_prism() +
  scale_x_continuous(breaks = c(0, 0.004, 0.04, 0.4, 4, 40, 400, 4000),
                     labels = label_comma(drop0trailing = TRUE),
                     trans = pseudo_log_trans(sigma = 10^-4, base = 10),
  ) +
  labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O[2]]/sec/CS))), x = "CSC Concentration (μg/mL)") +
  annotate('text', x = 25, y = 10.85, label = 'CSC Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 25, y = 10.35, label = 'Tissue Main Effect: p < 0.001', hjust = 0, fontface = 2) +
  annotate('text', x = 25, y = 9.85, label = 'Interaction Effect: p < 0.001', hjust = 0, fontface = 2) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 11, 1)) +
  coord_cartesian(ylim = c(0, 11))

CS_individual <- data.frame(matrix(ncol = 9, nrow = 0))


for(z in unique(line_log$Subject)){
  for(y in unique(line_log$Muscle)){
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      my.lm3 <- lm(data = subset(subset(line_log, Muscle %in% paste(y)), Subject %in% paste(z)), Rate_CS ~ Concentration)
      
      summary <- summary(my.lm3)
      
      bp3 <- segmented(my.lm3,
                       seg.Z =  ~Concentration)
      
      fit <- fitted(bp3)
      
      b0 <- coef(bp3)[[1]]
      b1 <- coef(bp3)[[2]]
      
      c1 <- coef(bp3)[[2]] + coef(bp3)[[3]]
      break1 <- bp3$psi[[2]]
      
      c0 <- b0 + b1 * break1 - c1 * break1
      
      lm_slope <- slope(bp3)
      
      adjrsquared <- summary$adj.r.squared
      
      rsquared <- summary$r.squared
      
      
      a <- cbind(z, y, bp3$psi[2], b0, b1, c0, c1, rsquared, adjrsquared)
      
      CS_individual <- rbind(CS_individual, a)}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    
  }
}

colnames(CS_individual) <- c("Subject", "Tissue", 'Break', "Y1", 'Slope1', "Y2", "Slope2", "R-Squared", "Adj R-Squared")

CS_individual$Break <- as.numeric(CS_individual$Break)

CS_individual$Slope2 <- as.numeric(CS_individual$Slope2)

CS_individual$Break <- 10^CS_individual$Break

CS_individual$Tissue <- factor(CS_individual$Tissue, levels = c("Heart", "Soleus", "Gastrocnemius", "Aorta"))

CS_individual$`Adj R-Squared` <- as.numeric(CS_individual$`Adj R-Squared`)

CS_individual$`R-Squared` <- as.numeric(CS_individual$`R-Squared`)

CS_breakpoint_pvals <- compare_means(Break~Tissue, data = CS_individual, p.adjust.method = "holm")

CS_slope_pvals <- compare_means(Slope2~Tissue, data = CS_individual, p.adjust.method = "holm")

CS_m_break <- lm(Break~Tissue, data = CS_individual)

CS_es_break <- eff_size(emmeans(CS_m_break, ~Tissue), sigma = sigma(CS_m_break), edf = df.residual(CS_m_break))

CS_m_slope <- lm(Slope2~Tissue, data = CS_individual)

CS_es_slope <- eff_size(emmeans(CS_m_slope, ~Tissue), sigma = sigma(CS_m_slope), edf = df.residual(CS_m_slope))

CS_individual_inflection <- ggplot(data = CS_individual, aes(x = Tissue, y = Break, color = Tissue)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.3) +
  stat_summary(geom = "bar", fun = "mean", size = .75, fill = 'white') +
  theme_prism() +
  scale_color_manual(values = c("green", "red", "blue", "purple")) +
  coord_cartesian(ylim = c(0, 1100), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 1100, 100)) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("[CSC] at Break (g/mL)") +
  geom_beeswarm(cex = 3, data = CS_individual, aes(x = Tissue, y = Break, pch = Tissue, color = Tissue),
                position = position_dodge(0.9), size = 2.5, stroke = 0.75, show.legend = FALSE, fill = "white") +
  scale_shape_manual(values = c(24, 22, 21, 23)) +
  annotate(geom = "text", x = 0.5, y = 1080, label = "Main Effect: p = 0.001", fontface = 2, hjust = 0) +
  #stat_compare_means(label.y = 870, label.x = .75, aes(label = paste0("Main Effect: p = ", ..p.format..)), fontface =2) +
  geom_bracket(xmin = "Heart", xmax = "Soleus", y.position = 460, label = "bold(paste('p = 0.025;', ~bolditalic(d),' = 0.7'))", type = 'expression', size = 1, fontface = 2, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Heart", xmax = "Aorta", y.position = 990, label = "bold(paste('p < 0.001;', ~bolditalic(d),' = 2.6'))", type = 'expression', size = 1, fontface = 2, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Aorta", xmax = "Soleus", y.position = 930, label = "bold(paste('p = 0.026;', ~bolditalic(d),' = 1.9'))", type = 'expression', size = 1, fontface = 2, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Aorta", xmax = "Gastrocnemius", y.position = 870, label = "bold(paste('p = 0.002;', ~bolditalic(d),' = 2.3'))", type = 'expression', size = 1, fontface = 2, tip.length = 0, inherit.aes = F)



CS_individual_slope <- ggplot(data = CS_individual, aes(x = Tissue, y = Slope2, color = Tissue)) +
  #stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.3) +
  stat_summary(geom = "bar", fun = "mean", size = .75, fill = "white") +
  theme_prism() +
  scale_color_manual(values = c("green", "red", "blue", "purple")) +
  coord_cartesian(ylim = c(-8, 0), clip = "off") +
  scale_y_continuous(expand = c(0,0), breaks = seq(-8, 0, 1)) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  scale_x_discrete(position = "top") +
  ylab("Slope After Break") +
  geom_beeswarm(cex = 3, data = CS_individual, aes(x = Tissue, y = Slope2, pch = Tissue, color = Tissue),
                position = position_dodge(0.9), size = 2.5, stroke = 0.75, show.legend = FALSE, fill = 'white') +
  scale_shape_manual(values = c(24, 22, 21, 23)) +
  #stat_compare_means(label.y = -290, label.x = .75, aes(label = paste0("Main Effect: p = ", ..p.format..)), fontface =2) +
  annotate(geom = "text", x = 0.5, y = -7.9, label = "Main Effect: p < 0.001", fontface = 2, hjust = 0) +
  geom_bracket(xmin = "Aorta", xmax = "Soleus", y.position = -7, label = "bold(paste('p < 0.001;', ~bolditalic(d),' = 2.2'))", type = 'expression', size = 1, fontface = 2, vjust = 1.5, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Aorta", xmax = "Heart", y.position = -7.5, label = "bold(paste('p = 0.021;', ~bolditalic(d),' = 1.6'))", type = 'expression', size = 1, fontface = 2, vjust = 1.5, tip.length = 0, inherit.aes = F) +
  geom_bracket(xmin = "Gastrocnemius", xmax = "Aorta", y.position = -6.5, label = "bold(paste('p < 0.001;', ~bolditalic(d),' = 2.6'))", type = 'expression', size = 1, fontface = 2, vjust = 1.5, tip.length = 0, inherit.aes = F)

for(y in unique(line$Muscle)){
  lineplots <- ggplot(subset(line, Muscle %in% paste(y)), aes(x = Concentration, y = Rate_CS)) +
    geom_point(color = "black", alpha = 0.5, fill = "white") +
    geom_line(stat = "summary", fun.y = "mean", size = 1) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, color = ifelse(y == "Gastrocnemius", "blue",
                                                                                   ifelse(y == "Soleus", "red",
                                                                                          ifelse(y == "Heart", "green", "purple")))) +
    geom_point(stat = "summary", fun.y = "mean", stroke = 0.5, fill = "white", size = ifelse(y == "Heart", 4,3),
               shape = ifelse(y == "Gastrocnemius", 21,
                              ifelse(y == "Soleus", 22,
                                     ifelse(y == "Heart", 24, 23))), color = ifelse(y == "Gastrocnemius", "blue",
                                                                                  ifelse(y == "Soleus", "red",
                                                                                         ifelse(y == "Heart", "green", "purple")))) +
    # geom_smooth(color = "black") +
    xlab("") +
    scale_x_continuous(breaks = c(0, 0.004, 0.04, 0.4, 4, 40, 400, 4000),
                       labels = label_comma(drop0trailing = TRUE),
                       trans = pseudo_log_trans(sigma = 10^-4, base = 10),
    ) +
    # scale_x_log10(breaks = c(0, 0.004, 0.04, 0.4, 4, 40, 400, 4000),
    #               labels = label_comma(drop0trailing = TRUE),
    # ) +
    theme_prism() +
    theme(axis.text.x = element_text(size = 10),
          legend.position = "none",
          strip.background = element_blank(),
          axis.title.y = element_text(size = 13)
    ) +
    labs(y=expression(bold(bolditalic(J)*O['2']~(pmol[O['2']]/sec/CS))), x = "CSC Concentration (μg/mL)")
  
  assign(paste(y, "CSpreline", sep = "_"), lineplots)
  
}


CSline_ag <- aggregate(line_log$Rate_CS, list(line_log$Concentration, line_log$Muscle), FUN = mean)

for(y in unique(line_ag$Group.2)){
  
  my.lm2 <- lm(data = subset(CSline_ag, Group.2 %in% paste(y)), x ~ Group.1)
  
  bp2 <- segmented(my.lm2,
                   seg.Z =  ~Group.1)
  
  assign(paste(y, "CSlm2", sep = "_"), my.lm2)
  
  assign(paste(y, "CSbreaks2", sep = "_"), bp2)
  
  fit <- fitted(bp2)
  
  assign(paste(y, "CSfitted2", sep = "_"), fit)
  
  b0 <- coef(bp2)[[1]]
  b1 <- coef(bp2)[[2]]
  
  c1 <- coef(bp2)[[2]] + coef(bp2)[[3]]
  break1 <- bp2$psi[[2]]
  
  c0 <- b0 + b1 * break1 - c1 * break1
  
  lm_slope <- slope(bp2)
  
  assign(paste(y, "CSb0", sep = "_"), b0)
  assign(paste(y, "CSb1", sep = "_"), b1)
  
  assign(paste(y, "CSc0", sep = "_"), c0)
  assign(paste(y, "CSc1", sep = "_"), c1)
  
  assign(paste(y, "CSslopes", sep = '_'), lm_slope)
}

CSindividual_describe <- describeBy(CS_individual, CS_individual$Tissue)

# Gastroc ----
Gastrocnemius_CSline <- Gastrocnemius_CSpreline +
  scale_y_continuous(breaks = seq(0,16,2), expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 16)) +
  geom_vline(xintercept = (10^Gastrocnemius_CSbreaks2$psi[2]), linetype = "dotted", size = .75) +
  geom_abline(intercept = Gastrocnemius_CSb0 + 0.25, slope = Gastrocnemius_CSb1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_abline(intercept = Gastrocnemius_CSc0*1.95, slope = Gastrocnemius_CSc1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_bracket(xmin = 1200, xmax = 2400, y.position = 50, label = "*", size = 1.25, fontface = 2, label.size = 8, tip.length = 0) +
  annotate(x = 0, y = 1, geom = "text", label = substitute(bold(paste('Adjusted ', R^'2', ' = ', x)), list(x = as.character(round(CSindividual_describe$Gastrocnemius[9,3], 3)))), fontface = 2, hjust = 0)


#Soleus ----
Soleus_CSline <- Soleus_CSpreline +
  scale_y_continuous(breaks = seq(0,12,2), expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 12)) +
  geom_vline(xintercept = (10^Soleus_CSbreaks2$psi[2]), linetype = "dotted", size = .75) +
  geom_abline(intercept = Soleus_CSb0 + .75, slope = Soleus_CSb1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_abline(intercept = Soleus_CSc0*2.025, slope = Soleus_CSc1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_bracket(xmin = 400, xmax = 2400, y.position = 100, label = "*", size = 1.25, fontface = 2, label.size = 8, tip.length = 0) +
  annotate(x = 0, y = .85, geom = "text", label = substitute(bold(paste('Adjusted ', R^'2', ' = ', x)), list(x = as.character(round(CSindividual_describe$Soleus[9,3], 3)))), fontface = 2, hjust = 0)


#Aorta ----
Aorta_CSline <- Aorta_CSpreline +
  scale_y_continuous(breaks = seq(0,2.5,.25), expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  geom_vline(xintercept = (10^Aorta_CSbreaks2$psi[2]), linetype = "dotted", size = .75) +
  geom_abline(intercept = Aorta_CSb0 + .15, slope = Aorta_CSb1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_abline(intercept = Aorta_CSc0*1.975, slope = Aorta_CSc1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_bracket(xmin = 2000, xmax = 2400, y.position = 30, label = "*", size = 1.25, fontface = 2, label.size = 8, tip.length = 0) +
  annotate(x = 0, y = .2, geom = "text", label = substitute(bold(paste('Adjusted ', R^'2', ' = ', x)), list(x = as.character(round(CSindividual_describe$Aorta[9,3], 3)))), fontface = 2, hjust = 0)


#Heart ----
Heart_CSline <- Heart_CSpreline +
  scale_y_continuous(breaks = seq(0,10,1), expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 10)) +
  geom_vline(xintercept = (10^Heart_CSbreaks2$psi[2]), linetype = "dotted", size = .75) +
  geom_abline(intercept = Heart_CSb0 + .3, slope = Heart_CSb1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_abline(intercept = Heart_CSc0*2.04, slope = Heart_CSc1, colour = "grey25", alpha = 0.75, linetype = "longdash", size = .75) +
  geom_bracket(xmin = 1600, xmax = 2400, y.position = 320, label = "*", size = 1.25, fontface = 2, label.size = 8, tip.length = 0)  +
  annotate(x = 0, y = .75, geom = "text", label = substitute(bold(paste('Adjusted ', R^'2', ' = ', x)), list(x = as.character(round(CSindividual_describe$Heart[9,3], 3)))), fontface = 2, hjust = 0)


anova(art(Rate_CS ~ factor(line$Concentration) * factor(line$Muscle), data = line))
effectsize::eta_squared(anova(art(Rate_CS ~ factor(line$Concentration) * factor(line$Muscle), data = line)))

kruskal.test(Break ~ Tissue, data = individual)
kruskal_effsize(Break ~ Tissue, data = individual)

kruskal.test(Slope2 ~ Tissue, data = individual)
kruskal_effsize(Slope2 ~ Tissue, data = individual)

kruskal.test(Break ~ Tissue, data = CS_individual)
kruskal_effsize(Break ~ Tissue, data = CS_individual)

kruskal.test(Slope2 ~ Tissue, data = CS_individual)
kruskal_effsize(Slope2 ~ Tissue, data = CS_individual)
