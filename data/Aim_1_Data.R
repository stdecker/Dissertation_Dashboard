library(readxl)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(plyr)
library(scales)
library(RColorBrewer)
library(EnvStats)
library(ggprism)
library(kableExtra)
library(segmented)
library(psych)

source("data/Aim1_CS.R")
#Load data ----
data <- read_excel('C:\\Users\\Stephen\\OneDrive - University of Massachusetts\\O2MDrive_Beta\\Projects\\CS-THNR\\O2K Analysis Files\\Project 1\\Analysis File1.xlsx', sheet = "Final Data")

data[data == 0] <- NA

data <- na.omit(data)

data$Muscle[data$Muscle == "Gastroc"] <- "Gastrocnemius"

names(data)[2] <- "Age"

data <- data |> 
  mutate(CS_activity = if_else(
    Muscle == "Gastrocnemius", mean_gastroc, ifelse(
      Muscle == "Soleus", mean_soleus, ifelse(
        Muscle == "Aorta", mean_aorta, mean_heart
      )
    )
  )
)



#Make long
data_long <- pivot_longer(data, 6:19, names_to = "State", values_to = "Rate")

data_long$State <- factor(data_long$State, levels = names(data[c(6:19)]))

data_long$Muscle <- factor(data_long$Muscle, levels = unique(data_long$Muscle))

data_long <- na.omit(data_long)


# Calculate Change from CytC ----

data <- data %>% 
  mutate_at(vars(8:19), list(Change_CytC = ~-(1-(./CytC))*100))

data <- data %>% 
  mutate_at(vars(20:31), list(log_change = ~log2(.)))

data <- data %>%
  mutate_at(vars(c(7, 9:19)), list(Change_GMDS = ~-(1-(./GMDS))*100))

change_c <- pivot_longer(data, 20:31, names_to = "Concentration", values_to = "Change")
log <- pivot_longer(data, 32:43, names_to = "Concentration", values_to = "log_Change")
change <- pivot_longer(data, 44:55, names_to = "Concentration", values_to = "Change")

change <- change[-c(6:43)]
change <- na.omit(change)
change$Concentration <- factor(change$Concentration, levels = names(data[c(44:55)]))

change_c <- change_c[-c(6:43)]
change_c <- na.omit(change_c)
change_c$Concentration <- factor(change_c$Concentration, levels = names(data[c(20:31)]))

data_long$Rate_CS <- (data_long$Mass*data_long$Rate)/data_long$CS_activity
  

# Main Line graph ----
line <- pivot_longer(data, c(6, 7:19), names_to = "Smoke", values_to = "Rate")

line <- line[-c(2,5,7:42)]

line <- na.omit(separate(line, Smoke, into = 'Concentration', sep = " ", remove = FALSE))

line$Concentration[line$Concentration == "GMDS"] <- 0

line$Concentration <- as.numeric(line$Concentration)

line$Rate_CS <- (line$Mass*line$Rate)/line$CS_activity

line <- na.omit(line)

line$Concentration<- as.numeric(line$Concentration)

scientific_10 <- function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))}


library(ARTool)

art_model <- art(Rate ~ State * Muscle, data = data_long)

art_results <- anova(art_model)

art_state_con <- art.con(art_model, "State")

art_muscle_anovas <- data.frame()

art_muscle_con_pvals <- data.frame()

for(y in unique(data_long$State)){

  model <- art(Rate ~ Muscle, data = subset(data_long, State %in% paste(y)))
  
  art_muscle <- anova(model)
  
  art_muscle$State <- paste(y)
  
  art_muscle_anovas <- rbind(art_muscle_anovas, art_muscle)
  
  art_muscle_con <- summary(art.con(model, "Muscle"))
  
  art_muscle_con$State <- rep(paste(y), times = length(art_muscle_con))
  
  art_muscle_con_pvals <- rbind(art_muscle_con_pvals, art_muscle_con)
}


for(y in unique(data_long$Muscle)){
  
  model <- art(Rate ~ State, data = subset(data_long, Muscle %in% paste(y)))
  
  art_muscle_state <- anova(model)
  
  art_muscle_state_con <- summary(art.con(model, "State"))
  
  assign(paste(y, "contrasts", sep = "_"), art_muscle_state_con)
}
