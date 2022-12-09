library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)

setwd("~/projects/MPE_impute_2022/simulation/validation_output/")

# Read results of simulation impute replicates
output.list <- list()
files <- list.files(path=".", pattern="*.out", full.names=TRUE, recursive=FALSE)
i<-1
for (f in files){
  dat <- read.table(f, header=F)
  colnames(dat)<-c("File", "Impute.Type", "Miss.Prop", "Miss.Type", "Rep", "Accuracy")
  dat$File <- basename(dat$File)
  dat <- dat %>% separate(File, c("Clade.Height", "Stem.Height", "Substitution.Model", "Other", "Impute.Method"), "_")
  dat$Clade.Height <- as.numeric(as.character(gsub("c","",dat$Clade.Height)))
  dat$Stem.Height <- as.numeric(as.character(gsub("s","",dat$Stem.Height)))
  dat$Accuracy <- as.numeric(as.character(dat$Accuracy))
  dat <- dat %>% 
    select(-c(Other, Impute.Type)) %>%
    mutate("Clade.Stem.Ratio" = Clade.Height / Stem.Height)
  output.list[[i]] <- dat
  i<-i+1
}
df <- rbind.fill(output.list)

# Plot 1: Accuracy across reps 
df <- df %>% filter(Substitution.Model=="gtrgamma")
pdf("validation_accuracy_v1.pdf")
ggplot(df, aes(x=Impute.Method, y=Accuracy, group=interaction(Impute.Method, Clade.Stem.Ratio))) + 
  theme_minimal() + 
  geom_boxplot(position=position_dodge(1), aes(color=Clade.Height)) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  facet_wrap(~Miss.Type, nrow=3)
dev.off()
  
  
  
  
  
  

