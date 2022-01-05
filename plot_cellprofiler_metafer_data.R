########################
#Script to process Cell profiler data
#Maarten Paul, 20-12-2021


# Specify data sources ----------------------------------------------------

folders <- c("/media/DATA/Maarten/Fatma/211129/data/nuclei")
#folders <- c("/home/maarten/Documents/Fatma/211202-U2OS 2nd repeat/data/nuclei")
#folders <- c("/home/maarten/Documents/Fatma/211221-U2OS 3rd repeat/data/nuclei")

conditions <- c("noIR","2Gy") #these should match the file name and defined groups in Cellprofiler
celllines <- c("WT","pMP101C5","pMP101C10","pMP101C22")

replicate_add <- c(0) #per folder add a number to have unique 

edu_threshold = 100 
dna_threshold = 1000

# Initialize dependencies -------------------------------------------------

library(tidyverse)
library(ggbeeswarm)

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#ffffff"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.margin = unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#ffffff",fill="#ffffff"),
            strip.text = element_text(face="bold")
    ))
  
}
scale_fill_Publication <- function(...){
  library(scales) 
  discrete_scale("fill","Publication",manual_pal(values = c("#c00000","#1F497D","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
} 
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#6599D9","#542788","#217D68","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462","#386cb0","#fdb462")), ...)
  
}


# Import data -------------------------------------------------------------

files <- c(list.files(folders,full.names = T,include.dirs = F,pattern = "*.csv")) 
filenames <- basename(dirname(files))
replicates <- c()
if(length(replicate_add)>1){
  for (i in 1:length(replicate_add)){
    replicates <- c(replicates,rep(replicate_add[i],table(filenames)[i]))
  } 
}
metadata <- tibble(filenames,)

data <- map2(files, filenames, ~read_tsv(.x) %>% mutate(exp_name = .y))

data2 <- data[[1]]
if(length(replicate_add)>1){
  data[[1]]$Metadata_replicate...10 <- data[[1]]$Metadata_replicate...10 + replicates[1]
}
for(i in 2:length(data)){
  if(length(replicate_add)>1){
    data[[i]]$Metadata_replicate...10 <- data[[i]]$Metadata_replicate...10 + replicates[i]
  } 
  data2<- bind_rows(data2,data[[i]])
  
}#


data <- data2
rm(data2)

data <- rename(data,Metadata_cellline =Metadata_cellline...8)
data <- rename(data,Metadata_replicate =Metadata_replicate...10)
data <- rename(data,Metadata_treatment =Metadata_treatment...12)
data$Metadata_treatment <- factor(data$Metadata_treatment,levels=conditions)
data$Metadata_cellline <- factor(data$Metadata_cellline,levels=celllines)

# Summarize data -----------------------------------------------------------
data %>%
  group_by(Metadata_cellline,Metadata_treatment,Metadata_replicate)%>%
  summarize(n=n())

data %>%
  filter(Intensity_IntegratedIntensity_EdU>edu_threshold)%>%
  group_by(Metadata_cellline)%>%
  summarise(n=n())


data %>%
  filter(Intensity_IntegratedIntensity_ImageAfterMathEdU) %>%
  group_by(Metadata_cellline,Metadata_treatment)%>%
  summarise(n=n())

data %>%
  filter(Intensity_IntegratedIntensity_EdU>edu_threshold)%>%
  group_by(Metadata_cellline,Metadata_treatment)%>%
  summarise(n=n())

data %>%
  group_by(Metadata_cellline,Metadata_treatment,Metadata_replicate)%>% 
  summarise(n=n())

data %>%
  filter(Intensity_IntegratedIntensity_EdU>edu_threshold)%>%
  group_by(Metadata_cellline,Metadata_treatment,Metadata_replicate)%>%
  summarise(n=n())


# Define thresholds -------------------------------------------------------
#EdU threshold 
p <- ggplot(data,aes(x=log10(Intensity_IntegratedIntensity_ImageAfterMathEdU),color=as.character(Metadata_replicate)))+geom_histogram()+
  geom_vline(xintercept=log10(edu_threshold),color= "black",linetype="dashed",size=0.3)+
scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=10)+facet_grid(Metadata_treatment~Metadata_cellline)
p

ggsave(p,filename = file.path(folders[1],"EdU intensity profile.pdf"))

# Make plots --------------------------------------------------------------


#plot RAD51 foci before after treatment
means <- data %>%
  group_by(Metadata_cellline,Metadata_treatment,Metadata_replicate)%>%
  summarise(mean_foci=mean(Children_IdentifyPrimaryObjects_RAD51_spot_Count))

p <- data %>%
  filter(Intensity_IntegratedIntensity_ImageAfterMathEdU>edu_threshold)%>%
#  filter(Intensity_IntegratedIntensity_ImageAfterMathDNA>dna_threshold)%>%
  ggplot()+geom_boxplot(aes(y=Children_IdentifyPrimaryObjects_RAD51_spot_Count,x=Metadata_treatment,fill="red"),outlier.shape = NA,notch=T)+xlab("")+ylab("RAD51 foci in EdU+ cells")+facet_grid(.~Metadata_cellline)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+ theme(legend.position = "none")+ylim(0,100)+geom_quasirandom(data=means,aes(y=mean_foci,x=Metadata_treatment))
p
p <- data %>%
  filter(Intensity_IntegratedIntensity_ImageAfterMathEdU>edu_threshold)%>%
  ggplot()+geom_violin(aes(y=Children_IdentifyPrimaryObjects_RAD51_spot_Count,x=Metadata_treatment,fill="red"),outlier.shape = NA,notch=T)+xlab("")+ylab("RAD51 foci in EdU+ cells")+facet_grid(.~Metadata_cellline)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+ ylim(0,100)+theme(legend.position = "none")+geom_quasirandom(data=means,aes(y=mean_foci,x=Metadata_treatment))

p
ggsave(p,filename = file.path(folders[1],"fociquantification_boxplot.pdf"))

celllines <- unique(data$Metadata_cellline)
for (i in 1:length(celllines)){
  p <- data %>%
    filter(Metadata_cellline==celllines[i])%>%
    ggplot(aes(x=log10(Intensity_IntegratedIntensity_DNA),y=log10(Intensity_IntegratedIntensity_ImageAfterMathEdU),color=Children_IdentifyPrimaryObjects_RAD51_spot_Count))+geom_point(alpha=0.3)+
    facet_grid(Metadata_treatment~Metadata_replicate)+scale_color_gradient(low="grey",high="red",name="RAD51 foci")+theme_Publication(base_size=16)+xlab("DNA")+ylab("EdU")
    print(p)
    ggsave(p,filename = file.path(folders[1],paste(celllines[i],"EdU_DNA_plot.pdf")))

}

