
# library(cowplot)
# library(gridExtra)
# library(grid)
# library(ggplot2)
# library(dplyr)
# library(edgeR)
# library(viridis)
# library(ggthemes)
# library(gtable)

Rrepos <- "https://cloud.r-project.org"
RlibPath <- "/wehisan/home/allstaff/q/quaglieri.a/R/x86_64-pc-linux-gnu-library/3.4"

list.of.packages <- c("cowplot","gridExtra","grid","ggplot2","ggthemes","dplyr","edgeR","viridis","gtable")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = Rrepos, lib = RlibPath)

lapply(list.of.packages, require, character.only = TRUE)

args = commandArgs(trailingOnly=T)
#args <- c("/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38/star_fusion",
# "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/samples_runs1.txt",
# "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_GRCh38/featureCounts",5)


fusion_dir <- args[1]
samples <- args[2]
dir_genes <- args[3]
align_dir <- paste0(strsplit(fusion_dir,split="/")[[1]][1:(length(strsplit(fusion_dir,split="/")[[1]])-1)],collapse="/")

## Functions
# function to add a columnd SID to every data.frame in the list = list of fusions

create_key_fusion <- function(fusion_table,sample_name){
  fusion_table <- data.frame(fusion_table)
  fusion_table$SID <- sample_name
  colnames(fusion_table)[1] <- "FusionName"
  return(fusion_table)
}

## Function to get the expression of the XIST gene to identify the SEX of the sample
get_gender <- function(sample,dir_gene=dir_genes){
  genes <- read.delim(file.path(dir_gene,sample),sep="\t",quote="",row.names = NULL,header = F,stringsAsFactors=F)[-c(1:2),]
  genes$V7 <- as.numeric(genes$V7)
  genes$V6 <- as.numeric(genes$V6)
  genes.rpkm <- data.frame(genes$V1,rpkm(genes$V7,gene.length=genes$V6))
  
  xist.count <- genes.rpkm[genes.rpkm[,1]=="XIST",2]
  data <- data.frame(sample,xist.count)
  colnames(data) <- c("RID","XIST")
  return(data)
}

# Combine R fusion output
#args[2]
sam_runs <- read.table(samples,stringsAsFactors=F)
colnames(sam_runs) <- c("SID","RID")

#############################
# Determine Gender of samples
#############################
#id XIST is expressed then it's a female
# XIST is expressed on the inactive Chromosoem X
 
if(file.exists(file.path(align_dir,"XIST_Gender_estimate.csv"))){
  print(paste(file.path(align_dir,"XIST_Gender_estimate.csv")," has already been estimated"))

  xist.all<-read.csv(file.path(align_dir,"XIST_Gender_estimate.csv"))

}else{
  xist.est <- do.call(rbind,apply(sam_runs,1,function(x){
  return(get_gender(x[2],dir_genes))}))
  xist.all <- merge(xist.est,sam_runs)

  write.csv(xist.all,file.path(align_dir,"XIST_estimate.csv"))

  xist.all$round_sex <-round(xist.all$XIST,2)
  est_sex <- xist.all %>%
      group_by(SID) %>%
      summarise_each(funs(mean),round_sex)
  est_sex$gender <- ifelse(est_sex$round_sex >= 1,"Female","Male")

  write.csv(est_sex,file.path(align_dir,"XIST_Gender_estimate.csv"))
  
}
 

#############################
#############################
## Get the fusions
filter_fusion <- as.numeric(as.character(args[4]))

print("here")

# name of list argument containing fusions
fusions <- paste("fus",unique(sam_runs[,1]),sep="_")

for(i in 1:length(fusions)){
  assign(fusions[i], read.delim(file.path(fusion_dir,unique(sam_runs[,1])[i],"star-fusion.fusion_predictions.abridged.tsv"),header=TRUE,
                                stringsAsFactors=F))
  print(i)
print(unique(sam_runs[,1])[i])}

fusions_list0 <- lapply(fusions,get)
names(fusions_list0) <- unique(sam_runs[,1])

# Add SampleName to every dataset
for(i in 1:length(fusions_list0)){
  fusions_list0[[i]] <- create_key_fusion(fusions_list0[[i]],names(fusions_list0[i]))}

# Summarise read for evry FusionName without caring of different breakpoint
fusions_list <- lapply(fusions_list0,function(x){
  y <- x %>% group_by(FusionName) %>% 
  summarise(JunctionReadCount=sum(JunctionReadCount),
            SpanningFragCount=sum(SpanningFragCount),
            SID=unique(SID),
            LeftBreakpoint = min(LeftBreakpoint), 
            RightBreakpoint = max(RightBreakpoint))
  return(y)
})

# Extract complete fusion dataset and reduced
fusions_data_complete <- do.call(rbind,fusions_list0)
fusions_data <- do.call(rbind,fusions_list)
# Summary fusion
data_fus <- fusions_data %>% group_by(SID) %>%
  summarise(NFus=length(FusionName),
    MedReads=median(JunctionReadCount+SpanningFragCount),
    MaxReads=max(JunctionReadCount+SpanningFragCount),
    MinReads=min(JunctionReadCount+SpanningFragCount))

g <- tableGrob(data_fus, rows = NULL,theme=ttheme_default(base_size = 15, base_colour = "black", 
parse = FALSE, padding = unit(c(6, 6), "mm")))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))



pdf(file.path(fusion_dir,paste("combine_fusions_",filter_fusion,".pdf",sep="")),width=15,height=10)
g1=ggplot(subset(fusions_data,JunctionReadCount + SpanningFragCount >= filter_fusion),aes(x=SID,y=FusionName,group=FusionName,fill=JunctionReadCount + SpanningFragCount)) +
geom_tile(color="white", size=0.1) +scale_fill_viridis(name="# Events") + theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(y=NULL)+
  ggtitle(paste("At least ", filter_fusion," reads/fragments supporting the fusion",sep=""))
g1
dev.off()

write.csv(fusions_data,file.path(fusion_dir,"fusion_data"),row.names=F)
write.csv(fusions_data_complete,file.path(fusion_dir,"fusion_data_complete"),row.names=F)





