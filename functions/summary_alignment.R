args <- commandArgs(trailingOnly=T)

#args = c("/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/aligned_pass1",
#  "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/quaglieri.a/samples.txt")

# Directories
inputdir=args[1]
samples_id=args[2]


#install.packages("splitstackshape")
library(dplyr)
library(ggplot2)
library(reshape2)

# Functions
read.output <- function(sample,
                        inputdir=inputdir){
log.out <- read.delim(file.path(inputdir,paste(sample,"Log.final.out",sep="")),sep="\t",header = F,stringsAsFactors = F)[-c(1:4,7,22,27,31),]
return(log.out$V2)
}

theme2<-function(angle.label=0,sizex = 10,sizey=10){
  theme(axis.text.x  = element_text(angle=angle.label, vjust=0.5, size=10),
        plot.title = element_text(size = 12, colour = "black",face="bold"),
        panel.background = element_rect(fill='white', colour='gray'),
        axis.title.x = element_text(colour="grey20",size=sizex,face="bold"),
        axis.title.y = element_text(colour="grey20",size=sizey,face="bold",angle=90),
        panel.grid.major = element_line(colour = "gray95"),
        complete = TRUE
  )
}

########
#### STAR output
## <=10 loci considered multireads, otherwise they are seen as mapped to too many loci


# Read samples' Ids - Each row has the prefix of the file
samples <- read.table(samples_id,stringsAsFactors = F)
samples <- data.frame(sample=samples[!duplicated(samples),])

# Analysis
out.align <- t(apply(samples,1,function(x){
  return(read.output(x,inputdir=inputdir))
}))

names.out <- c("input.reads","ave.readlen","uniquely.mapped.num","uniquely.mapped.perc","ave.mapped.len",
               "n.splices.tot","n.splices.annot","n.splices.gt.ag","n.splices.gc.ag","n.splices.at.ac","n.splices.non.can",
               "mismatch.rate.base","del.rate.base","del.ave.len","inser.rate.base","inser.ave.len",
               "multireads","multireads.prop","mapped.too.many","mapped.too.many.prop","unmapped.mismat","unmapped.short","unmapped.other","chimeric","chimeric.prop")
colnames(out.align) <- names.out
out.align <- as.data.frame(out.align)
out.align$uniquely.mapped.perc <- substr(out.align$uniquely.mapped.perc,1,5)
out.align$mismatch.rate.base <- substr(out.align$mismatch.rate.base,1,4)
out.align$del.rate.base <- substr(out.align$del.rate.base,1,4)
out.align$inser.rate.base <- substr(out.align$inser.rate.base,1,4)
out.align$multireads.prop <- substr(out.align$multireads.prop,1,4)
out.align$mapped.too.many.prop <- substr(out.align$mapped.too.many.prop,1,4)
out.align$unmapped.mismat <- substr(out.align$unmapped.mismat,1,4)
out.align$unmapped.short <- substr(out.align$unmapped.short,1,4)
out.align$unmapped.other <- substr(out.align$unmapped.other,1,4)
out.align$chimeric.prop <- substr(out.align$chimeric.prop,1,4)

out.align$SID <- samples$sample
write.csv(out.align,file.path(inputdir,"summary.plots.data.csv"),row.names=FALSE)

# Plot 
melt.align <- melt(out.align,id.vars=c("SID"),measure.vars = c("uniquely.mapped.perc","multireads.prop","unmapped.short","chimeric.prop"))


pdf(file.path(inputdir,"summary.plots.align1.pdf"),width = 12,height = 15,onefile = T)
print(ggplot(melt.align) + geom_bar(aes(y=as.numeric(as.character(value)),x=SID,fill=variable),stat = "identity",position=position_dodge())+
  theme_bw() + coord_flip()+labs(x="Sample IDs",y="Proportion of reads") +scale_y_continuous(limits=c(0,100)))
dev.off()

pdf(file.path(inputdir,"summary.plots.align.pdf"),width = 11,height = 10,onefile = T)
for(i in names.out){
g1=ggplot(out.align) + geom_bar(aes(y=as.numeric(as.character(out.align[,i])),x=SID),stat = "identity",position=position_dodge())+
theme2(angle.label = 45) + ggtitle(i)
print(g1)
}
dev.off()



