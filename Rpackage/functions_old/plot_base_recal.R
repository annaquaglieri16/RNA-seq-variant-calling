## get arguments from command line
args <- commandArgs(trailingOnly=TRUE)

### Plot recalibration stuff
recal_dir <- args[1] # directory to BaseRecal folder
sample <- args[2]

library(ggplot2)


theme2<-function(angle.label=0,sizex = 10,sizey=10){
  theme(axis.text.x  = element_text(angle=angle.label, vjust=0.5, size=10),
        plot.title = element_text(size = 12, colour = "black",face="bold"),
        panel.background = element_rect(fill='white', colour='black'),
        axis.title.x = element_text(colour="grey20",size=sizex,face="bold"),
        axis.title.y = element_text(colour="grey20",size=sizey,face="bold",angle=90),
        panel.grid.major = element_line(colour = "gray95"),
        complete = TRUE
  )
}


recal <- read.csv(file.path(recal_dir,sample,paste(sample,"_recalibration_plots.csv",sep="")),stringsAsFactor=F)
recal$dif <- recal$EmpiricalQuality - recal$AverageReportedQuality

pdf(file.path(recal_dir,sample,paste(sample,"_RBQS_plots.pdf",sep="")),width=12,height=7,onefile=TRUE)
print(ggplot(recal) + geom_point(aes(x=AverageReportedQuality,y=EmpiricalQuality),size=2) + facet_wrap(~Recalibration,nrow=2)+
theme2()+geom_line(aes(x=EmpiricalQuality,y=EmpiricalQuality,linetype="dotted"))+
ggtitle("Reported Quality vs Empirical Quality"))

print(ggplot(recal) + geom_histogram(aes(x=AverageReportedQuality)) + facet_wrap(~Recalibration,nrow=2)+
theme2()+ggtitle("Reported Quality before and \n after correction"))

nucleot <- recal[recal$CovariateName == "Context",]
nucleot$len <- nchar(nucleot$CovariateValue)

print(ggplot(nucleot[nucleot$len==2,]) + geom_point(aes(x=CovariateValue,y=dif),size=2) + 
geom_hline(yintercept=0,linetype="dotted")+
facet_wrap(~Recalibration)+theme2()+
labs(y="(Empirical - Observed) Quality Score")+
ggtitle("Reported Quality by dinucleotide \n (before and after correction)"))

print(ggplot(recal[recal$CovariateName == "Cycle",]) + geom_point(aes(x=as.numeric(CovariateValue),y=as.numeric(dif))) + 
facet_wrap(~Recalibration)+theme2()+labs(y="(Empirical - Observed) Quality Score")+
ggtitle("Reported Quality by cycle \n (before and after correction)"))

dev.off()
