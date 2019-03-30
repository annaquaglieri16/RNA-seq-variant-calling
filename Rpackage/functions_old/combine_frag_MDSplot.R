#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args <- c("/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_red/CollectMultipleMetrics",
#  "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_red/featureCounts",
#  "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/samples_runs.txt",
#  "all_runs")

# args <- c("/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_trimmed/CollectMultipleMetrics",
#  "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/aligned_pass2_trimmed/featureCounts",
#   "/wehisan/general/user_managed/grpu_majewski_3/AML_RNA/cbf_aml_agrf/samples_runs.txt",
#   "all_runs")


input_metrics <- args[1]
input_feat_counts <- args[2]
sampledir <- args[3]
plot_id <- args[4]

library(dplyr)
library(ggplot2)
library(limma)
library(edgeR)

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

#################
#### Datasets
#################
# Datasets that I ran
samples <- read.table(args[3],stringsAsFactors = F)
colnames(samples) <- c("SID","RID")
# sample and run ID all samples
sra_path <- read.table(args[3],stringsAsFactors=F,header=F)
colnames(sra_path) <- c("SID","RID")
# Subset of the original samples including only the sample that I have aligned
insert.fin <- merge(samples,sra_path,all.x = T)

#################
#### Functions
#################

# function to combine the fragment plots for each samples from each run
read.fragdist <- function(file.dir,run){
  log.out<-read.delim(file.dir,stringsAsFactors = F,quote="",row.names = NULL)
  line.read <- which(log.out[,2]=="All_Reads.fr_count")+1
  frag.hist <- log.out[line.read:nrow(log.out),]
  frag.hist[,1] <- as.numeric(as.character(frag.hist[,1]))
  frag.hist[,2] <- as.numeric(as.character(frag.hist[,2]))
  colnames(frag.hist) <- c("insert.size","n.reads")
  frag.hist$ID <- run
  return(frag.hist)
}                                  


read.cmm.output <- function(samples_runs=insert.fin,
                            inputdir=inputdir,
                            plot_id=plot_id){
  sid <- unique(insert.fin$SID) # unique samples that I ran
  dir.create(file.path(inputdir,"combined_frag_plots"),showWarnings = FALSE,recursive = TRUE)
  
  for(i in 1:length(sid)){
    subs <- subset(samples_runs,SID%in%sid[i])
    
    pdf(file.path(inputdir,"combined_frag_plots",paste(plot_id,"_",sid[i],"_combined_frag_plots.pdf",sep="")))
   
    runs <- subs[,"RID"]
    # Read and assign cross-corr plot for each run
    list.fragdist <- list()
    for(j in 1:length(runs)){
      list.fragdist[[j]]  <- read.fragdist(file.path(inputdir,paste(runs[j],".insert_size_metrics",sep="")),
                                           run=runs[j])}
    
    combined <- do.call(rbind,list.fragdist)
    
    print(ggplot(combined) + geom_density(aes(x=insert.size,y=n.reads,fill=ID),col="white",stat="identity",alpha=0.6)+theme2() +
            ggtitle(paste("Fragment size distribution ",sid[i],sep=""))) 
                
    dev.off()
  }
  
}

## FeatureCounts
read.gene.counts <- function(file.dir){
  # Read feature counts output for every run for every sample and convert it to RPKM
  log.out<-read.delim(file.dir,stringsAsFactors = F,sep="\t",quote="",row.names = NULL,header = F)[-c(1:2),c(1,6,7)]
  log.out[,2] <- as.numeric(as.character(log.out[,2]))
  log.out[,3] <- as.numeric(as.character(log.out[,3]))
  return(log.out)
}                                  

read.fc.output <- function(samples_runs=insert.fin,
                           inputdir=input_feat_counts,
                           plot_id=plot_id){
  
  
  # Read and assign gene counts 
  list.gene.counts <- list()
  for(i in 1:nrow(samples_runs)){
    list.gene.counts[[i]]  <- read.gene.counts(file.dir=file.path(inputdir,samples_runs[i,"RID"]))}
  
  combined.gene.counts <- do.call(cbind,sapply(list.gene.counts,function(x)x[3]))

  colnames(combined.gene.counts) <- paste(samples_runs[,1],samples_runs[,2],sep="_")
  write.csv(combined.gene.counts,file.path(inputdir,paste(plot_id,"_","combined_gene_counts_data.csv",sep="")))

  pdf(file.path(inputdir,paste(plot_id,"_","mds_plot_top500genes_dge.pdf",sep="")))
  plotMDS(DGEList(combined.gene.counts),col=as.numeric(as.factor(samples_runs[,1])),xlab="Dim1",ylab="Dim2",main="MDS with top 500 gene counts")
  legend("topleft",legend=c(as.character(unique(as.factor(samples_runs[,1])))),pch=rep(16,length(unique(as.factor(samples_runs[,1])))),
    col=unique(as.numeric(as.factor(samples_runs[,1]))))
  dev.off()

  # Combine counts with info
  gene.details  <- read.delim(file.path(inputdir,samples_runs[1,"RID"]),stringsAsFactors = F,sep="\t",quote="",row.names = NULL,header = F)[-c(1:2),c(1:6)]
  combined.gene.counts.details <- cbind(gene.details,combined.gene.counts)
  colnames(combined.gene.counts.details)[1:6] <- c("Symbol","Chr","Start","End","Strand","Length")
  write.csv(combined.gene.counts.details,file.path(inputdir,paste(plot_id,"_","combined_gene_counts_data_geneIDs.csv",sep="")), row.names = FALSE)

  #### PCA
  combined_nozero <-  combined.gene.counts[-which(rowSums(combined.gene.counts)==0),]
  combined_nozero <- log2(combined_nozero+0.5)
  pcac <- prcomp(t(combined_nozero))
  res <- pcac$x
  res <- data.frame(res)
  cols <- rep(1:6,each=2)

  write.csv(res,file.path(inputdir,paste(plot_id,"_","pca_all_non_zero_genes_data.csv",sep="")))

  pdf(file.path(inputdir,paste(plot_id,"_","pca_all_non_zero_genes.pdf",sep="")))
  plot(res$PC1,res$PC2,type="n")
  text(res$PC1,res$PC2,col=cols,labels=rownames(res))
  dev.off()

 pdf(file.path(inputdir,paste(plot_id,"_","pca_all_non_zero_genes_PC1_PC3.pdf",sep="")))
  plot(res$PC1,res$PC3,type="n")
  text(res$PC1,res$PC3,col=cols,labels=rownames(res))
  dev.off()

  pdf(file.path(inputdir,paste(plot_id,"_","pca_all_non_zero_genes_PC2_PC3.pdf",sep="")))
  plot(res$PC2,res$PC3,type="n")
  text(res$PC2,res$PC3,col=cols,labels=rownames(res))
  dev.off()
  
}


read.cmm.output(samples_runs=insert.fin,
                inputdir=input_metrics,
                plot_id=plot_id)


# MDS with gene counts
read.fc.output(samples_runs=insert.fin,
               inputdir=input_feat_counts,
               plot_id=plot_id)








