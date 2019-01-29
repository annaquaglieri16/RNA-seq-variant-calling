zero_prop <- function(x){
  sum(x == 0)/length(x)
}

boxplotSample <- function(count_matrix=counts_only,RLE=TRUE,raw=FALSE){
  
  library(edgeR)
  library(dplyr)
  library(tidyr)
  
  combined_nozero_filter <- count_matrix[apply(count_matrix,1,zero_prop) < 0.7,]
  combined_nozero_filter <- as.matrix(combined_nozero_filter)
  log2CPM <- cpm(combined_nozero_filter,log = TRUE)
  
  if(RLE){
    # Compute RLE
    log2RLE <- t(apply(log2CPM,1,function(z){z - median(z)}))
    log2CPM_reshape <- data.frame(log2CPM=c(log2RLE),SampleID=rep(colnames(log2RLE),each=nrow(log2RLE)))
   log2CPM_reshape$flowcell <- do.call(c,lapply(strsplit(as.character(log2CPM_reshape$SampleName),split="_"),function(x) x[[3]] ))
  }else{
    log2CPM_reshape <- data.frame(log2CPM=c(log2CPM),SampleID=rep(colnames(log2CPM),each=nrow(log2CPM)))
    log2CPM_reshape$flowcell <- do.call(c,lapply(strsplit(as.character(log2CPM_reshape$SampleName),split="_"),function(x) x[[3]] ))
  }
    
  ave_table1 <- data.frame(log2CPM_reshape %>% group_by(SampleID) %>% 
                             
                          dplyr::summarise(Median=round(median(log2CPM),3),
                                          Mean=round(mean(log2CPM),3),
                                          FirstQuartile=round(quantile(log2CPM,c(0.25)),3),
                                          ThirdQuartile=round(quantile(log2CPM,c(0.75)),3),
                                          Max=round(max(log2CPM),3),
                                          Min=round(min(log2CPM),3),
                                          StDev=round(sd(log2CPM),3),
                                          IQR_log2=IQR(log2CPM),
                                          Min_tukey=min(log2CPM[log2CPM >= (FirstQuartile - 1.5*IQR_log2)]),
                                          Max_tukey=max(log2CPM[log2CPM <= (ThirdQuartile + 1.5*IQR_log2)])))
    
  ave_table1 <- ave_table1[order(ave_table1$SampleID),]
  ave_table1$order_sample <- rep(seq(1,nrow(ave_table1)))
  ave_table1$order_sample <- factor(ave_table1$order_sample,labels=ave_table1$SampleID)
  
  return(ave_table1)
}


viridis_pal <- function(alpha=1) {
  function(n) {
    viridis(n, alpha)
  }
}


scale_color_viridis <- function(..., alpha=1, discrete=TRUE) {
  if (discrete) {
    discrete_scale("colour", "viridis", viridis_pal(alpha), ...)
  } else {
    scale_color_gradientn(colours = viridis(256, alpha), ...)
  }
}

scale_fill_viridis <- function (..., alpha=1, discrete=TRUE) {
  if (discrete) {
    discrete_scale("fill", "viridis", viridis_pal(alpha), ...)
  } else {
    scale_fill_gradientn(colours = viridis(256, alpha), ...)
  }
}

# Function extract gene ids from pathway results
get_geneID <- function(pathway_match, list_GeneID, ncbi = ncbi, kegg = FALSE, output){
  
  x <- strsplit(as.character(pathway_match[2]), split = ",")[[1]]
  gene_id <- x[x %in% as.character(list_GeneID)]
  names(gene_id) <- ncbi$Symbol[match(gene_id,as.character(ncbi$GeneID))]
  
  if(!kegg){
    
    if(length(gene_id)==0){
      data_return <- data.frame(GeneID = "noID", 
                                Symbol ="noID", 
                                GOterm = pathway_match[1], Ont = pathway_match[3])} else {
                                  
                                  data_return <- data.frame(GeneID = gene_id, 
                                                            Symbol = names(gene_id), 
                                                            GOterm = pathway_match[1], Ont = pathway_match[3])}
  }else{
    
    data_return <- data.frame(GeneID = gene_id, 
                              Symbol = names(gene_id), KEGGpath = pathway_match[1])
  }
}

 

# Function extract gene ids from pathway results
get_geneID <- function(pathway_match, list_GeneID, ncbi = ncbi, kegg = FALSE, output){
  
  x <- strsplit(as.character(pathway_match[2]), split = ",")[[1]]
  gene_id <- x[x %in% as.character(list_GeneID)]
  names(gene_id) <- ncbi$Symbol[match(gene_id,as.character(ncbi$GeneID))]
  
  if(!kegg){
    
    if(length(gene_id)==0){
      data_return <- data.frame(GeneID = "noID", 
                                Symbol ="noID", 
                                GOterm = pathway_match[1], Ont = pathway_match[3])} else {
                                  
                                  data_return <- data.frame(GeneID = gene_id, 
                                                            Symbol = names(gene_id), 
                                                            GOterm = pathway_match[1], Ont = pathway_match[3])}
  }else{
    
    data_return <- data.frame(GeneID = gene_id, 
                              Symbol = names(gene_id), KEGGpath = pathway_match[1])
  }
}

