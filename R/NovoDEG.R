
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)
library("data.table")
dge <- "../NoveSmart_US/NovoSmart_US.nr/data/DGE/ToDP00vsToDV00.DEG.xls"

hello <- function() {
  print("I am NovoDEG.")
}

# Obtain expr matrix
nd_deg <- function(DGE_file,log2fc=0,include="a",pval=1,padj=1,desc=T){
  dge <- fread(DGE_file,sep = "\t",data.table = F)
  dge <- dge[which(dge$pval < pval & dge$padj < padj),]
  dge$lg2fc <- unlist(log2((dge[,2]+1)/(dge[,3]+1)))
  #names(gene_list) <- dge$gene_id
  dge <- dge[order(dge$lg2fc,decreasing = desc),]
  if (include == "a") {
    dge <- dge[which(abs(dge$lg2fc) > log2fc),]
  }else if (include =="u") {
    dge <- dge[which(dge$lg2fc > log2fc),]
  }else if (include == "d"){
    dge <- dge[which(dge$lg2fc < -log2fc),]
  }else{
    dge <- dge[which(abs(dge$lg2fc) > log2fc),]
  }
  return(dge)
}

if(F){
  nd_deg(dge)
  nd_deg(DGE_file = dge,log2fc = 0,include = "a",pval=1,padj=1,desc=T)
  nd_deg(DGE_file = dge,log2fc = 10,include = "d",pval = 0.05)
}

# Volcano graph
nd_volcano <- function(dge,top_n=10,gene_list=c(),fc=1,pt=0.05) {
  if (length(gene_list)){
    gene_list <- gene_list
    dge$sign <- ifelse(dge$gene_id %in% gene_list,dge$gene_id,NA)
  }else if (top_n){
    gene_list =dge[order(-log10(dge$padj),decreasing = T),]$gene_id[1:top_n]
    dge$sign <- ifelse(dge$gene_id %in% gene_list,dge$gene_id,NA)
  }else{
    dge$sign <- rep(NA,nrow(dge))  }
  dge$change <-  as.factor(ifelse(dge$padj < pt & dge$pval < pt & abs(dge$lg2fc) > fc,ifelse(dge$lg2fc > fc,'UP_reg','DOWN_reg'),'no_Sig'))

  p2 <- ggplot(dge,aes(x=lg2fc,y=-log10(padj),fill=change)) +
    geom_point(shape=21,size=2,colour=rgb(0,0,0,0.2))+
    scale_fill_manual(name = "", values = c(rgb(1,0,0,0.7), rgb(0,1,0,0.7), rgb(0.7,0.7,0.7,0.7)), limits = c('UP_regulation','DOWN_regulation','no_SigChange'))+
    geom_vline(xintercept = c(-fc,0,fc),linetype =c("dotted","dashed","dotted"),color=c("green","black","red"))+
    geom_hline(yintercept = c(-log10(pt)),linetype ="dotted")+
    geom_text_repel(aes(label = sign), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)+
    theme_bw()
  return(list(graph=p2,table=dge))
}

if(F){
  nd_volcano(nd_deg(dge,include = "a"),pt = 0.0001,fc = 7,top_n = 1)$graph
}

# Heatmap

