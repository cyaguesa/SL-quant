# ANALYZE of SL transplicing events.

#________________________________________________________________________________________________________________
# 1- FUNCTIONS
#________________________________________________________________________________________________________________

return_pos_in_operon<-function(my_op, annot=SL_quant){
  # gives you the positional order of the genes of an operon
  n=length(my_op$genes[[1]])
  res=rep("NA",n)
  if (my_op$strand=="+"){
    for (j in c(1:n)){
      res[j]=as.numeric(annot$start[which(annot$GeneID == my_op$genes[[1]][j])])
    }
    res=order(res)
  }
  else if (my_op$strand=="-"){
    for (j in c(1:n)){
      res[j]=as.numeric(annot$start[which(annot$GeneID == my_op$genes[[1]][j])])
    }
    res=order(res, decreasing=T)
  }
  return (res)
}

calculate_roc <- function(df=blu_trans,truth=13, predictor=10,cost_of_fp=1, cost_of_fn=1, n=100) {
  # Return dataframe with true positive rate and true negative rate for the predicition of thruth.
  # [truth]=column number of df containing the truth.
  # [predictor]=column number of df containing the predictor.
  tpr <- function(df, threshold,predictor, truth) {
    sum(df[,predictor] >= threshold & df[,truth]) / sum(df[,truth])
  }
  fpr <- function(df, threshold,predictor, truth) {
    sum(df[,predictor] >= threshold & !df[,truth]) / sum(!df[,truth])
  }
  cost <- function(df, threshold, cost_of_fp, cost_of_fn,predictor, truth) {
    sum(df[,predictor] >= threshold & !df[,truth]) * cost_of_fp + 
      sum(df[,predictor] < threshold & df[,truth]) * cost_of_fn
  }
  roc <- data.frame(threshold = seq(min(df[,predictor]),max(df[,predictor]),length.out=n), tpr=NA, fpr=NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th,predictor, truth))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th,predictor, truth))
  roc$cost <- sapply(roc$threshold, function(th) cost(df, th, cost_of_fp, cost_of_fn,predictor, truth))
  return(roc)
}


#________________________________________________________________________________________________________________
# 2- LOAD & pre-PROCESS data
#________________________________________________________________________________________________________________

# set directory to the home SL-quant directory
setwd("SL-quant")

# SL-quant data
SL_quant=read.table("SL-quant_results/SRR1585277_counts.txt", header=T)
SL_quant_single=read.table("SL-quant_results/SRR1585277_single_counts.txt", header=T)
SL_quant_modENCODE=read.table("SL-quant_results/modENCODE_4594_counts.txt", header=T)
SL_quant=merge(merge(SL_quant,SL_quant_single),SL_quant_modENCODE)
colnames(SL_quant)=c("GeneID","Chr","Start","End", "Strand", "Length","SL1", "SL2", "SL1_single", "SL2_single","SL1_modENCODE","SL2_modENCODE")

# prepare operon table
operon<-read.table("data/operon.gff3", sep="\t", header=F)
operon=operon[which(operon$V2=="operon"),]
operon$ID<-as.character(sapply(as.character(operon$V9),function(x) strsplit(strsplit(x, ";")[[1]][1],"Name=")[[1]][2]))
operon$genes<-as.vector(sapply(as.character(operon$V9),function(x) strsplit(strsplit(x, ';genes=')[[1]][2],',')[[1]]))
operon=operon[c(1,4,5,7,10,11)];names(operon)[1:4]=c("chrom", "start", "end", "strand")

# Manually curate operon data
operon=operon[which(sapply(operon$genes, function(g) (length(g)> 1))),] # remove operons with only one gene because they do not make sense
unlist(operon$genes)[-which((unlist(operon$genes)) %in% SL_quant$GeneID)]
operon=operon[-which(operon$ID=="CEOP4633"),] #remove CEOP4633 as it has only two genes and the first one is now considered a transposon

# order genes in operon table
genes_order=list()
for (i in c(1:length(operon$start))){
  genes_order[[i]]=(return_pos_in_operon(operon[i,]))  ;print(i)
}
operon$genes_order=genes_order
operon$genes_in_order=apply(operon, 1, function(x) x$genes[x$genes_order])

# Define the position of the genes realtive to the operons : 
# Isolated (not in operons), First_in_operon (SL1 trasnpliced), Seconds_in_operons (mostly SL2 transpliced)
SL_quant$Position="Isolated"
SL_quant$Position[which(SL_quant$GeneID %in% unlist(operon$genes))]="Seconds_in_operon"
SL_quant$Position[which(SL_quant$GeneID %in% unlist(sapply(operon$genes_in_order, function(x) x[[1]])))]="First_in_operon"
SL_quant$Position=as.factor(SL_quant$Position);summary(SL_quant$Position)

#________________________________________________________________________________________________________________
# 3- QUALITY CONTROL on blast results
#________________________________________________________________________________________________________________

# load it
blast_data=read.table("SL-quant_results/SRR1585277_blasted.txt",comment.char = "", col.names = c('query','subject','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))
blast_data_single=read.table("SL-quant_results/SRR1585277_single_blasted.txt",comment.char = "", col.names = c('query','subject','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))
blast_data_modENCODE=read.table("SL-quant_results/modENCODE_4594_blasted.txt",comment.char = "", col.names = c('query','subject','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))

#   - Length of the match (enriched at 10-11 nt length, probably due to RT drop off when approaching the cap on the SL)
plot(as.numeric(names(summary(as.factor(blast_data$length)))),summary(as.factor(blast_data$length)), type="b", ylim=c(0,185000),ylab="# blasted reads",xlab="length of the alignment (nt)", col="blue", frame=F, xlim=c(8,22))

#   - Idem comparing significant vs not significant alignments
toPlot=blast_data[which(blast_data$evalue<0.05),];plot(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), type="b", ylim=c(0,80000),ylab="# blasted reads",xlab="length of the match (nt)", col=rgb(0.9,0,0,0.8), frame=F, xlim=c(8,22), lwd=1.5)
toPlot=blast_data[which(blast_data$evalue>=0.05),];lines(as.numeric(names(summary(as.factor(toPlot$length))))[1:14],summary(as.factor(toPlot$length))[1:14], col=rgb(0,0,0.9,0.8), type="b", lwd=1.5)

#   - Idem, but only significant, comparing SL1 vs SL2
toPlot=blast_data[which(blast_data$evalue<0.05 & grepl("SL1", blast_data$subject)),];plot(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), type="b", ylim=c(0,35000),ylab="number of SL-containing reads",xlab="length of the alignment (nt)", col=rgb(0.9,0,0,0.8), frame=F, xlim=c(8,22), lwd=1.5)
toPlot=blast_data[which(blast_data$evalue<0.05 & grepl("SL2", blast_data$subject)),];lines(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), col=rgb(0,0,0.9,0.8), type="b", lwd=1.5)
legend("topright",legend=c("SL1", "SL2"),pch=c(1,1), lty=c(1,1) ,col=rep(c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)), 1),bty="n", lwd=1.5)

#   - In single mode, more enriched for 8 nt length ---> because more reads are blasted on SL sequences, we expect more aspecific matches
toPlot=blast_data_single[which(blast_data_single$evalue<0.05),];lines(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), type="b", col=rgb(0.9,0,0,0.8),lty=2, pch=3, lwd=1.5)
toPlot=blast_data_single[which(blast_data_single$evalue>0.05),];lines(as.numeric(names(summary(as.factor(toPlot$length))))[1:14],summary(as.factor(toPlot$length))[1:14], col=rgb(0,0,0.9,0.8), type="b", lty=2, pch=3, lwd=1.5)
legend("topright",legend=c("paired-end, sig.", "paired-end, NS","single-end, sig.", "single-end, NS"),pch=c(1,1,3,3), lty=c(1,1,2,2) ,col=rep(c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)), 2),bty="n", lwd=1.5)
sum(blast_data[which(blast_data$evalue<0.05),]$qstart==1)/sum(blast_data$evalue<0.05);sum(blast_data[which(blast_data$evalue>=0.05),]$qstart==1)/sum(blast_data$evalue>=0.05)
sum(blast_data_single[which(blast_data_single$evalue<0.05),]$qstart==1)/sum(blast_data_single$evalue<0.05);sum(blast_data_single[which(blast_data_single$evalue>=0.05),]$qstart==1)/sum(blast_data_single$evalue>=0.05)
sum(blast_data_modENCODE[which(blast_data_modENCODE$evalue<0.05),]$qstart==1)/sum(blast_data_modENCODE$evalue<0.05);sum(blast_data_modENCODE[which(blast_data_modENCODE$evalue>=0.05),]$qstart==1)/sum(blast_data_modENCODE$evalue>=0.05)

#   - Idem in log
toPlot=blast_data[which(blast_data$evalue<0.05),];plot(c(8:22),log2(1+c(0,0,summary(as.factor(toPlot$length))[1:13])), type="b", ylim=c(0,22),ylab="# blasted reads",xlab="length of the match (nt)", col="blue", frame=F, xlim=c(8,22))
toPlot=blast_data[which(blast_data$evalue>=0.05),];lines(as.numeric(names(summary(as.factor(toPlot$length))))[1:14],log2(1+summary(as.factor(toPlot$length)))[1:14], col="red", type="b")
toPlot=blast_data_single[which(blast_data_single$evalue<0.05),];lines(c(8:22),log2(1+c(0,0,summary(as.factor(toPlot$length))))[1:15], type="b", col="blue",lty=2, pch=3)
toPlot=blast_data_single[which(blast_data_single$evalue>0.05),];lines(as.numeric(names(summary(as.factor(toPlot$length))))[1:14],log2(1+summary(as.factor(toPlot$length)))[1:14], col="red", type="b", lty=2, pch=3)

#   - With the modENCODE dataset (true single-end data), the enrichment is sharpest at 10 nt. Perhaps because much more SL1 ! and different lib prep
plot(as.numeric(names(summary(as.factor(blast_data_modENCODE$length)))),summary(as.factor(blast_data_modENCODE$length)), type="b", ylim=c(0,185000),ylab="# blasted reads",xlab="length of the match", col="blue", frame=F, xlim=c(8,22)); abline(v=9.8, lty=2, col="darkgrey")
toPlot=blast_data_modENCODE[which(blast_data_modENCODE$evalue<0.05 & grepl("SL1", blast_data_modENCODE$subject)),];plot(c(8:22),c(0,summary(as.factor(toPlot$length))[1:14]), type="b", ylim=c(0,120000),ylab="# blasted reads",xlab="length of the match (nt)", col=rgb(0.9,0,0,0.8), frame=F, xlim=c(8,22), lwd=1.5)
toPlot=blast_data_modENCODE[which(blast_data_modENCODE$evalue<0.05 & grepl("SL2", blast_data_modENCODE$subject)),];lines(c(8:22),c(0,summary(as.factor(toPlot$length))[1:14]), col=rgb(0,0,0.9,0.8), type="b", lwd=1.5)

#________________________________________________________________________________________________________________
# 3- CAN SL-QUANT RESULTS PREDICT POSITION IN OPERONS ?
#________________________________________________________________________________________________________________

colSums(SL_quant[,c(7:12)]) # total number of SL-transplicing events : 10% Less are found in single mode for the same dataset.

# SL_quant : paired
plot(log2(SL_quant$SL1), log2(SL_quant$SL2),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,11.5), xlim=c(0,11.5))#SL1 trans-splicing events (log2)
axis(1, at=c(0,2,4,6,8,10), cex=1.5);axis(2, at=c(0,2,4,6,8,10), cex=1.5)
legend("top",legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = T)
#SL_quant[which(log2(SL_quant$SL1)>10 & SL_quant$Position=="Seconds_in_operon"),]

# SL_quant : single
plot(log2(SL_quant$SL1_single), log2(SL_quant$SL2_single),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,11), xlim=c(0,11.5))#SL1 trans-splicing events (log2)
axis(1, at=c(0,2,4,6,8,10), cex=1.5);axis(2, at=c(0,2,4,6,8,10), cex=1.5)
legend("top",legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = T)


# SL_quant : modENCODE
plot(log2(SL_quant$SL1_modENCODE), log2(SL_quant$SL2_modENCODE),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,12), xlim=c(0,12))#SL1 trans-splicing events (log2)
axis(1, at=seq(0,12,4), cex=1.5);axis(2, at=seq(0,12,4), cex=1.5)
legend(x=1, y=13,legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = F)

# ROC curves : More than 90% TPR ar 5% FDR threshold
SL_quant$SL_ratio=SL_quant$SL2/(SL_quant$SL1+SL_quant$SL2)
SL_quant$truth=SL_quant$Position=="Seconds_in_operon"
SL_quant$SL_ratio_single=SL_quant$SL2_single/(SL_quant$SL1_single+SL_quant$SL2_single)
SL_quant$SL_ratio_modENCODE=SL_quant$SL2_modENCODE/(SL_quant$SL1_modENCODE+SL_quant$SL2_modENCODE)

roc_ratio=calculate_roc(SL_quant[-which(is.na(rowSums(SL_quant[,c(16,18,19)]))),],17,16,1,1,1000)
roc_ratio_single=calculate_roc(SL_quant[-which(is.na(rowSums(SL_quant[,c(16,18,19)]))),],17,18,1,1,1000)
roc_ratio_modENCODE=calculate_roc(SL_quant[-which(is.na(rowSums(SL_quant[,c(16,18,19)]))),],17,19,1,1,1000)

plot(c(roc_ratio$fpr), c(roc_ratio$tpr), type="l", ylim=c(0.5,1), xlim=c(0,0.5), frame=F, col="blue", xlab="FPR", ylab="TPR", lwd=2)
lines(c(roc_ratio_single$fpr), c(roc_ratio_single$tpr), type="l", col="red", lty=2, lwd=2)
lines(c(roc_ratio_modENCODE$fpr), c(roc_ratio_modENCODE$tpr),lty=4, type="l", col="darkgreen", lwd=2)
abline(v=0.05, lty=2, col="darkgrey")
legend("bottomright", legend =c("SRR1585277(paired)","SRR1585277(single)","modENCODE_4594"),lty=c(1,2,4), col=c("blue", "red","darkgreen"), bty="n", lwd=2)

# Area under the curbe (AUC)
install.packages("DescTools")
library("DescTools")
AUC(c(roc_ratio$fpr), c(roc_ratio$tpr)) #0.9540
AUC(c(roc_ratio_single$fpr), c(roc_ratio_single$tpr)) #0.9536
AUC(c(roc_ratio_modENCODE$fpr), c(roc_ratio_modENCODE$tpr)) #0.9490
