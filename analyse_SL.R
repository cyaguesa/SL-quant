# ANALYZE of SL transplicing events.

#–––––––––––––––––––––––––––––––––––––––––––––––––––
# 1- FUNCTIONS
#–––––––––––––––––––––––––––––––––––––––––––––––––––

return_pos_in_operon<-function(my_op, annot=SL_quant){
  # gives you the positional order of the genes of an operon
  n=length(my_op$genes[[1]])
  res=rep("NA",n)
  if (my_op$strand=="+"){
    for (j in c(1:n)){
      res[j]=as.numeric(annot$Start[which(annot$GeneID == my_op$genes[[1]][j])])
    }
    res=order(res)
  }
  else if (my_op$strand=="-"){
    for (j in c(1:n)){
      res[j]=as.numeric(annot$Start[which(annot$GeneID == my_op$genes[[1]][j])])
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

log_reg <- function(df, size=10,truth=13,predictors=c(9:11)) {
  # Make a logistic regression
  # [size]
  df=df[,c(truth,predictors)]
  N <- nrow(df)
  size=10
  df <- df[sample(N),]
  num <- floor(N/size)
  rest <- N - num * size
  ncv <- cumsum(c(rep(size,num), rest))
  predictions <- data.frame(truth = df[,1], pred = NA)
  for(n in ncv) {
    v <- rep(TRUE, N)
    v[(n-size+1):n] <- FALSE
    lr <- glm(truth ~ ., data = df, family = binomial(logit))
    predictions[!v,"pred"] <- predict(lr, newdata=df[!v,], type="response")
  }
  return(predictions)
}


#–––––––––––––––––––––––––––––––––––––––––––––––––––
# 2- LOAD & pre-PROCESS data
#–––––––––––––––––––––––––––––––––––––––––––––––––––

# set directory to the home SL-quant directory
setwd("SL-quant/")

# SL-quant data
SL_quant=read.table("SL-quant_results/SRR1585277_paired/_counts.txt", header=T)
SL_quant_single=read.table("SL-quant_results/SRR1585277_single/_counts.txt", header=T)
SL_quant_paired_s=read.table("SL-quant_results/SRR1585277_paired_s/_counts.txt", header=T)
SL_quant_single_s=read.table("SL-quant_results/SRR1585277_single_s/_counts.txt", header=T)
SL_quant_single_cutadapt=read.table("SL-quant_results/SRR1585277_single_cutadapt/_counts.txt", header=T)
SL_quant_modENCODE=read.table("SL-quant_results/modENCODE_4594/_counts.txt", header=T)
SL_quant_modENCODE_s=read.table("SL-quant_results/modENCODE_4594_s/_counts.txt", header=T)
SL_quant_modENCODE_cutadapt=read.table("SL-quant_results/modENCODE_4594_cutadapt/_counts.txt", header=T)
SL_quant=cbind(SL_quant,SL_quant_single[,c(7,8)],SL_quant_paired_s[,c(7,8)],SL_quant_single_s[,c(7,8)],SL_quant_single_cutadapt[,c(7,8)],SL_quant_modENCODE[,c(7,8)],SL_quant_modENCODE_s[,c(7,8)],SL_quant_modENCODE_cutadapt[,c(7,8)])
colnames(SL_quant)=c("GeneID","Chr","Start","End", "Strand", "Length","SL1_paired", "SL2_paired", "SL1_single", "SL2_single","SL1_paired_s", "SL2_paired_s", "SL1_single_s", "SL2_single_s", "SL1_cutadapt", "SL2_cutadapt","SL1_modENCODE","SL2_modENCODE","SL1_modENCODE_s","SL2_modENCODE_s","SL1_modENCODE_cutadapt","SL2_modENCODE_cutadapt")

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

#–––––––––––––––––––––––––––––––––––––––––––––––––––
# 3- QUALITY CONTROL on blast results
#–––––––––––––––––––––––––––––––––––––––––––––––––––

# load it
blast_data=read.table("SL-quant_results/SRR1585277_paired/_blasted.txt",comment.char = "", col.names = c('query','subject','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))
blast_data_single=read.table("SL-quant_results/SRR1585277_single/_blasted.txt",comment.char = "", col.names = c('query','subject','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))
blast_data_modENCODE=read.table("SL-quant_results/modENCODE_4594/_blasted.txt",comment.char = "", col.names = c('query','subject','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))
blast_data_random=read.table("SL-quant_results/random/_blasted.txt",comment.char = "", col.names = c('query','subject','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))

#   - Length of the match (enriched at 10-11 nt length, probably due to RT drop off when approaching the cap on the SL)
plot(as.numeric(names(summary(as.factor(blast_data$length)))),summary(as.factor(blast_data$length)), type="b", ylim=c(0,185000),ylab="# blasted reads",xlab="length of the alignment (nt)", col="blue", frame=F, xlim=c(8,22))

#   - Idem comparing significant vs not significant alignments
toPlot=blast_data[which(blast_data$evalue<0.05),];plot(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), type="b", ylim=c(0,50000),ylab="number of blasted reads",xlab="length of the match (nt)", col=rgb(0.9,0,0,0.8), frame=F, xlim=c(8,22), lwd=2.5)
toPlot=blast_data[which(blast_data$evalue>=0.05),];lines(as.numeric(names(summary(as.factor(toPlot$length))))[1:14],summary(as.factor(toPlot$length))[1:14], col=rgb(0,0,0.9,0.8), type="b", lwd=2.5)
legend("topright",legend=c("sig. al. paired mode", "NS al. paired mode"),pch=c(1,1), lty=c(1,1) ,col=rep(c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)), 1),bty="n", lwd=2.5)

#   - In single mode, more enriched for 8 nt length ---> because more reads are blasted on SL sequences, we expect more aspecific matches
toPlot=blast_data_single[which(blast_data_single$evalue<0.05),];lines(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), type="b", col=rgb(0.6,0,0,0.6),lty=2, pch=3, lwd=2.5)
toPlot=blast_data_single[which(blast_data_single$evalue>0.05),];lines(as.numeric(names(summary(as.factor(toPlot$length))))[1:14],summary(as.factor(toPlot$length))[1:14], col=rgb(0,0,0.6,0.6), type="b", lty=2, pch=3, lwd=2.5)
legend("topright",legend=c("sig. al. paired mode", "NS al. paired mode", "sig. al. single mode","NS al. single mode"),pch=c(1,1), lty=c(1,1,2,2) ,col=rep(c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8),rgb(0.6,0,0,0.6),rgb(0,0,0.6,0.6)), 1),bty="n", lwd=2.5)

#   - Idem, but only significant, comparing SL1 vs SL2
toPlot=blast_data[which(blast_data$evalue<0.05 & grepl("SL1", blast_data$subject)),];plot(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), type="b", ylim=c(0,35000),ylab="number of SL-containing reads",xlab="length of the alignment (nt)", col=rgb(0.9,0,0,0.8), frame=F, xlim=c(8,22), lwd=2.5)
toPlot=blast_data[which(blast_data$evalue<0.05 & grepl("SL2", blast_data$subject)),];lines(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), col=rgb(0,0,0.9,0.8), type="b", lwd=2.5)
legend("topright",legend=c("SL1", "SL2"),pch=c(1,1), lty=c(1,1) ,col=rep(c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)), 1),bty="n", lwd=2.5)
toPlot=blast_data_single[which(blast_data_single$evalue<0.05 & grepl("SL1", blast_data_single$subject)),];lines(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), type="b", col=rgb(0.6,0,0,0.6),lty=2, pch=3, lwd=2.5)
toPlot=blast_data_single[which(blast_data_single$evalue<0.05 & grepl("SL2", blast_data_single$subject)),];lines(as.numeric(names(summary(as.factor(toPlot$length))))[1:14],summary(as.factor(toPlot$length))[1:14], col=rgb(0,0,0.6,0.6), type="b", lty=2, pch=3, lwd=2.5)
legend("topright",legend=c("SL1", "SL2","SL1 single mode", "SL2 single mode"),pch=c(1,1), lty=c(1,1,2,2) ,col=rep(c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8),rgb(0.6,0,0,0.6),rgb(0,0,0.6,0.6)), 1),bty="n", lwd=2.5)

# - clean figure for manuscript : SL1-containing reads vs SL2, SRR dataset, single-end mode
A=blast_data_single[which(blast_data_single$evalue<0.05 & grepl("SL1", blast_data_single$subject) & blast_data_single$qstart==1  & blast_data_single$send==22),]
B=blast_data_single[which(blast_data_single$evalue<0.05 & grepl("SL2", blast_data_single$subject) & blast_data_single$qstart==1  & blast_data_single$send==22),]
align_length=data.frame("reads"=c(0,0,summary(as.factor(A$length))[1:13],0,0,summary(as.factor(B$length))[1:13]))
align_length$SL=rep(c("SL1","SL2"), each=15)
align_length$length=rep(c(8:22), 2)
library(ggplot2)
ggplot(align_length, aes(x=length, y=reads, group=SL)) +
  geom_line(aes(color=SL))+ geom_point(aes(color=SL))+
  theme_minimal()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=15),legend.title=element_text(size=0))

#   - Proportion of significant/NS alignments starting at qstart==1 and ending at send==22.
sum(blast_data[which(blast_data$evalue<0.05),]$qstart==1 & blast_data[which(blast_data$evalue<0.05),]$send==22)/sum(blast_data$evalue<0.05);sum(blast_data[which(blast_data$evalue>=0.05),]$qstart==1 & blast_data[which(blast_data$evalue>=0.05),]$send==22)/sum(blast_data$evalue>=0.05)
sum(blast_data_single[which(blast_data_single$evalue<0.05),]$qstart==1 & blast_data_single[which(blast_data_single$evalue<0.05),]$send==22)/sum(blast_data_single$evalue<0.05);sum(blast_data_single[which(blast_data_single$evalue>=0.05),]$qstart==1 & blast_data_single[which(blast_data_single$evalue>=0.05),]$send==22)/sum(blast_data_single$evalue>=0.05)
sum(blast_data_modENCODE[which(blast_data_modENCODE$evalue<0.05),]$qstart==1 & blast_data_modENCODE[which(blast_data_modENCODE$evalue<0.05),]$send==22)/sum(blast_data_modENCODE$evalue<0.05);sum(blast_data_modENCODE[which(blast_data_modENCODE$evalue>=0.05),]$qstart==1 & blast_data_modENCODE[which(blast_data_modENCODE$evalue>=0.05),]$send==22)/sum(blast_data_modENCODE$evalue>=0.05)
sum(blast_data_random[which(blast_data_random$evalue<0.05),]$qstart==1 & blast_data_random[which(blast_data_random$evalue<0.05),]$send==22)/sum(blast_data_random$evalue<0.05);sum(blast_data_random[which(blast_data_random$evalue>=0.05),]$qstart==1 & blast_data_random[which(blast_data_random$evalue>=0.05),]$send==22)/sum(blast_data_random$evalue>=0.05)

proper_orientation=function(dat=blast_data){
  return(sum(dat[which(dat$evalue<0.05),]$qstart==1 & dat[which(dat$evalue<0.05),]$send==22))
}
wrong_orientation=function(dat=blast_data){
  return(sum(dat$evalue<0.05)-sum(dat[which(dat$evalue<0.05),]$qstart==1 & dat[which(dat$evalue<0.05),]$send==22))
}

orientation=data.frame("alignments"=c(proper_orientation(blast_data_modENCODE),wrong_orientation(blast_data_modENCODE),proper_orientation(blast_data_single),wrong_orientation(blast_data_single),proper_orientation(),wrong_orientation(),proper_orientation(blast_data_random),wrong_orientation(blast_data_random)))
orientation$method=factor(rep(c("modENCODE_4594", "single", "paired", "random"), each=2), levels=c("modENCODE_4594", "single", "paired","random"))
orientation$type=factor(rep(c("proper orientation", "spurious"), 4), levels=c("spurious","proper orientation"))
ggplot(orientation[which(orientation$method != "modENCODE_4594"),]) + geom_bar(aes(y =alignments, x = method, fill = type), stat="identity")+
  theme_minimal()+scale_fill_brewer(palette="Paired")+ theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=15),legend.title=element_text(size=0))


#   - With the modENCODE dataset (true single-end data), the enrichment is sharpest at 10 nt. Perhaps because much more SL1 and/or and different lib prep
plot(as.numeric(names(summary(as.factor(blast_data_modENCODE$length)))),summary(as.factor(blast_data_modENCODE$length)), type="b", ylim=c(0,185000),ylab="# blasted reads",xlab="length of the match", col="blue", frame=F, xlim=c(8,22)); abline(v=9.8, lty=2, col="darkgrey")

toPlot=blast_data_modENCODE[which(blast_data_modENCODE$evalue<0.05),];plot(c(8:22),c(0,0,summary(as.factor(toPlot$length))[1:13]), type="b", ylim=c(0,150000),ylab="number of blasted reads",xlab="length of the match (nt)", col=rgb(0.9,0,0,0.8), frame=F, xlim=c(8,22), lwd=2.5)
toPlot=blast_data_modENCODE[which(blast_data_modENCODE$evalue>=0.05),];lines(as.numeric(names(summary(as.factor(toPlot$length))))[1:14],summary(as.factor(toPlot$length))[1:14], col=rgb(0,0,0.9,0.8), type="b", lwd=2.5)
legend("topright",legend=c("sig. alignments", "NS alignments"),pch=c(1,1), lty=c(1,1) ,col=rep(c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)), 1),bty="n", lwd=2.5)

toPlot=blast_data_modENCODE[which(blast_data_modENCODE$evalue<0.05 & grepl("SL1", blast_data_modENCODE$subject)),];plot(c(8:22),c(0,summary(as.factor(toPlot$length))[1:14]), type="b", ylim=c(0,120000),ylab="# blasted reads",xlab="length of the match (nt)", col=rgb(0.9,0,0,0.8), frame=F, xlim=c(8,22), lwd=2.5)
toPlot=blast_data_modENCODE[which(blast_data_modENCODE$evalue<0.05 & grepl("SL2", blast_data_modENCODE$subject)),];lines(c(8:22),c(0,summary(as.factor(toPlot$length))[1:14]), col=rgb(0,0,0.9,0.8), type="b", lwd=2.5)
legend("topright",legend=c("SL1", "SL2"),pch=c(1,1), lty=c(1,1) ,col=rep(c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)), 1),bty="n", lwd=2.5)

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3- CAN SL-QUANT RESULTS PREDICT POSITION IN OPERONS ?
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––

colSums(SL_quant[,c(7:18)]) # total number of SL-transplicing events on genes.

# SL_quant : paired
plot(log2(SL_quant$SL1_paired), log2(SL_quant$SL2_paired),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,11.5), xlim=c(0,11.5))#SL1 trans-splicing events (log2)
axis(1, at=c(0,2,4,6,8,10), cex=1.5);axis(2, at=c(0,2,4,6,8,10), cex=1.5)
legend("top",legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = T)

# SL_quant : single
plot(log2(SL_quant$SL1_single), log2(SL_quant$SL2_single),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,11), xlim=c(0,11.5))#SL1 trans-splicing events (log2)
axis(1, at=c(0,2,4,6,8,10), cex=1.5);axis(2, at=c(0,2,4,6,8,10), cex=1.5)
legend("top",legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = T)

# SL_quant : single merged datasets
plot(log2(SL_quant$SL1_single+SL_quant$SL1_modENCODE), log2(SL_quant$SL2_single+SL_quant$SL2_modENCODE),axes=F, xlab="", ylab="",pch=19, cex=0.8, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,13), xlim=c(0,13))#SL1 trans-splicing events (log2)
axis(1, at=c(0,4,8,12), cex=1.5);axis(2, at=c(0,4,8,12), cex=1.5)
legend(0,12,legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = F)

# SL_quant : cutadapt
plot(log2(SL_quant$SL1_cutadapt), log2(SL_quant$SL2_cutadapt),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,11), xlim=c(0,11.5))#SL1 trans-splicing events (log2)
axis(1, at=c(0,2,4,6,8,10), cex=1.5);axis(2, at=c(0,2,4,6,8,10), cex=1.5)
legend("top",legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = T)

# SL_quant : single --sensitive
plot(log2(SL_quant$SL1_single_s), log2(SL_quant$SL2_single_s),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,11), xlim=c(0,11.5))#SL1 trans-splicing events (log2)
axis(1, at=c(0,2,4,6,8,10), cex=1.5);axis(2, at=c(0,2,4,6,8,10), cex=1.5)
legend("top",legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = T)

# SL_quant : single --sensitive merged datasets
plot(log2(SL_quant$SL1_single_s+SL_quant$SL1_modENCODE_s), log2(SL_quant$SL2_single_s+SL_quant$SL2_modENCODE_s),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,12), xlim=c(0,12))#SL1 trans-splicing events (log2)
axis(1, at=c(0,2,4,6,8,10), cex=1.5);axis(2, at=c(0,2,4,6,8,10), cex=1.5)
legend("top",legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = T)

# SL_quant : paired --sensitive
plot(log2(SL_quant$SL1_paired_s), log2(SL_quant$SL2_paired_s),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,11), xlim=c(0,11.5))#SL1 trans-splicing events (log2)
axis(1, at=c(0,2,4,6,8,10), cex=1.5);axis(2, at=c(0,2,4,6,8,10), cex=1.5)
legend("top",legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = T)

# SL_quant : modENCODE
plot(log2(SL_quant$SL1_modENCODE), log2(SL_quant$SL2_modENCODE),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,12.5), xlim=c(0,12))#SL1 trans-splicing events (log2)
axis(1, at=seq(0,12,4), cex=1.5);axis(2, at=seq(0,12,4), cex=1.5)
legend(x=1, y=13,legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = F)

# SL_quant : modENCODE_s
plot(log2(SL_quant$SL1_modENCODE_s), log2(SL_quant$SL2_modENCODE_s),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,12.5), xlim=c(0,12))#SL1 trans-splicing events (log2)
axis(1, at=seq(0,12,4), cex=1.5);axis(2, at=seq(0,12,4), cex=1.5)
legend(x=1, y=13,legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = F)

# SL_quant : modENCODE_cutadapt
plot(log2(SL_quant$SL1_modENCODE_cutadapt), log2(SL_quant$SL2_modENCODE_cutadapt),axes=F, xlab="", ylab="",pch=19, cex=1, col=ifelse(SL_quant$Position=="Seconds_in_operon",rgb(0.9,0,0,0.5),rgb(0,0,0.9,0.5)), frame=F, ylim=c(0,12.5), xlim=c(0,12))#SL1 trans-splicing events (log2)
axis(1, at=seq(0,12,4), cex=1.5);axis(2, at=seq(0,12,4), cex=1.5)
legend(x=1, y=13,legend=c("downstream in operon      ", "other genes"),pch=19 ,col=c(rgb(0.9,0,0,0.8),rgb(0,0,0.9,0.8)),bty="n",cex=1,seg.len=0, lwd=4, horiz = F)



# ROC curves

SL_quant$truth=SL_quant$Position=="Seconds_in_operon"
SL_quant$SL_ratio_paired=SL_quant$SL2_paired/(SL_quant$SL1_paired+SL_quant$SL2_paired)
SL_quant$SL_ratio_single=SL_quant$SL2_single/(SL_quant$SL1_single+SL_quant$SL2_single)
SL_quant$SL_ratio_paired_s=SL_quant$SL2_paired_s/(SL_quant$SL1_paired_s+SL_quant$SL2_paired_s)
SL_quant$SL_ratio_single_s=SL_quant$SL2_single_s/(SL_quant$SL1_single_s+SL_quant$SL2_single_s)
SL_quant$SL_ratio_modENCODE=SL_quant$SL2_modENCODE/(SL_quant$SL1_modENCODE+SL_quant$SL2_modENCODE)
SL_quant$SL_ratio_modENCODE_s=SL_quant$SL2_modENCODE_s/(SL_quant$SL1_modENCODE_s+SL_quant$SL2_modENCODE_s)
SL_quant$SL_ratio_single_merged=(SL_quant$SL2_single+SL_quant$SL2_modENCODE)/(SL_quant$SL1_single+SL_quant$SL2_single+SL_quant$SL1_modENCODE+SL_quant$SL2_modENCODE)
SL_quant$SL_ratio_single_s_merged=(SL_quant$SL2_single_s+SL_quant$SL2_modENCODE_s)/(SL_quant$SL1_single_s+SL_quant$SL2_single_s+SL_quant$SL2_modENCODE_s+SL_quant$SL1_modENCODE_s)
SL_quant$SL_ratio_cutadapt=(SL_quant$SL2_cutadapt)/(SL_quant$SL2_cutadapt+SL_quant$SL1_cutadapt)
SL_quant$SL_ratio_modENCODE_cutadapt=(SL_quant$SL2_modENCODE_cutadapt)/(SL_quant$SL2_modENCODE_cutadapt+SL_quant$SL1_modENCODE_cutadapt)

# paired-end dataset

considered_genes=which(!is.na(rowSums(SL_quant[,c(25:28,33)])))
roc_ratio_paired=calculate_roc(SL_quant[considered_genes,],24,25,1,1,1000)
roc_ratio_single=calculate_roc(SL_quant[considered_genes,],24,26,1,1,1000)
roc_ratio_paired_s=calculate_roc(SL_quant[considered_genes,],24,27,1,1,1000)
roc_ratio_single_s=calculate_roc(SL_quant[considered_genes,],24,28,1,1,1000)
roc_ratio_cutadapt=calculate_roc(SL_quant[considered_genes,],24,33,1,1,1000)

plot(c(roc_ratio_single$fpr), c(roc_ratio_single$tpr), type="l", ylim=c(0.5,1), xlim=c(0,0.5), frame=F, col="blue", xlab="FPR", ylab="TPR", lwd=2.5)
lines(c(roc_ratio_single_s$fpr), c(roc_ratio_single_s$tpr), type="l", col="red", lty=2, lwd=2.5)
lines(c(roc_ratio_paired$fpr), c(roc_ratio_paired$tpr), type="l", col="darkblue", lty=1, lwd=2.5)
lines(c(roc_ratio_paired_s$fpr), c(roc_ratio_paired_s$tpr), type="l", col="darkred", lty=2, lwd=2.5)
abline(v=0.05, lty=2, col="darkgrey");abline(v=0.10, lty=2, col="darkgrey");abline(h=0.90, lty=2, col="darkgrey")
legend("bottomright", legend =c(  "single","sensitive", "paired","paired, sensitive"),lty=c(1,2), col=c("blue", "red","darkblue", "darkred"), bty="n", lwd=2.5)

#paired-end dataset zoomed in
plot(c(roc_ratio_single$fpr), c(roc_ratio_single$tpr), type="l", ylim=c(0.8,1), xlim=c(0,0.2), frame=F, col="blue", xlab="FPR", ylab="TPR", lwd=2.5)
lines(c(roc_ratio_single_s$fpr), c(roc_ratio_single_s$tpr), type="l", col="red", lty=2, lwd=2.5)
lines(c(roc_ratio_paired$fpr), c(roc_ratio_paired$tpr), type="l", col="darkblue", lty=1, lwd=2.5)
lines(c(roc_ratio_paired_s$fpr), c(roc_ratio_paired_s$tpr), type="l", col="darkred", lty=2, lwd=2.5)
abline(v=0.05, lty=2, col="darkgrey");abline(v=0.10, lty=2, col="darkgrey");abline(h=0.90, lty=2, col="darkgrey")
legend(0.12,0.875, legend =c("SL-quant","SL-quant -s", "SL-quant -p","SL-quant -s -p"),lty=c(1,2), col=c("blue", "red","darkblue", "darkred"), bty="n", lwd=2.5, cex=1.2)


# single end dataset
considered_genes=which(!is.na(rowSums(SL_quant[,c(29:30,34)])))
roc_ratio_modENCODE=calculate_roc(SL_quant[considered_genes,],24,29,1,1,1000)
roc_ratio_modENCODE_s=calculate_roc(SL_quant[considered_genes,],24,30,1,1,1000)
roc_ratio_modENCODE_cutadapt=calculate_roc(SL_quant[considered_genes,],24,34,1,1,1000)

plot(c(roc_ratio_modENCODE$fpr), c(roc_ratio_modENCODE$tpr), type="l", ylim=c(0.5,1), xlim=c(0,0.5), frame=F, col="blue", xlab="FPR", ylab="TPR", lwd=2.5)
lines(c(roc_ratio_modENCODE_s$fpr), c(roc_ratio_modENCODE_s$tpr), type="l", col="red", lty=2, lwd=2.5)
abline(v=0.05, lty=2, col="darkgrey");abline(h=0.90, lty=2, col="darkgrey")
legend("bottomright", legend =c("modENCODE_4594 (normal)","modENCODE_4594 (sensitive)"),lty=c(1,2), col=c("blue", "red"), bty="n", lwd=2.5)


# Area under the curbe (AUC)
#install.packages("DescTools")
library("DescTools")
AUC(c(roc_ratio_single$fpr), c(roc_ratio_single$tpr)) # 0.9475906
AUC(c(roc_ratio_paired$fpr), c(roc_ratio_paired$tpr)) #0.9446432
AUC(c(roc_ratio_single_s$fpr), c(roc_ratio_single_s$tpr)) # 0.9578731
AUC(c(roc_ratio_paired_s$fpr), c(roc_ratio_paired_s$tpr)) #0.955828
AUC(c(roc_ratio_modENCODE$fpr), c(roc_ratio_modENCODE$tpr)) #0.9350319
AUC(c(roc_ratio_modENCODE_s$fpr), c(roc_ratio_modENCODE_s$tpr)) #0.949761


#–––––––––––––––––––––––––––––––––––
# 4- "AG" consensus in SL sites ?
#–––––––––––––––––––––––––––––––––––

# Run SL_sites.sh -> create dataframe with the results
# paired-end dataset

#paired-end dataset
AG_usage=data.frame("sites"=c(8372,9154-8372,8254,8770-8254,7957,8436-7957,6402,6539-6402,6149,6301-6149),"method"=as.factor(rep(c("Tourasse et al.", "SL-quant -s","SL-quant -s -p", "SL-quant", "SL-quant -p"), each=2)),"type"=rep(c("AG", "spurious"), 5))
AG_usage$type=as.factor(AG_usage$type);AG_usage$type=relevel(AG_usage$type,"spurious")
AG_usage$method=factor(AG_usage$method, levels=(levels(AG_usage$method))[c(5,3,4,1,2)])
ggplot(AG_usage[which(AG_usage$method != "'Ecutadapt'"),]) + geom_bar(aes(y = sites, x = method, fill = type), stat="identity")+
  theme_minimal()+scale_fill_brewer(palette="Paired")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=15),legend.title=element_text(size=0))

ggplot(AG_usage[which(AG_usage$method %in% c("Tourasse et al.", "SL-quant -s", "SL-quant")),]) + geom_bar(aes(y = sites, x = method, fill = type), stat="identity")+
  theme_minimal()+scale_fill_brewer(palette="Paired")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=15),legend.title=element_text(size=0))

#single-end dataset
AG_usage=data.frame("sites"=c(9953,11155-9953,9948,10735-9948,8081,8247-8081),"method"=as.factor(rep(c("Tourasse et al.", "SL-quant -s", "SL-quant"), each=2)),"type"=rep(c("AG", "spurious"), 3))
AG_usage$type=as.factor(AG_usage$type);AG_usage$type=relevel(AG_usage$type,"spurious")
AG_usage$method=factor(AG_usage$method, levels=rev(levels(AG_usage$method)))
ggplot(AG_usage[which(AG_usage$method != "'Ecutadapt'"),]) + geom_bar(aes(y = sites, x = method, fill = type), stat="identity")+
  theme_minimal()+scale_fill_brewer(palette="Paired")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=15),legend.title=element_text(size=0))

