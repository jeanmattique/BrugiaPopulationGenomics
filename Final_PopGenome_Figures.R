if(!require(RColorBrewer)) {
  install.packages("RColorBrewer",repos = "http://cran.us.r-project.org")
}
library("RColorBrewer")
if(!require(gplots)) {
  install.packages("gplots",repos = "http://cran.us.r-project.org")
}
library("gplots")
if(!require(png)) {
  install.packages("png",repos = "http://cran.us.r-project.org")
}
library(png)
if(!require(grid)) {
  install.packages("grid",repos = "http://cran.us.r-project.org")
}
library(grid)
if(!require(scatterplot3d)) {
  install.packages("scatterplot3d",repos = "http://cran.us.r-project.org")
}
library(scatterplot3d)
if(!require(plotly)) {
  install.packages("plotly",repos = "http://cran.us.r-project.org")
}
library(plotly)
if(!require(ggplot2)) {
  install.packages("ggplot2",repos = "http://cran.us.r-project.org")
}
library("ggplot2")
if(!require(jpeg)) {
  install.packages("jpeg",repos = "http://cran.us.r-project.org")
}
library(jpeg)
if(!require(data.table)) {
  install.packages("data.table",repos = "http://cran.us.r-project.org")
}
library(data.table)
if(!require(ggpubr)) {
  install.packages("ggpubr",repos = "http://cran.us.r-project.org")
}
library(ggpubr)
if(!require(ggsignif)) {
  install.packages("ggsignif",repos = "http://cran.us.r-project.org")
}
library(ggsignif)
#if(!require(Biostrings)) {
#	install.packages("Biostrings",repos = "http://cran.us.r-project.org")
#}
#library("Biostrings")
if(!require(ape)) {
  install.packages("ape",repos = "http://cran.us.r-project.org")
}
library("ape")
if(!require(phytools)) {
  install.packages("phytools",repos = "http://cran.us.r-project.org")
}
library("phytools")

if(!require(seqinr)) {
  install.packages("seqinr",repos = "http://cran.us.r-project.org")
}
library(seqinr)
if(!require(ggplot2)) {
  install.packages("ggplot2",repos = "http://cran.us.r-project.org")
}
library("ggplot2")
if(!require(jpeg)) {
  install.packages("jpeg",repos = "http://cran.us.r-project.org")
}
library(jpeg)
if(!require(data.table)) {
  install.packages("data.table",repos = "http://cran.us.r-project.org")
}
library(data.table)
if(!require(ggpubr)) {
  install.packages("ggpubr",repos = "http://cran.us.r-project.org")
}
library(ggpubr)
if(!require(ggsignif)) {
  install.packages("ggsignif",repos = "http://cran.us.r-project.org")
}
if(!require(ggtree)) {
  install.packages("ggtree",repos = "http://cran.us.r-project.org")
}
library(ggtree)
library("ape")
library("PopGenome")
library("pegas")
library("tidyverse")



library(ggsignif)
library(gggenes)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(scales)


BmFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BMalayi",pattern=".tsv")
BpFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BPahangi",pattern=".tsv")
WbFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/WBancrofti",pattern=".tsv")
BtFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BTimori",pattern=".tsv")

#BpFileList<-FileList[grep("Bp",FileList)]
#BmFileList<-FileList[grep("Bp",FileList,invert = T)]

Bp.res<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BPahangi.res.txt",header = FALSE,stringsAsFactors = FALSE))
Bm.res<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BMalayi.res.txt",header = FALSE,stringsAsFactors = FALSE))
Wb.res<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/WBancrofti.res.txt",header = FALSE,stringsAsFactors = FALSE))
BT.res<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/b.timori.res.txt",header = FALSE,stringsAsFactors = FALSE))

Bp.chr.res<-Bp.res[grep("Chr",Bp.res$V1),1:2]
Bm.chr.res<-Bm.res[grep("Chr",Bm.res$V1),1:2]
Wb.res<-Wb.res[,1:2]
BT.res<-BT.res[,1:2]


###Wuch Chromosome Contigs
BMRelevants<-Bm.chr.res$V1

if(1 == 2){
  Matches<-data.frame()
  pdf("/Users/jmattick/Documents/Wuch_vs_BMRef_ChrMatches_NEW.pdf")
  for (a in 1:length(BMRelevants))
  {
    
    ##Get BP Contigs Attached to BM via Promer
    sub_Bm_vs_newBP<-Bm_vs_newBP[Bm_vs_newBP$V9 == BMRelevants[a],]
    chrLength<-as.numeric(unique(sub_Bm_vs_newBP$V7))
    contigs<-unique(sub_Bm_vs_newBP$V10)
    HighLengthContigs<-data.frame(matrix(ncol=2,nrow=length(contigs)),stringsAsFactors = FALSE)
    colnames(HighLengthContigs)<-c("Contig","Length")
    for (b in 1:length(contigs))
    {
      totallen<-sum(sub_Bm_vs_newBP[sub_Bm_vs_newBP$V10 == contigs[b],]$V6)
      HighLengthContigs[b,1]<-contigs[b]
      HighLengthContigs[b,2]<-totallen
    }
    HighLengthContigs$Length<-HighLengthContigs$Length/chrLength
    HighLengthContigsFiltered<-HighLengthContigs[HighLengthContigs$Length > 0.001,]
    sumcovered<-round(sum(HighLengthContigsFiltered$Length)*100,digits=1)
    
    ##Subset Data Frames
    HighLengthContigsFiltered<-HighLengthContigs
    
    
    
    signif_sub_Bm_vs_newBP<-sub_Bm_vs_newBP[sub_Bm_vs_newBP$V10 %in% HighLengthContigsFiltered$Contig,]
    colnames(signif_sub_Bm_vs_newBP)<-c("BMStart","BMEnd","BPStart","BPEnd","BMLength","BPLength","BMSize","BPSize","BMContig","BPContig")
    if(length(signif_sub_Bm_vs_newBP$BMStart) > 0){
      print(ggplot() + geom_segment(data=signif_sub_Bm_vs_newBP,mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd)) + ggtitle(paste("Promer match for ",BMRelevants[a]," vs Wuch Contigs\nCovering ",sumcovered,"%",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[a],y="WuchContigs"))
      
      MatchesSpecific<-HighLengthContigsFiltered
      MatchesSpecific$BMChr<-BMRelevants[a]
      Matches<-rbind(Matches,MatchesSpecific)
      #+ facet_grid(V10~.,scales="free",space="free")
      #,color=BPContig
      # + facet_grid(BPContig~.,scales="free",space="free")
    }
  }
  
  dev.off()
}


####Wuch Contig IDs...
Bm_vs_Wuch <- as.data.frame(read.delim("/Users/jmattick/Documents/NucmerBmtoWuch.coords", sep="\t", header=FALSE,stringsAsFactors = F))

WuchAssignments<-data.frame(matrix(ncol=3,nrow=length(Wb.res$V1)),stringsAsFactors = FALSE)
colnames(WuchAssignments)<-c("Contig","Chr","Start")
for (a in 1:length(Wb.res$V1))
{
  subBm_vs_Wuch<-Bm_vs_Wuch[Bm_vs_Wuch$V10 == as.character(Wb.res[a,]$V1),]
  if(length(subBm_vs_Wuch$V1) > 0) {
    ChrList<-unique(subBm_vs_Wuch$V9)
    PickingFrame<-data.frame(matrix(ncol=3,nrow=length(ChrList)),stringsAsFactors = FALSE)
    colnames(PickingFrame)<-c("Chr","Sum","Start")
    for (b in 1:length(ChrList))
    {
      sumvalue<-sum(subBm_vs_Wuch[subBm_vs_Wuch$V9==ChrList[b],]$V5)
      subBm_vs_Wuch<-subBm_vs_Wuch[order(subBm_vs_Wuch$V5,decreasing = TRUE),]
      startvalue<-as.numeric(subBm_vs_Wuch[1,]$V1)
      PickingFrame[b,1]<-ChrList[b]
      PickingFrame[b,2]<-sumvalue
      PickingFrame[b,3]<-startvalue
      
    }
    PickingFrame<-PickingFrame[order(PickingFrame$Sum,decreasing = TRUE),]
    WuchAssignments[a,1]<-as.character(Wb.res[a,]$V1)
    WuchAssignments[a,2]<-PickingFrame[1,]$Chr
    WuchAssignments[a,3]<-PickingFrame[1,]$Start
    
  } else {
    WuchAssignments[a,1]<-as.character(Wb.res[a,]$V1)
    WuchAssignments[a,2]<-"None"
    WuchAssignments[a,3]<-"NA"
    
  }
}

WuchAssignments_Chrs<-WuchAssignments[grep("Chr",WuchAssignments$Chr),]


Bm_vs_Tim <- as.data.frame(read.delim("/Users/jmattick/Documents/BMvsBTCURRENT.coords", sep="\t", header=FALSE,stringsAsFactors = F))

TimAssignments<-data.frame(matrix(ncol=3,nrow=length(BT.res$V1)),stringsAsFactors = FALSE)
colnames(TimAssignments)<-c("Contig","Chr","Start")
for (a in 1:length(BT.res$V1))
{
  subBm_vs_Tim<-Bm_vs_Tim[Bm_vs_Tim$V10 == as.character(BT.res[a,]$V1),]
  if(length(subBm_vs_Tim$V1) > 0) {
    ChrList<-unique(subBm_vs_Tim$V9)
    PickingFrame<-data.frame(matrix(ncol=3,nrow=length(ChrList)),stringsAsFactors = FALSE)
    colnames(PickingFrame)<-c("Chr","Sum","Start")
    for (b in 1:length(ChrList))
    {
      sumvalue<-sum(subBm_vs_Tim[subBm_vs_Tim$V9==ChrList[b],]$V5)
      subBm_vs_Tim<-subBm_vs_Tim[order(subBm_vs_Tim$V5,decreasing = TRUE),]
      startvalue<-as.numeric(subBm_vs_Tim[1,]$V1)
      PickingFrame[b,1]<-ChrList[b]
      PickingFrame[b,2]<-sumvalue
      PickingFrame[b,3]<-startvalue
      
    }
    PickingFrame<-PickingFrame[order(PickingFrame$Sum,decreasing = TRUE),]
    TimAssignments[a,1]<-as.character(BT.res[a,]$V1)
    TimAssignments[a,2]<-PickingFrame[1,]$Chr
    TimAssignments[a,3]<-PickingFrame[1,]$Start
    
  } else {
    TimAssignments[a,1]<-as.character(BT.res[a,]$V1)
    TimAssignments[a,2]<-"None"
    TimAssignments[a,3]<-"NA"
    
  }
}

TimAssignments_Chrs<-TimAssignments[grep("Chr",TimAssignments$Chr),]



#pdf("/Users/jmattick/Documents/Bruga.malayi.pahangi.pi.pdf",width = 6,height = 15)
#pdf("/Users/jmattick/Documents/Pop.Genome.Figure2.pdf",width = 9,height = 15,useDingbats=FALSE)
pdf("/Users/jmattick/Documents/BPahangiAll.pdf",width = 9,height = 15,useDingbats=FALSE)
#pdf("/Users/jmattick/Documents/BPahangiAll_het.pdf",width = 9,height = 15,useDingbats=FALSE)


TotalMaxPreCalc<-0.02452273


###BP Plotting
#BPGroups<-c("MultiAF","Clinical","FR3","AM")
BPGroups<-c("All")
TotalMax<-c()
for (a in 1:length(BPGroups))
{
  Graphing<-data.frame()
  if(BPGroups[a] == "MultiAF"){
    subBPGroups<-BpFileList[1]
  }else if(BPGroups[a] == "Clinical"){
    subBPGroups<-BpFileList[2:4]
  }else if(BPGroups[a] == "FR3"){
    subBPGroups<-BpFileList[10:13]
  }else if(BPGroups[a] == "AM"){
    subBPGroups<-BpFileList[7:9]
  }else if(BPGroups[a] == "All"){
    subBPGroups<-BpFileList[c(1:4,7:13)]
  }
  
  TotalSNPs<-data.frame()
  SampleFactor<-length(subBPGroups)
  for (b in 1:SampleFactor)
  {
    TempSNP<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BPahangi/",subBPGroups[b],sep=""), sep="\t", header=FALSE,stringsAsFactors = FALSE))
    colnames(TempSNP)<-c("Chr","Pos","Zyg")
    TempSNP<-TempSNP[grep("Chr",TempSNP$Chr),]
    TotalSNPs<-rbind(TempSNP,TotalSNPs)
  }
  #TotalSNPs<-TotalSNPs[TotalSNPs$Zyg == "het",]
  
  ChrList<-Bm.chr.res$V1
  
  for (b in 1:length(ChrList))
  {
    ChrName<-strsplit(Bm.chr.res[b,]$V1,"_")[[1]][3]
    subChromosomes<-unique(TotalSNPs[grep(ChrName,TotalSNPs$Chr),]$Chr)
    subChromosomes<-subChromosomes[order(subChromosomes)]
    
    for (c in 1:length(subChromosomes))
    {
      
      
      ChrSize<-as.numeric(Bp.res[Bp.res$V1 == subChromosomes[c],]$V2)
      subSNPs<-TotalSNPs[TotalSNPs$Chr == subChromosomes[c],]
      
      PartialGraphing<-data.frame(matrix(ncol=4,nrow=ceiling(ChrSize/10000)),stringsAsFactors = FALSE)
      colnames(PartialGraphing)<-c("Position","SNVDensity","Chr","subChr")
      for (d in 1:ceiling(ChrSize/10000))
      {
        if(d<ceiling(ChrSize/10000)){
          SNVDensity<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/(20000*SampleFactor)
          PartialGraphing[d,]$Position<-(d*10000)
          PartialGraphing[d,]$SNVDensity<-SNVDensity
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-subChromosomes[c]
        } else {
          SNVDensity<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= ChrSize & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/((ChrSize-((d*10000)-10000))*SampleFactor)
          PartialGraphing[d,]$Position<-ChrSize
          PartialGraphing[d,]$SNVDensity<-SNVDensity
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-subChromosomes[c]
        }
      }
      Graphing<-rbind(Graphing,PartialGraphing)
    }
  }
  SNVDensity_ymax<-max(Graphing$SNVDensity)
  TotalMax<-c(TotalMax,SNVDensity_ymax)
  Graphing$Position<-Graphing$Position/1000000
  g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  
  #print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
  gl<-list(g1,g2,g3,g4,g5)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
  BPahangiGraphing<-Graphing
  
}
dev.off()



g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("Chr1") + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("Chr2") + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("Chr3") + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("Chr4") + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("ChrX") + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))


pdf("/Users/jmattick/Documents/BPChr1.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g1)
dev.off()
pdf("/Users/jmattick/Documents/BPChr2.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g2)
dev.off()
pdf("/Users/jmattick/Documents/BPChr3.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g3)
dev.off()
pdf("/Users/jmattick/Documents/BPChr4.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g4)
dev.off()
#pdf("/Users/jmattick/Documents/BPChrX.pdf",width = 20,height = 10,useDingbats=FALSE)
pdf("/Users/jmattick/Documents/BPChrX_het.pdf",width = 20,height = 10,useDingbats=FALSE)

print(g5)
dev.off()

###BM Plotting
#BMGroups<-c("Lucknow","Malaysia","FR3","TRS","WashU","Liverpool")
#pdf("/Users/jmattick/Documents/Pop.Genome.Figure1.pdf",width = 9,height = 15,useDingbats=FALSE)
pdf("/Users/jmattick/Documents/BMalayi_all.pdf",width = 9,height = 15,useDingbats=FALSE)
#pdf("/Users/jmattick/Documents/BMalayi_subsetted.pdf",width = 9,height = 15,useDingbats=FALSE)

BMGroups<-c("All")
#BMGroups<-c("Lucknow","Malaysia","FR3","TRS","WashU","Liverpool")

#BMGroups<-c("Malaysia","Malaysia_Masked")
TotalMaxPreCalc<-0.02
TotalMax<-0

for (a in 1:length(BMGroups))
{
  Graphing<-data.frame()
  if(BMGroups[a] == "Lucknow"){
    subBMGroups<-BmFileList[1:4]
  }else if(BMGroups[a] == "Malaysia"){
    subBMGroups<-BmFileList[5:8]
  }else if(BMGroups[a] == "FR3"){
    subBMGroups<-BmFileList[9:12]
  }else if(BMGroups[a] == "TRS"){
    subBMGroups<-BmFileList[c(13:16)]
  }else if(BMGroups[a] == "WashU"){
    subBMGroups<-BmFileList[21:26]
  }else if(BMGroups[a] == "Liverpool"){
    subBMGroups<-BmFileList[17:20]
  }else if(BMGroups[a] == "All"){
    subBMGroups<-BmFileList
  }
  
  TotalSNPs<-data.frame()
  SampleFactor<-length(subBMGroups)
  for (b in 1:SampleFactor)
  {
    TempSNP<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BMalayi/",subBMGroups[b],sep=""), sep="\t", header=FALSE,stringsAsFactors = FALSE))
    colnames(TempSNP)<-c("Chr","Pos","Zyg")
    TempSNP<-TempSNP[grep("Chr",TempSNP$Chr),]
    TotalSNPs<-rbind(TempSNP,TotalSNPs)
  }
  
  
  ChrList<-Bm.chr.res$V1
  
  for (b in 1:length(ChrList))
  {
    ChrName<-strsplit(Bm.chr.res[b,]$V1,"_")[[1]][3]
    subChromosomes<-unique(TotalSNPs[grep(ChrName,TotalSNPs$Chr),]$Chr)
    subChromosomes<-subChromosomes[order(subChromosomes)]
    
    for (c in 1:length(subChromosomes))
    {
      
      
      ChrSize<-as.numeric(Bm.res[Bm.res$V1 == subChromosomes[c],]$V2)
      subSNPs<-TotalSNPs[TotalSNPs$Chr == subChromosomes[c],]
      
      PartialGraphing<-data.frame(matrix(ncol=4,nrow=ceiling(ChrSize/10000)),stringsAsFactors = FALSE)
      colnames(PartialGraphing)<-c("Position","SNVDensity","Chr","subChr")
      for (d in 1:ceiling(ChrSize/10000))
      {
        if(d<ceiling(ChrSize/10000)){
          SNVDensity<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/(20000*SampleFactor)
          PartialGraphing[d,]$Position<-(d*10000)
          PartialGraphing[d,]$SNVDensity<-SNVDensity
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-subChromosomes[c]
        } else {
          SNVDensity<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= ChrSize & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/((ChrSize-((d*10000)-10000))*SampleFactor)
          PartialGraphing[d,]$Position<-ChrSize
          PartialGraphing[d,]$SNVDensity<-SNVDensity
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-subChromosomes[c]
        }
      }
      Graphing<-rbind(Graphing,PartialGraphing)
    }
  }
  SNVDensity_ymax<-max(Graphing$SNVDensity)
  TotalMax<-max(c(TotalMax,SNVDensity_ymax))
  Graphing$Position<-Graphing$Position/1000000
  #g6<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  #g7<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  #g8<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  #g9<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  #g10<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  
  
  g6<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":A",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g7<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":B",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g8<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":C",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g9<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":D",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g10<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=SNVDensity),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":E",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  
  #print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
  gl<-list(g6,g7,g8,g9,g10)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
  BMalayiGraphing<-Graphing
  #BMalayiGraphingStorage<-BMalayiGraphing
  
}
dev.off()

g6<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":Chr1",sep="")) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
g7<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":Chr2",sep="")) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
g8<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":Chr3",sep="")) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
g9<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":Chr4",sep="")) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
g10<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=Pi),size=1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":ChrX",sep="")) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))


pdf("/Users/jmattick/Documents/BMChr1.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g6)
dev.off()
pdf("/Users/jmattick/Documents/BMChr2.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g7)
dev.off()
pdf("/Users/jmattick/Documents/BMChr3.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g8)
dev.off()
pdf("/Users/jmattick/Documents/BMChr4.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g9)
dev.off()
pdf("/Users/jmattick/Documents/BMChrX.pdf",width = 20,height = 10,useDingbats=FALSE)
print(g10)
dev.off()

#BMalayiGraphing<-BMalayiGraphingStorage

###Wuchereria Plotting
#pdf("/Users/jmattick/Documents/WBancroftiPlottingAll.pdf",width = 9,height = 15,useDingbats=FALSE)
pdf("/Users/jmattick/Documents/WBancroftiPlotting_lowmax.pdf",width = 20,height = 10,useDingbats=FALSE)


#TotalMaxPreCalc<-0.075
TotalMaxPreCalc<-0.025


###WB Plotting
#WBGroups<-c("Haiti","Kenya","Mali","PNG")
WBGroups<-c("All")
TotalMax<-c()
for (a in 1:length(WBGroups))
{
  Graphing<-data.frame()
  if(WBGroups[a] == "Haiti"){
    subWBGroups<-WbFileList[1:4]
  }else if(WBGroups[a] == "Kenya"){
    subWBGroups<-WbFileList[5:8]
  }else if(WBGroups[a] == "Mali"){
    subWBGroups<-WbFileList[9:12]
  }else if(WBGroups[a] == "PNG"){
    subWBGroups<-WbFileList[13:16]
  }else if(WBGroups[a] == "All"){
    subWBGroups<-WbFileList
  }
  
  TotalSNPs<-data.frame()
  SampleFactor<-length(subWBGroups)
  for (b in 1:SampleFactor)
  {
    TempSNP<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/WBancrofti/",subWBGroups[b],sep=""), sep="\t", header=FALSE,stringsAsFactors = FALSE))
    colnames(TempSNP)<-c("Chr","Pos","Zyg")
    TotalSNPs<-rbind(TempSNP,TotalSNPs)
  }
  
  
  ChrList<-Bm.chr.res$V1
  
  for (b in 1:length(ChrList))
  {
    ChrName<-strsplit(Bm.chr.res[b,]$V1,"_")[[1]][3]
    subContigs<-WuchAssignments_Chrs[grep(ChrName,WuchAssignments_Chrs$Chr),]
    subContigs<-subContigs[order(subContigs$Start),]
    sumpos<-0
    for (c in 1:length(subContigs$Contig))
    {
      
      
      ChrSize<-as.numeric(Wb.res[Wb.res$V1 == subContigs[c,]$Contig,]$V2)
      subSNPs<-TotalSNPs[TotalSNPs$Chr == subContigs[c,]$Contig,]
      
      PartialGraphing<-data.frame(matrix(ncol=5,nrow=ceiling(ChrSize/10000)),stringsAsFactors = FALSE)
      colnames(PartialGraphing)<-c("Position","SumPosition","Pi","Chr","subChr")
      for (d in 1:ceiling(ChrSize/10000))
      {
        if(d<ceiling(ChrSize/10000)){
          pi<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/(20000*SampleFactor)
          PartialGraphing[d,]$Position<-(d*10000)
          PartialGraphing[d,]$SumPosition<-(d*10000)+sumpos
          PartialGraphing[d,]$Pi<-pi
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-as.character(subContigs[c,]$Contig)
        } else {
          pi<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= ChrSize & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/((ChrSize-((d*10000)-10000))*SampleFactor)
          PartialGraphing[d,]$Position<-ChrSize
          PartialGraphing[d,]$SumPosition<-ChrSize+sumpos
          PartialGraphing[d,]$Pi<-pi
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-as.character(subContigs[c,]$Contig)
        }
      }
      sumpos<-sumpos+ChrSize
      Graphing<-rbind(Graphing,PartialGraphing)
      if(length(Graphing[is.na(Graphing$Pi),]$Pi) > 0){
        stop()
      }
      
    }
  }
  pi_ymax<-max(Graphing$Pi)
  TotalMax<-c(TotalMax,pi_ymax)
  Graphing$Position<-Graphing$Position/1000000
  Graphing$SumPosition<-Graphing$SumPosition/1000000
  
  g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  
  #print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
  gl<-list(g1,g2,g3,g4,g5)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
  WBancroftiGraphing<-Graphing
  
}
dev.off()


pdf("/Users/jmattick/Documents/BTimoriPlotting.pdf",width = 20,height = 10,useDingbats=FALSE)

TotalMaxPreCalc<-0.025

TotalMax<-c()
for (a in 1:length(WBGroups))
{
  Graphing<-data.frame()
  
  subBTGroups<-BtFileList
  
  
  TotalSNPs<-data.frame()
  SampleFactor<-length(subBTGroups)
  for (b in 1:SampleFactor)
  {
    TempSNP<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BTimori/",subBTGroups[b],sep=""), sep="\t", header=FALSE,stringsAsFactors = FALSE))
    colnames(TempSNP)<-c("Chr","Pos","Zyg")
    TotalSNPs<-rbind(TempSNP,TotalSNPs)
  }
  
  
  ChrList<-Bm.chr.res$V1
  
  for (b in 1:length(ChrList))
  {
    ChrName<-strsplit(Bm.chr.res[b,]$V1,"_")[[1]][3]
    subContigs<-TimAssignments_Chrs[grep(ChrName,TimAssignments_Chrs$Chr),]
    subContigs<-subContigs[order(subContigs$Start),]
    sumpos<-0
    for (c in 1:length(subContigs$Contig))
    {
      
      
      ChrSize<-as.numeric(BT.res[BT.res$V1 == subContigs[c,]$Contig,]$V2)
      subSNPs<-TotalSNPs[TotalSNPs$Chr == subContigs[c,]$Contig,]
      
      PartialGraphing<-data.frame(matrix(ncol=5,nrow=ceiling(ChrSize/10000)),stringsAsFactors = FALSE)
      colnames(PartialGraphing)<-c("Position","SumPosition","Pi","Chr","subChr")
      for (d in 1:ceiling(ChrSize/10000))
      {
        if(d<ceiling(ChrSize/10000)){
          pi<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/(20000*SampleFactor)
          PartialGraphing[d,]$Position<-(d*10000)
          PartialGraphing[d,]$SumPosition<-(d*10000)+sumpos
          PartialGraphing[d,]$Pi<-pi
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-as.character(subContigs[c,]$Contig)
        } else {
          pi<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= ChrSize & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/((ChrSize-((d*10000)-10000))*SampleFactor)
          PartialGraphing[d,]$Position<-ChrSize
          PartialGraphing[d,]$SumPosition<-ChrSize+sumpos
          PartialGraphing[d,]$Pi<-pi
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-as.character(subContigs[c,]$Contig)
        }
      }
      sumpos<-sumpos+ChrSize
      Graphing<-rbind(Graphing,PartialGraphing)
      if(length(Graphing[is.na(Graphing$Pi),]$Pi) > 0){
        stop()
      }
      
    }
  }
  pi_ymax<-max(Graphing$Pi)
  TotalMax<-c(TotalMax,pi_ymax)
  TotalMaxPreCalc<-pi_ymax
  Graphing$Position<-Graphing$Position/1000000
  Graphing$SumPosition<-Graphing$SumPosition/1000000
  
  g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=SumPosition,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  
  #print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
  gl<-list(g1,g2,g3,g4,g5)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
  BTimoriGraphing<-Graphing
  
}
dev.off()




####G-Pho Sequence File Sorting
BMalayiGraphing_GPHO<-BMalayiGraphing
BMalayiGraphing_GPHO$Position<-BMalayiGraphing_GPHO$Position*1000000
BMalayiGraphing_GPHO<-BMalayiGraphing_GPHO[order(-BMalayiGraphing_GPHO$Pi),]
BMalayiGraphing_GPHO<-BMalayiGraphing_GPHO[BMalayiGraphing_GPHO$Chr != "ChrX",]
BMalayiGraphing_GPHO<-BMalayiGraphing_GPHO[1:1000,]

GPho_Coords<-data.frame(BMalayiGraphing_GPHO$Chr,(BMalayiGraphing_GPHO$Position-9999),BMalayiGraphing_GPHO$Position)
colnames(GPho_Coords)<-c("Chr","Start","Stop")
GPho_Coords$Start <- format(GPho_Coords$Start, scientific = FALSE)
GPho_Coords$Stop <- format(GPho_Coords$Stop, scientific = FALSE)


write.table(GPho_Coords,"/Users/jmattick/Documents/GPhoCoords.bed",row.names = FALSE,col.names = FALSE,sep = "\t",quote = FALSE)

###PCA Plots
##BP

Bp.vec<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/PCA/BP.plink2.eigenvec",header = TRUE,stringsAsFactors = FALSE))
Bp.val<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/PCA/BP.plink2.eigenval",header = FALSE,stringsAsFactors = FALSE))

Bp.vec.subset<-Bp.vec[,1:4]
Bp.vec.subset$SampleGroup<-c("MultiAF",rep("FR3",3),rep("Clinical",3),rep("FR3",4))
#print(ggplot(gene345) + geom_point(aes(x=miRNA,y=mRNA,color=Sample,size=5)) + theme_bw() + ggtitle("WB_new_34 miRNA Depth vs mRNA depth") + theme(plot.title = element_text(size=14)) + xlim(0,7000) + ylim(0,2000) + scale_color_manual(values = colors)+ theme(legend.text=element_text(size=12),legend.key.size = unit(12,"point")))

p1<-ggplot(Bp.vec.subset) + geom_point(aes(x=PC1,y=PC2,shape=SampleGroup,color=SampleGroup),size=3) + theme_bw() + xlab(paste("PC1: ",round(as.numeric(Bp.val[1,1]),1),"%",sep="")) + ylab(paste("PC2: ",round(as.numeric(Bp.val[2,1]),1),"%",sep=""))+ theme(legend.text=element_text(size=12),legend.key.size = unit(12,"point"),axis.text=element_text(size=16),axis.title=element_text(size=16))

pdf("/Users/jmattick/Documents/BPPCa.pdf")
print(p1)
dev.off()
Bp.vec.noclin<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/PCA/BP.plink2.noclin.eigenvec",header = TRUE,stringsAsFactors = FALSE))
Bp.val.noclin<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/PCA/BP.plink2.noclin.eigenval",header = FALSE,stringsAsFactors = FALSE))

Bp.vec.subset.noclin<-Bp.vec.noclin[,1:4]
Bp.vec.subset.noclin$SampleGroup<-c("MultiAF",rep("AM",3),rep("FR3",4))

p2<-ggplot(Bp.vec.subset.noclin) + geom_point(aes(x=PC1,y=PC2,shape=SampleGroup,color=SampleGroup),size=4) + theme_bw() + xlab(paste("PC1: ",round(as.numeric(Bp.val.noclin[1,1]),1),"%",sep="")) + ylab(paste("PC2: ",round(as.numeric(Bp.val.noclin[2,1]),1),"%",sep=""))


Bm.vec<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/PCA/BM.plink2.eigenvec",header = TRUE,stringsAsFactors = FALSE))
Bm.val<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Brugia_SNP_Tables/PCA/BM.plink2.eigenval",header = FALSE,stringsAsFactors = FALSE))

Bm.vec.subset<-Bm.vec[,1:6]
Bm.vec.subset$SampleGroup<-c(rep("Lucknow",4),rep("TRS",2),rep("Thailand",4),rep("FR3",4),rep("TRS",2),rep("Liverpool",4),rep("WashU",6))

p3<-ggplot(Bm.vec.subset) + geom_point(aes(x=PC1,y=PC2,shape=SampleGroup,color=SampleGroup),size=3) + theme_bw() + xlab(paste("PC1: ",round(as.numeric(Bm.val[1,1]),1),"%",sep="")) + ylab(paste("PC2: ",round(as.numeric(Bm.val[2,1]),1),"%",sep=""))+ theme(legend.text=element_text(size=12),legend.key.size = unit(12,"point"),axis.text=element_text(size=16),axis.title=element_text(size=16))
p4<-ggplot(Bm.vec.subset) + geom_point(aes(x=PC2,y=PC3,shape=SampleGroup,color=SampleGroup),size=1) + theme_bw() + xlab(paste("PC2: ",round(as.numeric(Bm.val[2,1]),1),"%",sep="")) + ylab(paste("PC3: ",round(as.numeric(Bm.val[3,1]),1),"%",sep=""))
p5<-ggplot(Bm.vec.subset) + geom_point(aes(x=PC3,y=PC4,shape=SampleGroup,color=SampleGroup),size=1) + theme_bw() + xlab(paste("PC3: ",round(as.numeric(Bm.val[3,1]),1),"%",sep="")) + ylab(paste("PC4: ",round(as.numeric(Bm.val[4,1]),1),"%",sep=""))


pdf("/Users/jmattick/Documents/BMPCa.pdf")
print(p3)
dev.off()

pdf("/Users/jmattick/Documents/Pop.Genome.Figure3.pdf",useDingbats=FALSE)
print(p1)
dev.off()
pdf("/Users/jmattick/Documents/Pop.Genome.Figure4.pdf",useDingbats=FALSE)
print(p3)
dev.off()

pdf("/Users/jmattick/Documents/Pop.Genome.All_BMPCAs.pdf",useDingbats=FALSE)
print(p3)
print(p4)
print(p5)

dev.off()


###Box Plots
##BP vs BM

BMalayiGraphing$Species<-"BMalayi"
BPahangiGraphing$Species<-"BPahangi"
TotalGraphing<-rbind(BMalayiGraphing,BPahangiGraphing)
pdf("/Users/jmattick/Documents/Pop.Genome.FigureE.pdf",useDingbats=FALSE)

ggplot(TotalGraphing, aes(x=Species, y=Pi,color=Species)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("SNV Density between ",italic("B. malayi")," and ",italic("B. pahangi"),sep=""))) + theme_bw()
dev.off()


TotalGraphing$Region<-"NonX"
TotalGraphing[TotalGraphing$Species == "BMalayi" & TotalGraphing$Chr == "ChrX" & TotalGraphing$Pos > 4 & TotalGraphing$Pos < 20,]$Region<-"CentralX"
TotalGraphing[TotalGraphing$Species == "BPahangi" & TotalGraphing$subChr == "BP_ChrX_b" & TotalGraphing$Pos > 1 & TotalGraphing$Pos < 17,]$Region<-"CentralX"

q1<-ggplot(TotalGraphing[TotalGraphing$Species == "BMalayi",], aes(x=Pi,color=Region)) + geom_line(stat="density") + ggtitle(expression(paste("SNV Density between non-X and X regions in ",italic("B. malayi"),sep=""))) + theme_bw() + xlim(0,0.0025) + ylab("10kb Windows")
q2<-ggplot(TotalGraphing[TotalGraphing$Species == "BPahangi",], aes(x=Pi,color=Region)) + geom_line(stat="density") + ggtitle(expression(paste("SNV Density between non-X and X regions in ",italic("B. pahangi"),sep=""))) + theme_bw() + xlim(0,0.01) + ylab("10kb Windows")
gl<-list(q1,q2)
pdf("/Users/jmattick/Documents/Pop.Genome.Density.pdf",useDingbats=FALSE)

print(grid.arrange(grobs = gl))

dev.off()
print(mean(TotalGraphing[TotalGraphing$Species == "BMalayi" & TotalGraphing$Region == "CentralX",]$Pi)/mean(TotalGraphing[TotalGraphing$Species == "BMalayi" & TotalGraphing$Region == "NonX",]$Pi))
print(mean(TotalGraphing[TotalGraphing$Species == "BPahangi" & TotalGraphing$Region == "CentralX",]$Pi)/mean(TotalGraphing[TotalGraphing$Species == "BPahangi" & TotalGraphing$Region == "NonX",]$Pi))

###X vs Non X Density


BMDensityGraphing<-BMalayiGraphing
BMDensityGraphing$Region<-"NonX"
BMDensityGraphing[BMDensityGraphing$Chr == "ChrX",]$Region<-"X"
q1<-ggplot(BMDensityGraphing, aes(x=Pi,color=Region)) + geom_line(stat="density") + ggtitle(expression(paste("SNV Density between non-X and X regions in ",italic("B. malayi"),sep=""))) + theme_bw() + xlim(0,0.0025) + ylab("10kb Windows")
pdf("/Users/jmattick/Documents/Pop.Genome.BM.Density.pdf",useDingbats=FALSE,height=4,width=8)
print(q1)
dev.off()

BPDensityGraphing<-BPahangiGraphing
BPDensityGraphing$Region<-"NonX"
BPDensityGraphing[BPDensityGraphing$Chr == "ChrX",]$Region<-"X"
q2<-ggplot(BPDensityGraphing, aes(x=Pi,color=Region)) + geom_line(stat="density") + ggtitle(expression(paste("SNV Density between non-X and X regions in ",italic("B. pahangi"),sep=""))) + theme_bw() + xlim(0,0.0025) + ylab("10kb Windows")
pdf("/Users/jmattick/Documents/Pop.Genome.BP.Density.pdf",useDingbats=FALSE,height=4,width=8)
print(q2)
dev.off()

WBDensityGraphing<-WBancroftiGraphing
WBDensityGraphing$Region<-"NonX"
WBDensityGraphing[WBDensityGraphing$Chr == "ChrX",]$Region<-"X"
q3<-ggplot(WBDensityGraphing, aes(x=Pi,color=Region)) + geom_line(stat="density") + ggtitle(expression(paste("SNV Density between non-X and X regions in ",italic("W. bancrofti"),sep=""))) + theme_bw() + xlim(0,0.0025) + ylab("10kb Windows")
pdf("/Users/jmattick/Documents/Pop.Genome.WB.Density.pdf",useDingbats=FALSE,height=4,width=8)
print(q3)
dev.off()

BTDensityGraphing<-BTimoriGraphing
BTDensityGraphing$Region<-"NonX"
BTDensityGraphing[BTDensityGraphing$Chr == "ChrX",]$Region<-"X"
q4<-ggplot(BTDensityGraphing, aes(x=Pi,color=Region)) + geom_line(stat="density") + ggtitle(expression(paste("SNV Density between non-X and X regions in ",italic("B. timori"),sep=""))) + theme_bw() + xlim(0,0.01) + ylab("10kb Windows")
pdf("/Users/jmattick/Documents/Pop.Genome.BT.Density.pdf",useDingbats=FALSE,height=4,width=8)
print(q4)
dev.off()



gl<-list(q1,q2,q3,q4)
pdf("/Users/jmattick/Documents/Pop.Genome.BMBPWBBT.Density.pdf",useDingbats=FALSE)

print(grid.arrange(grobs = gl))

dev.off()



######FASTA FILE ANALYSIS: GC Content
library("seqinr")

brugia.unmasked <- read.fasta(file = "/Volumes/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Bm.v4.all.fa")
brugia.masked <- read.fasta(file = "/Volumes/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/brugia_malayi.PRJNA10729.WBPS14.genomic_masked.fa")
brugia.genemasked <- read.fasta(file = "/Volumes/projects-t3/EBMAL/JMattick_miRNA_Wolbachia/Brugia.malayi.WGS.genemasked.fa")


AllNames<-getName(brugia.unmasked)

ChrNameVector<-grep("Chr",AllNames)
Graphing<-data.frame()

for (a in ChrNameVector)
{
  subFASTA<-brugia.unmasked[[a]]
  

  ChrSize<-length(getSequence(subFASTA))
  ChrFullName<-getName(subFASTA)
  ChrName<-strsplit(ChrFullName,"_")[[1]][3]
  
  submaskedFASTA<-brugia.masked[[grep(ChrFullName,getName(brugia.masked))]]
  
  subgenemaskedFASTA<-brugia.genemasked[[grep(ChrFullName,getName(brugia.genemasked))]]
  
  
  PartialGraphing<-data.frame(matrix(ncol=5,nrow=ceiling(ChrSize/10000)),stringsAsFactors = FALSE)
  colnames(PartialGraphing)<-c("Chr","Position","GC","Repeats","Genes")
  for (d in 1:ceiling(ChrSize/10000))
  {
    if(d<ceiling(ChrSize/10000)){
      subGC<-GC(getFrag(subFASTA,((d*10000)-9999),(d*10000)))
      if(is.na(subGC)){
        subGC<-0
        subRepeatPerc<-0
        subGenePerc<-0
        
      }else{
        subRepeatTable<-as.data.frame(table(getFrag(submaskedFASTA,((d*10000)-9999),(d*10000))[1:length(getFrag(submaskedFASTA,((d*10000)-9999),(d*10000)))]))
        if(length(subRepeatTable[subRepeatTable$Var1 == "n",]$Freq) > 0){
          subRepeatPerc<-100*as.numeric(subRepeatTable[subRepeatTable$Var1 == "n",]$Freq)/10000
        }else{
          subRepeatPerc<-0
        }
        subGeneTable<-as.data.frame(table(getFrag(subgenemaskedFASTA,((d*10000)-9999),(d*10000))[1:length(getFrag(subgenemaskedFASTA,((d*10000)-9999),(d*10000)))]))
        if(length(subGeneTable[subGeneTable$Var1 == "n",]$Freq) > 0){
          subGenePerc<-100*as.numeric(subGeneTable[subGeneTable$Var1 == "n",]$Freq)/10000
        }else{
          subGenePerc<-0
        }
      }
      PartialGraphing[d,]$Position<-(d*10000)
      PartialGraphing[d,]$Chr<-ChrName
      PartialGraphing[d,]$GC<-100*subGC
      PartialGraphing[d,]$Repeats<-subRepeatPerc
      PartialGraphing[d,]$Genes<-subGenePerc
      
      
      
    } else {
      subGC<-GC(getFrag(subFASTA,((d*10000)-9999),ChrSize))
      if(is.na(subGC)){
        subGC<-0        
        subRepeatPerc<-0
        subGenePerc<-0
        
      }else{
        subRepeatTable<-as.data.frame(table(getFrag(submaskedFASTA,((d*10000)-9999),ChrSize)[1:length(getFrag(submaskedFASTA,((d*10000)-9999),ChrSize))]))
        if(length(subRepeatTable[subRepeatTable$Var1 == "n",]$Freq) > 0){
          subRepeatPerc<-100*as.numeric(subRepeatTable[subRepeatTable$Var1 == "n",]$Freq)/(ChrSize-((d*10000)-10000))
        }else{
          subRepeatPerc<-0
        }
        subGeneTable<-as.data.frame(table(getFrag(subgenemaskedFASTA,((d*10000)-9999),ChrSize)[1:length(getFrag(subgenemaskedFASTA,((d*10000)-9999),ChrSize))]))
        if(length(subGeneTable[subGeneTable$Var1 == "n",]$Freq) > 0){
          subGenePerc<-100*as.numeric(subGeneTable[subGeneTable$Var1 == "n",]$Freq)/(ChrSize-((d*10000)-10000))
        }else{
          subGenePerc<-0
        }
      }
      PartialGraphing[d,]$Position<-ChrSize
      PartialGraphing[d,]$Chr<-ChrName
      PartialGraphing[d,]$GC<-100*subGC
      PartialGraphing[d,]$Repeats<-subRepeatPerc
      PartialGraphing[d,]$Genes<-subGenePerc
      
    }
  }
  
  Graphing<-rbind(Graphing,PartialGraphing)
}

GC_Graphing<-Graphing
write.table(GC_Graphing,"/Users/jmattick/Documents/GC_Graphing.tsv",row.names = FALSE,col.names = T,sep = "\t",quote = FALSE)

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

pdf("/Users/jmattick/Documents/BMalayi.ChromosomeStats.pdf",useDingbats=FALSE)
####GC
g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=GC),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC Plots: A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,50))
g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=GC),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC Plots: B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,50))
g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=GC),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC Plots: C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,50))
g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=GC),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC Plots: D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,50))
g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=GC),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC Plots: E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,50))

#print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))


####Repeat Density
g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=Repeats),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat Plots: A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))
g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=Repeats),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat Plots: B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))
g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=Repeats),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat Plots: C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))
g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=Repeats),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat Plots: D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))
g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=Repeats),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat Plots: E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))

#print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))


####Repeat Density
g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=Genes),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene Plots: A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))
g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=Genes),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene Plots: B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))
g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=Genes),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene Plots: C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))
g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=Genes),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene Plots: D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))
g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=Genes),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene Plots: E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,105))

#print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))


dev.off()

pdf("/Users/jmattick/Documents/BMalayi.ChromosomeStats.vsSNVs.pdf",useDingbats=FALSE)

Combined.GC.SNV.frame<-GC_Graphing
Combined.GC.SNV.frame$Pi<-BMalayiGraphing$Pi

####GC
g1<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr1",]) + geom_point(aes(x=GC,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC vs SNV Plots: A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent GC")+scale_x_continuous(expand = c(0,0),limits=c(0,50))+scale_y_continuous(expand = c(0,0))
g2<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr2",]) + geom_point(aes(x=GC,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC vs SNV Plots: B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent GC")+scale_x_continuous(expand = c(0,0),limits=c(0,50))+scale_y_continuous(expand = c(0,0))
g3<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr3",]) + geom_point(aes(x=GC,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC vs SNV Plots: C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent GC")+scale_x_continuous(expand = c(0,0),limits=c(0,50))+scale_y_continuous(expand = c(0,0))
g4<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr4",]) + geom_point(aes(x=GC,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC vs SNV Plots: D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent GC")+scale_x_continuous(expand = c(0,0),limits=c(0,50))+scale_y_continuous(expand = c(0,0))
g5<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "ChrX",]) + geom_point(aes(x=GC,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("GC vs SNV Plots: E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent GC")+scale_x_continuous(expand = c(0,0),limits=c(0,50))+scale_y_continuous(expand = c(0,0))

#print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,NA))))


####Repeat Density
g1<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr1",]) + geom_point(aes(x=Repeats,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat vs SNV Plots: A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Repeats")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
g2<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr2",]) + geom_point(aes(x=Repeats,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat vs SNV Plots: B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Repeats")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
g3<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr3",]) + geom_point(aes(x=Repeats,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat vs SNV Plots: C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Repeats")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
g4<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr4",]) + geom_point(aes(x=Repeats,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat vs SNV Plots: D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Repeats")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
g5<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "ChrX",]) + geom_point(aes(x=Repeats,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Repeat vs SNV Plots: E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Repeats")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))

#print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,NA))))


####Repeat Density
g1<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr1",]) + geom_point(aes(x=Genes,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene vs SNV Plots: A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Genes")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
g2<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr2",]) + geom_point(aes(x=Genes,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene vs SNV Plots: B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Genes")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
g3<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr3",]) + geom_point(aes(x=Genes,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene vs SNV Plots: C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Genes")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
g4<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "Chr4",]) + geom_point(aes(x=Genes,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene vs SNV Plots: D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Genes")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
g5<-ggplot(Combined.GC.SNV.frame[Combined.GC.SNV.frame$Chr == "ChrX",]) + geom_point(aes(x=Genes,y=Pi),size=0.1) + theme_bw() + facet_grid(~Chr,scales="free",space="free") + ggtitle("Gene vs SNV Plots: E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Percent Genes")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))

#print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,NA))))


dev.off()


###B Malayi Non Coding RNA in 3' Regions
Bm.res<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BMalayi.res.txt",header = FALSE,stringsAsFactors = FALSE))

brugia.gtf <- as.data.frame(read.delim(file = "/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BMalayi_DirectRNA/b_malayi.PRJNA10729.WS259.canonical_geneset.gtf", sep="\t", header=FALSE,stringsAsFactors = F))
brugia.f.depth <- as.data.frame(read.delim(file = "/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BMalayi_DirectRNA/bmalayi_minion_f.depth", sep="\t", header=FALSE,stringsAsFactors = F))
brugia.r.depth <- as.data.frame(read.delim(file = "/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BMalayi_DirectRNA/bmalayi_minion_r.depth", sep="\t", header=FALSE,stringsAsFactors = F))
Illumina.depth <- as.data.frame(read.delim(file = "/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BMalayi_DirectRNA/AF_Illumina_Combined.depth", sep="\t", header=FALSE,stringsAsFactors = F))
brugia.mpileup <- as.data.frame(read.delim(file = "/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BMalayi_DirectRNA/Longread_mpileup_Combined.mpileup", sep="\t", quote="", header=FALSE,stringsAsFactors = F))

brugia.sub.gtf<-brugia.gtf[brugia.gtf$V3 == "gene" | brugia.gtf$V3 == "exon",]

brugia.genelist<-brugia.sub.gtf[brugia.sub.gtf$V3 == "gene",]$V9
brugia.genelist<-brugia.genelist[grep("protein_coding",brugia.genelist)]
brugia.genelist<-gsub(";.*","",brugia.genelist)
brugia.genelist<-gsub("gene_id ","",brugia.genelist)


brugia.split.depth<-brugia.f.depth
brugia.split.depth$V4<-brugia.r.depth$V3


chrlist<-unique(brugia.gtf[brugia.gtf$V3 == "gene",]$V1)
subdepthframe<-data.frame()

for (i in 1:length(chrlist)){
  geneframe<-do.call(rbind, apply(brugia.gtf[brugia.gtf$V3 == "gene" & brugia.gtf$V1 == chrlist[i],], 1, function(x) {
    data.frame(geneloc=seq(as.numeric(x[4])-3000, as.numeric(x[5])+3000))}))
  geneframe$chr<-chrlist[i]
  geneframe<-geneframe[order(geneframe$geneloc),]
  geneframe<-geneframe[geneframe$geneloc > 0,]
  
  genedepth<-brugia.split.depth[(brugia.split.depth$V2 %in% geneframe$geneloc) & brugia.split.depth$V1 == chrlist[i],]
  subdepthframe<-rbind(subdepthframe,genedepth)
  

}


######3' UTR Discovery

GeneUTRFrame<-data.frame(matrix(ncol=9,nrow=length(brugia.genelist)),stringsAsFactors = FALSE)
colnames(GeneUTRFrame)<-c("Gene","Chr","Start","Stop","Orient","ExonLength","ExonDepth","UTR.nonzero.Length","UTR.nonzero.Depth")
count<-1
for (j in 1:length(chrlist)){
  
  brugia.chr.gtf<-brugia.sub.gtf[brugia.sub.gtf$V1 == chrlist[j],]
  
  chr.genelist<-brugia.chr.gtf[brugia.chr.gtf$V3 == "gene",]$V9
  chr.genelist<-chr.genelist[grep("protein_coding",chr.genelist)]
  if(length(chr.genelist) > 0){
    chr.genelist<-gsub(";.*","",chr.genelist)
    chr.genelist<-gsub("gene_id ","",chr.genelist)
    subdepth<-subdepthframe[subdepthframe$V1 == chrlist[j],]
    chrlength<-as.numeric(Bm.res[Bm.res$V1 == chrlist[j],]$V2)
    
    for (p in 1:length(chr.genelist)){
      brugia.gene.gtf<-brugia.chr.gtf[grep(chr.genelist[p],brugia.chr.gtf$V9),]
      orient<-as.character(unique(brugia.gene.gtf$V7))
      
        
      exonframe<-do.call(rbind, apply(brugia.gene.gtf[brugia.gene.gtf$V3 == "exon",], 1, function(x) {
         data.frame(exon=seq(x[4], x[5]))}))
      exonframe$chr<-chrlist[j]
      exonframe<-exonframe[order(exonframe$exon),]
        
      exondepth<-subdepth[subdepth$V2 %in% exonframe$exon,]
        
      if(orient == "+"){
        exonmean<-mean(exondepth$V3)
        exonlength<-length(exondepth$V3)
        UTR.range<-seq(brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V5,brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V5+2000)
        UTR.range<-UTR.range[UTR.range <= chrlength]
        UTRdepth<-subdepth[subdepth$V2 %in% UTR.range,]
        UTR.nonzero.depth<-UTRdepth[UTRdepth$V3 > 0,]
        if(length(UTR.nonzero.depth$V3) > 0){
          utrmean<-mean(UTR.nonzero.depth$V3)
        }else{
          utrmean<-0
        }
        utr.nonzero.length<-as.numeric(length(UTR.nonzero.depth$V3))
          
      }else{
        exonmean<-mean(exondepth$V4)
        exonlength<-length(exondepth$V4)
        UTR.range<-seq(brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V4-2000,brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V4)
        UTR.range<-UTR.range[UTR.range >= 1]
        UTRdepth<-subdepth[subdepth$V2 %in% UTR.range,]
        UTR.nonzero.depth<-UTRdepth[UTRdepth$V4 > 0,]
        if(length(UTR.nonzero.depth$V4) > 0){
          utrmean<-mean(UTR.nonzero.depth$V4)
        }else{
          utrmean<-0
        }
        utr.nonzero.length<-as.numeric(length(UTR.nonzero.depth$V4))
          
      }
      GeneUTRFrame[count,1]<-chr.genelist[p]
      GeneUTRFrame[count,2]<-chrlist[j]
      GeneUTRFrame[count,3]<-as.numeric(brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V4)
      GeneUTRFrame[count,4]<-as.numeric(brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V5)
      GeneUTRFrame[count,5]<-orient
      GeneUTRFrame[count,6]<-exonlength
      GeneUTRFrame[count,7]<-exonmean
      GeneUTRFrame[count,8]<-utr.nonzero.length
      GeneUTRFrame[count,9]<-utrmean
      count<-count+1
    }
  }
}

write.table(GeneUTRFrame,"/Users/jmattick/Documents/GeneUTRFrames.tsv",row.names = FALSE,col.names = T,sep = "\t",quote = FALSE)

FilteredGeneUTRFrame<-GeneUTRFrame

FilteredGeneUTRFrame$Ratio<-(FilteredGeneUTRFrame$UTR.nonzero.Depth+1)/(FilteredGeneUTRFrame$ExonDepth+1)

FilteredGeneUTRFrame<-FilteredGeneUTRFrame[FilteredGeneUTRFrame$Ratio > 5 & FilteredGeneUTRFrame$UTR.nonzero.Depth > 20,]

filteredchrlist<-unique(FilteredGeneUTRFrame$Chr)


####### 5' UTR Discovery
GeneUTRFrame<-data.frame(matrix(ncol=9,nrow=length(brugia.genelist)),stringsAsFactors = FALSE)
colnames(GeneUTRFrame)<-c("Gene","Chr","Start","Stop","Orient","ExonLength","ExonDepth","UTR.nonzero.Length","UTR.nonzero.Depth")
count<-1
for (j in 1:length(chrlist)){
  
  brugia.chr.gtf<-brugia.sub.gtf[brugia.sub.gtf$V1 == chrlist[j],]
  
  chr.genelist<-brugia.chr.gtf[brugia.chr.gtf$V3 == "gene",]$V9
  chr.genelist<-chr.genelist[grep("protein_coding",chr.genelist)]
  if(length(chr.genelist) > 0){
    chr.genelist<-gsub(";.*","",chr.genelist)
    chr.genelist<-gsub("gene_id ","",chr.genelist)
    subdepth<-subdepthframe[subdepthframe$V1 == chrlist[j],]
    chrlength<-as.numeric(Bm.res[Bm.res$V1 == chrlist[j],]$V2)
    
    for (p in 1:length(chr.genelist)){
      brugia.gene.gtf<-brugia.chr.gtf[grep(chr.genelist[p],brugia.chr.gtf$V9),]
      orient<-as.character(unique(brugia.gene.gtf$V7))
      
      
      exonframe<-do.call(rbind, apply(brugia.gene.gtf[brugia.gene.gtf$V3 == "exon",], 1, function(x) {
        data.frame(exon=seq(x[4], x[5]))}))
      exonframe$chr<-chrlist[j]
      exonframe<-exonframe[order(exonframe$exon),]
      
      exondepth<-subdepth[subdepth$V2 %in% exonframe$exon,]
      
      if(orient == "+"){
        exonmean<-mean(exondepth$V3)
        exonlength<-length(exondepth$V3)
        UTR.range<-seq(brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V4-2000,brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V4)
        UTR.range<-UTR.range[UTR.range <= chrlength]
        UTRdepth<-subdepth[subdepth$V2 %in% UTR.range,]
        UTR.nonzero.depth<-UTRdepth[UTRdepth$V3 > 0,]
        if(length(UTR.nonzero.depth$V3) > 0){
          utrmean<-mean(UTR.nonzero.depth$V3)
        }else{
          utrmean<-0
        }
        utr.nonzero.length<-as.numeric(length(UTR.nonzero.depth$V3))
        
      }else{
        exonmean<-mean(exondepth$V4)
        exonlength<-length(exondepth$V4)
        UTR.range<-seq(brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V5,brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V5+2000)
        UTR.range<-UTR.range[UTR.range >= 1]
        UTRdepth<-subdepth[subdepth$V2 %in% UTR.range,]
        UTR.nonzero.depth<-UTRdepth[UTRdepth$V4 > 0,]
        if(length(UTR.nonzero.depth$V4) > 0){
          utrmean<-mean(UTR.nonzero.depth$V4)
        }else{
          utrmean<-0
        }
        utr.nonzero.length<-as.numeric(length(UTR.nonzero.depth$V4))
        
      }
      GeneUTRFrame[count,1]<-chr.genelist[p]
      GeneUTRFrame[count,2]<-chrlist[j]
      GeneUTRFrame[count,3]<-as.numeric(brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V4)
      GeneUTRFrame[count,4]<-as.numeric(brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]$V5)
      GeneUTRFrame[count,5]<-orient
      GeneUTRFrame[count,6]<-exonlength
      GeneUTRFrame[count,7]<-exonmean
      GeneUTRFrame[count,8]<-utr.nonzero.length
      GeneUTRFrame[count,9]<-utrmean
      count<-count+1
    }
  }
}

write.table(GeneUTRFrame,"/Users/jmattick/Documents/5PrimeGeneUTRFrames.tsv",row.names = FALSE,col.names = T,sep = "\t",quote = FALSE)

Filtered5PrimeGeneUTRFrame<-as.data.frame(read.delim(file = "/Users/jmattick/Documents/5PrimeGeneUTRFrames.tsv", sep="\t", header=T,stringsAsFactors = F))

Filtered5PrimeGeneUTRFrame$Ratio<-(Filtered5PrimeGeneUTRFrame$UTR.nonzero.Depth+1)/(Filtered5PrimeGeneUTRFrame$ExonDepth+1)

Filtered5PrimeGeneUTRFrame<-Filtered5PrimeGeneUTRFrame[Filtered5PrimeGeneUTRFrame$Ratio > 5 & Filtered5PrimeGeneUTRFrame$UTR.nonzero.Depth > 20,]

filteredchrlist<-unique(Filtered5PrimeGeneUTRFrame$Chr)




##5' VERSION
#pdf("/Users/jmattick/Documents/BMalayi.NonCodingUTRs.pdf",useDingbats=FALSE)
pdf("/Users/jmattick/Documents/BMalayi.5Prime.NonCodingUTRs.AllIncluded.pdf",useDingbats=FALSE)

###FIX GENE OVERLAP

GeneOverlap.FilteredGeneUTRFrame<-data.frame()
GeneOverlap.f.bedfile<-data.frame()
GeneOverlap.r.bedfile<-data.frame()



for (j in 1:length(filteredchrlist)){
  brugia.chr.gtf<-brugia.sub.gtf[brugia.sub.gtf$V1 == filteredchrlist[j],]
  chr.subdepth<-brugia.split.depth[brugia.split.depth$V1 == filteredchrlist[j],]
  chr.submpileup<-brugia.mpileup[brugia.mpileup$V1 == filteredchrlist[j],]
  chr.subilluminadepth<-Illumina.depth[Illumina.depth$V1 == filteredchrlist[j],]
  
  chrfiltered.UTR<-Filtered5PrimeGeneUTRFrame[Filtered5PrimeGeneUTRFrame$Chr == filteredchrlist[j],]
  genelist<-chrfiltered.UTR$Gene
  for (p in 1:length(genelist)){
    brugia.gene.gtf<-brugia.chr.gtf[grep(genelist[p],brugia.chr.gtf$V9),]
    brugia.geneonly.gtf<-brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]
    brugia.exononly.gtf<-brugia.gene.gtf[brugia.gene.gtf$V3 == "exon",]
    brugia.exononly.gtf$Gene<-genelist[p]
    
    chrfiltered.gene.UTR<-chrfiltered.UTR[chrfiltered.UTR$Gene == genelist[p],]
    orient<-as.character(unique(chrfiltered.gene.UTR$Orient))
    if(orient == "+"){
      start<-as.numeric(brugia.geneonly.gtf$V4-2050)
      end<-as.numeric(brugia.geneonly.gtf$V5+50)
      
      gene.start<-as.numeric(brugia.geneonly.gtf$V4)-1
      
      


      newergenes<-c(brugia.chr.gtf[brugia.chr.gtf$V4 >= start & brugia.chr.gtf$V4 <= gene.start,]$V4,brugia.chr.gtf[brugia.chr.gtf$V5 >= start & brugia.chr.gtf$V5 <= gene.start,]$V5)
      if(length(newergenes) > 0){
        start<-max(newergenes)+1
      }
      
      exonstart<-brugia.chr.gtf[brugia.chr.gtf$V4 >= start & brugia.chr.gtf$V4 <= end,]$V5
      exonend<-brugia.chr.gtf[brugia.chr.gtf$V5 >= start & brugia.chr.gtf$V5 <= end,]$V4
      
      blankframe<-data.frame(matrix(ncol=3,nrow=1),stringsAsFactors = FALSE)
      blankframe[1,1]<-start
      blankframe[1,2]<-end
      blankframe[1,3]<-genelist[p]
      colnames(blankframe)<-c("V4","V5","Gene")
      
      genedepth<-chr.subdepth[chr.subdepth$V2 >= start & chr.subdepth$V2 <= end,]
      genempileup<-chr.submpileup[chr.submpileup$V2 >= start & chr.submpileup$V2 <= end,]
      geneilluminadepth<-chr.subilluminadepth[chr.subilluminadepth$V2 >= start & chr.subilluminadepth$V2 <= end,]
      
      
      colnames(genedepth)<-c("Chr","Pos","Depth_forward","Depth_reverse")
      colnames(genempileup)<-c("Chr","Pos","Ref","Depth_forward","Base_forward","Qual_forward","Depth_reverse","Base_reverse","Qual_reverse")
      colnames(geneilluminadepth)<-c("Chr","Pos","f_1","f_2","f_3","f_4","r_1","r_2","r_3","r_4")
      genempileup$ReadStarts<-0
      genempileup[grep("\\^",genempileup$Base_forward),]$ReadStarts<-nchar(gsub("[^\\^]+", "",genempileup[grep("\\^",genempileup$Base_forward),]$Base_forward))
      
      genempileup$ReadEnds<-0
      genempileup[grep("\\$",genempileup$Base_forward),]$ReadEnds<-nchar(gsub("[^\\$]+", "",genempileup[grep("\\$",genempileup$Base_forward),]$Base_forward))
      
      geneilluminadepth$ForwardDepth<-rowMeans(geneilluminadepth[,3:6])
      if(length(genedepth[genedepth$Pos >= start & genedepth$Pos <= gene.start & genedepth$Depth_forward > 0,]$Depth_forward) > 0){
        UTRDepth<-mean(genedepth[genedepth$Pos >= start & genedepth$Pos <= gene.start & genedepth$Depth_forward > 0,]$Depth_forward)
        
        if(UTRDepth > 20){
          p1<-ggplot(genedepth, aes(x=Pos, y=Depth_forward)) + geom_bar(stat="identity", fill="black") + ggtitle(genelist[p]) + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          p2<-ggplot(genedepth, aes(x=Pos, y=Depth_reverse)) + geom_bar(stat="identity", fill="black") + ggtitle("Opposite Strand Depth") + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          p3<-ggplot(genempileup) + geom_bar(aes(x=Pos, y=ReadStarts),stat="identity", fill="blue",width=5) + geom_bar(aes(x=Pos, y=ReadEnds),stat="identity", fill="red",width=5) + ggtitle("Read Starts and Stops") + theme_bw()
          p4<-ggplot(geneilluminadepth, aes(x=Pos, y=ForwardDepth)) + geom_bar(stat="identity", fill="black") + ggtitle("Illumina Adult Female Depth") + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          
          p5<-ggplot(brugia.exononly.gtf, aes(xmin = V4, xmax = V5, y = Gene)) + geom_gene_arrow() + theme_bw() + theme(axis.text.y = element_text(angle = 90,hjust=0)) + geom_blank(data=blankframe)
          gl<-list(p1,p2,p3,p4,p5)
          print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,1),c(1,1),c(1,1),c(2,2),c(2,2),c(2,2),c(3,3),c(3,3),c(3,3),c(4,4),c(4,4),c(4,4),c(5,5))))
          
          GeneOverlap.FilteredGeneUTRFrame<-rbind(GeneOverlap.FilteredGeneUTRFrame,chrfiltered.gene.UTR)
          
          GeneOverlap.f.sub.bedfile<-data.frame(matrix(ncol=3,nrow=1),stringsAsFactors = FALSE)
          
          GeneOverlap.f.sub.bedfile[1,1]<-filteredchrlist[j]
          GeneOverlap.f.sub.bedfile[1,2]<-gene.start
          GeneOverlap.f.sub.bedfile[1,3]<-end
          
          colnames(GeneOverlap.f.sub.bedfile)<-c("Chr","Start","End")
          
          GeneOverlap.f.bedfile<-rbind(GeneOverlap.f.bedfile,GeneOverlap.f.sub.bedfile)
          
        }
      }
    }else{
      start<-as.numeric(brugia.geneonly.gtf$V4)-50
      end<-as.numeric(brugia.geneonly.gtf$V5)+2000
      gene.end<-as.numeric(brugia.geneonly.gtf$V5)+1
      
      
      
      
      
      
      newergenes<-c(brugia.chr.gtf[brugia.chr.gtf$V4 >= gene.end & brugia.chr.gtf$V4 <= end,]$V4,brugia.chr.gtf[brugia.chr.gtf$V5 >= gene.end & brugia.chr.gtf$V5 <= end,]$V5)
      if(length(newergenes) > 0){
        end<-min(newergenes)-1
      }
      exonstart<-brugia.chr.gtf[brugia.chr.gtf$V4 >= start & brugia.chr.gtf$V4 <= end,]$V4
      exonend<-brugia.chr.gtf[brugia.chr.gtf$V5 >= start & brugia.chr.gtf$V5 <= end,]$V5
      
      blankframe<-data.frame(matrix(ncol=3,nrow=1),stringsAsFactors = FALSE)
      blankframe[1,1]<-start
      blankframe[1,2]<-end
      blankframe[1,3]<-genelist[p]
      
      colnames(blankframe)<-c("V4","V5","Gene")
      
      genedepth<-chr.subdepth[chr.subdepth$V2 >= start & chr.subdepth$V2 <= end,]
      genempileup<-chr.submpileup[chr.submpileup$V2 >= start & chr.submpileup$V2 <= end,]
      geneilluminadepth<-chr.subilluminadepth[chr.subilluminadepth$V2 >= start & chr.subilluminadepth$V2 <= end,]
      
      colnames(genedepth)<-c("Chr","Pos","Depth_forward","Depth_reverse")
      colnames(genempileup)<-c("Chr","Pos","Ref","Depth_forward","Base_forward","Qual_forward","Depth_reverse","Base_reverse","Qual_reverse")
      colnames(geneilluminadepth)<-c("Chr","Pos","f_1","f_2","f_3","f_4","r_1","r_2","r_3","r_4")
      genempileup$ReadStarts<-0
      genempileup[grep("\\^",genempileup$Base_reverse),]$ReadStarts<-nchar(gsub("[^\\^]+", "",genempileup[grep("\\^",genempileup$Base_reverse),]$Base_reverse))
      genempileup$ReadEnds<-0
      genempileup[grep("\\$",genempileup$Base_reverse),]$ReadEnds<-nchar(gsub("[^\\$]+", "",genempileup[grep("\\$",genempileup$Base_reverse),]$Base_reverse))
      
      geneilluminadepth$ReverseDepth<-rowMeans(geneilluminadepth[,7:10])
      
      if(length(genedepth[genedepth$Pos >= gene.end & genedepth$Pos <= end & genedepth$Depth_reverse > 0,]$Depth_reverse) > 0){
        UTRDepth<-mean(genedepth[genedepth$Pos >= gene.end & genedepth$Pos <= end & genedepth$Depth_reverse > 0,]$Depth_reverse)
        if(UTRDepth > 20){
          p1<-ggplot(genedepth, aes(x=Pos, y=Depth_reverse)) + geom_bar(stat="identity", fill="black") + ggtitle(genelist[p]) + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          p2<-ggplot(genedepth, aes(x=Pos, y=Depth_forward)) + geom_bar(stat="identity", fill="black") + ggtitle("Opposite Strand Depth") + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          p3<-ggplot(genempileup) + geom_bar(aes(x=Pos, y=ReadStarts),stat="identity", fill="blue",width=5) + geom_bar(aes(x=Pos, y=ReadEnds),stat="identity", fill="red",width=5) + ggtitle("Read Starts and Stops") + theme_bw()
          p4<-ggplot(geneilluminadepth, aes(x=Pos, y=ReverseDepth)) + geom_bar(stat="identity", fill="black") + ggtitle("Illumina Adult Female Depth") + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          
          p5<-ggplot(brugia.exononly.gtf, aes(xmin = V4, xmax = V5, y = Gene,forward=F)) + geom_gene_arrow() + theme_bw() + theme(axis.text.y = element_text(angle = 90,hjust=0)) + geom_blank(data=blankframe)
          gl<-list(p1,p2,p3,p4,p5)
          print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,1),c(1,1),c(1,1),c(2,2),c(2,2),c(2,2),c(3,3),c(3,3),c(3,3),c(4,4),c(4,4),c(4,4),c(5,5))))
          
          
          GeneOverlap.r.sub.bedfile<-data.frame(matrix(ncol=3,nrow=1),stringsAsFactors = FALSE)
          
          GeneOverlap.r.sub.bedfile[1,1]<-filteredchrlist[j]
          GeneOverlap.r.sub.bedfile[1,2]<-start
          GeneOverlap.r.sub.bedfile[1,3]<-gene.end
          
          colnames(GeneOverlap.r.sub.bedfile)<-c("Chr","Start","End")
          
          GeneOverlap.r.bedfile<-rbind(GeneOverlap.r.bedfile,GeneOverlap.r.sub.bedfile)
          
          
          GeneOverlap.FilteredGeneUTRFrame<-rbind(GeneOverlap.FilteredGeneUTRFrame,chrfiltered.gene.UTR)
        }
      }
    }
  }
}

dev.off()

write.table(GeneOverlap.f.bedfile,"/Users/jmattick/Documents/GeneOverlap.f.5prime.bed",row.names = FALSE,col.names = F,sep = "\t",quote = FALSE)
write.table(GeneOverlap.r.bedfile,"/Users/jmattick/Documents/GeneOverlap.r.5prime.bed",row.names = FALSE,col.names = F,sep = "\t",quote = FALSE)




###3' VERSION

pdf("/Users/jmattick/Documents/BMalayi.NonCodingUTRs.AllIncluded.pdf",useDingbats=FALSE)

###FIX GENE OVERLAP

GeneOverlap.FilteredGeneUTRFrame<-data.frame()
GeneOverlap.f.bedfile<-data.frame()
GeneOverlap.r.bedfile<-data.frame()



for (j in 1:length(filteredchrlist)){
  brugia.chr.gtf<-brugia.sub.gtf[brugia.sub.gtf$V1 == filteredchrlist[j],]
  chr.subdepth<-brugia.split.depth[brugia.split.depth$V1 == filteredchrlist[j],]
  chr.submpileup<-brugia.mpileup[brugia.mpileup$V1 == filteredchrlist[j],]
  chr.subilluminadepth<-Illumina.depth[Illumina.depth$V1 == filteredchrlist[j],]
  
  chrfiltered.UTR<-FilteredGeneUTRFrame[FilteredGeneUTRFrame$Chr == filteredchrlist[j],]
  genelist<-chrfiltered.UTR$Gene
  for (p in 1:length(genelist)){
    brugia.gene.gtf<-brugia.chr.gtf[grep(genelist[p],brugia.chr.gtf$V9),]
    brugia.geneonly.gtf<-brugia.gene.gtf[brugia.gene.gtf$V3 == "gene",]
    brugia.exononly.gtf<-brugia.gene.gtf[brugia.gene.gtf$V3 == "exon",]
    brugia.exononly.gtf$Gene<-genelist[p]
    
    chrfiltered.gene.UTR<-chrfiltered.UTR[chrfiltered.UTR$Gene == genelist[p],]
    orient<-as.character(unique(chrfiltered.gene.UTR$Orient))
    if(orient == "+"){
      start<-as.numeric(brugia.geneonly.gtf$V4)-50
      end<-as.numeric(brugia.geneonly.gtf$V5)+2000
      gene.end<-as.numeric(brugia.geneonly.gtf$V5)+1
      
      
      
      newergenes<-c(brugia.chr.gtf[brugia.chr.gtf$V4 >= gene.end & brugia.chr.gtf$V4 <= end,]$V4,brugia.chr.gtf[brugia.chr.gtf$V5 >= gene.end & brugia.chr.gtf$V5 <= end,]$V5)
      if(length(newergenes) > 0){
        end<-min(newergenes)-1
      }
      exonstart<-brugia.chr.gtf[brugia.chr.gtf$V4 >= start & brugia.chr.gtf$V4 <= end,]$V4
      exonend<-brugia.chr.gtf[brugia.chr.gtf$V5 >= start & brugia.chr.gtf$V5 <= end,]$V5
      
      blankframe<-data.frame(matrix(ncol=3,nrow=1),stringsAsFactors = FALSE)
      blankframe[1,1]<-start
      blankframe[1,2]<-end
      blankframe[1,3]<-genelist[p]
      
      colnames(blankframe)<-c("V4","V5","Gene")
      
      genedepth<-chr.subdepth[chr.subdepth$V2 >= start & chr.subdepth$V2 <= end,]
      genempileup<-chr.submpileup[chr.submpileup$V2 >= start & chr.submpileup$V2 <= end,]
      geneilluminadepth<-chr.subilluminadepth[chr.subilluminadepth$V2 >= start & chr.subilluminadepth$V2 <= end,]
      
      
      colnames(genedepth)<-c("Chr","Pos","Depth_forward","Depth_reverse")
      colnames(genempileup)<-c("Chr","Pos","Ref","Depth_forward","Base_forward","Qual_forward","Depth_reverse","Base_reverse","Qual_reverse")
      colnames(geneilluminadepth)<-c("Chr","Pos","f_1","f_2","f_3","f_4","r_1","r_2","r_3","r_4")
      genempileup$ReadStarts<-0
      genempileup[grep("\\^",genempileup$Base_forward),]$ReadStarts<-nchar(gsub("[^\\^]+", "",genempileup[grep("\\^",genempileup$Base_forward),]$Base_forward))
      
      genempileup$ReadEnds<-0
      genempileup[grep("\\$",genempileup$Base_forward),]$ReadEnds<-nchar(gsub("[^\\$]+", "",genempileup[grep("\\$",genempileup$Base_forward),]$Base_forward))
      
      geneilluminadepth$ForwardDepth<-rowMeans(geneilluminadepth[,3:6])
      if(length(genedepth[genedepth$Pos >= gene.end & genedepth$Pos <= end & genedepth$Depth_forward > 0,]$Depth_forward) > 0){
        UTRDepth<-mean(genedepth[genedepth$Pos >= gene.end & genedepth$Pos <= end & genedepth$Depth_forward > 0,]$Depth_forward)
        if(UTRDepth > 20){
          p1<-ggplot(genedepth, aes(x=Pos, y=Depth_forward)) + geom_bar(stat="identity", fill="black") + ggtitle(genelist[p]) + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          p2<-ggplot(genedepth, aes(x=Pos, y=Depth_reverse)) + geom_bar(stat="identity", fill="black") + ggtitle("Opposite Strand Depth") + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          p3<-ggplot(genempileup) + geom_bar(aes(x=Pos, y=ReadStarts),stat="identity", fill="blue",width=5) + geom_bar(aes(x=Pos, y=ReadEnds),stat="identity", fill="red",width=5) + ggtitle("Read Starts and Stops") + theme_bw()
          p4<-ggplot(geneilluminadepth, aes(x=Pos, y=ForwardDepth)) + geom_bar(stat="identity", fill="black") + ggtitle("Illumina Adult Female Depth") + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          
          p5<-ggplot(brugia.exononly.gtf, aes(xmin = V4, xmax = V5, y = Gene)) + geom_gene_arrow() + theme_bw() + theme(axis.text.y = element_text(angle = 90,hjust=0)) + geom_blank(data=blankframe)
          gl<-list(p1,p2,p3,p4,p5)
          print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,1),c(1,1),c(1,1),c(2,2),c(2,2),c(2,2),c(3,3),c(3,3),c(3,3),c(4,4),c(4,4),c(4,4),c(5,5))))
          
          GeneOverlap.FilteredGeneUTRFrame<-rbind(GeneOverlap.FilteredGeneUTRFrame,chrfiltered.gene.UTR)
          
          GeneOverlap.f.sub.bedfile<-data.frame(matrix(ncol=3,nrow=1),stringsAsFactors = FALSE)
          
          GeneOverlap.f.sub.bedfile[1,1]<-filteredchrlist[j]
          GeneOverlap.f.sub.bedfile[1,2]<-gene.end
          GeneOverlap.f.sub.bedfile[1,3]<-end
          
          colnames(GeneOverlap.f.sub.bedfile)<-c("Chr","Start","End")
          
          GeneOverlap.f.bedfile<-rbind(GeneOverlap.f.bedfile,GeneOverlap.f.sub.bedfile)
          
        }
      }
    }else{
      start<-as.numeric(brugia.geneonly.gtf$V4-2050)
      end<-as.numeric(brugia.geneonly.gtf$V5+50)
      
      gene.start<-as.numeric(brugia.geneonly.gtf$V4)-1
      
      
      
      newergenes<-c(brugia.chr.gtf[brugia.chr.gtf$V4 >= start & brugia.chr.gtf$V4 <= gene.start,]$V4,brugia.chr.gtf[brugia.chr.gtf$V5 >= start & brugia.chr.gtf$V5 <= gene.start,]$V5)
      if(length(newergenes) > 0){
        start<-max(newergenes)+1
      }
      
      exonstart<-brugia.chr.gtf[brugia.chr.gtf$V4 >= start & brugia.chr.gtf$V4 <= end,]$V5
      exonend<-brugia.chr.gtf[brugia.chr.gtf$V5 >= start & brugia.chr.gtf$V5 <= end,]$V4
      
      blankframe<-data.frame(matrix(ncol=3,nrow=1),stringsAsFactors = FALSE)
      blankframe[1,1]<-start
      blankframe[1,2]<-end
      blankframe[1,3]<-genelist[p]
      colnames(blankframe)<-c("V4","V5","Gene")
      
      genedepth<-chr.subdepth[chr.subdepth$V2 >= start & chr.subdepth$V2 <= end,]
      genempileup<-chr.submpileup[chr.submpileup$V2 >= start & chr.submpileup$V2 <= end,]
      geneilluminadepth<-chr.subilluminadepth[chr.subilluminadepth$V2 >= start & chr.subilluminadepth$V2 <= end,]
      
      colnames(genedepth)<-c("Chr","Pos","Depth_forward","Depth_reverse")
      colnames(genempileup)<-c("Chr","Pos","Ref","Depth_forward","Base_forward","Qual_forward","Depth_reverse","Base_reverse","Qual_reverse")
      colnames(geneilluminadepth)<-c("Chr","Pos","f_1","f_2","f_3","f_4","r_1","r_2","r_3","r_4")
      genempileup$ReadStarts<-0
      genempileup[grep("\\^",genempileup$Base_reverse),]$ReadStarts<-nchar(gsub("[^\\^]+", "",genempileup[grep("\\^",genempileup$Base_reverse),]$Base_reverse))
      genempileup$ReadEnds<-0
      genempileup[grep("\\$",genempileup$Base_reverse),]$ReadEnds<-nchar(gsub("[^\\$]+", "",genempileup[grep("\\$",genempileup$Base_reverse),]$Base_reverse))
      
      geneilluminadepth$ReverseDepth<-rowMeans(geneilluminadepth[,7:10])
      
      if(length(genedepth[genedepth$Pos >= start & genedepth$Pos <= gene.start & genedepth$Depth_reverse > 0,]$Depth_reverse) > 0){
        UTRDepth<-mean(genedepth[genedepth$Pos >= start & genedepth$Pos <= gene.start & genedepth$Depth_reverse > 0,]$Depth_reverse)
        if(UTRDepth > 20){
          p1<-ggplot(genedepth, aes(x=Pos, y=Depth_reverse)) + geom_bar(stat="identity", fill="black") + ggtitle(genelist[p]) + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          p2<-ggplot(genedepth, aes(x=Pos, y=Depth_forward)) + geom_bar(stat="identity", fill="black") + ggtitle("Opposite Strand Depth") + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          p3<-ggplot(genempileup) + geom_bar(aes(x=Pos, y=ReadStarts),stat="identity", fill="blue",width=5) + geom_bar(aes(x=Pos, y=ReadEnds),stat="identity", fill="red",width=5) + ggtitle("Read Starts and Stops") + theme_bw()
          p4<-ggplot(geneilluminadepth, aes(x=Pos, y=ReverseDepth)) + geom_bar(stat="identity", fill="black") + ggtitle("Illumina Adult Female Depth") + theme_bw() + geom_vline(xintercept = exonstart, linetype="dotted",color = "blue") + geom_vline(xintercept = exonend, linetype="dotted",color = "red")
          
          p5<-ggplot(brugia.exononly.gtf, aes(xmin = V4, xmax = V5, y = Gene,forward=F)) + geom_gene_arrow() + theme_bw() + theme(axis.text.y = element_text(angle = 90,hjust=0)) + geom_blank(data=blankframe)
          gl<-list(p1,p2,p3,p4,p5)
          print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,1),c(1,1),c(1,1),c(2,2),c(2,2),c(2,2),c(3,3),c(3,3),c(3,3),c(4,4),c(4,4),c(4,4),c(5,5))))
          
          
          GeneOverlap.r.sub.bedfile<-data.frame(matrix(ncol=3,nrow=1),stringsAsFactors = FALSE)
          
          GeneOverlap.r.sub.bedfile[1,1]<-filteredchrlist[j]
          GeneOverlap.r.sub.bedfile[1,2]<-start
          GeneOverlap.r.sub.bedfile[1,3]<-gene.start
          
          colnames(GeneOverlap.r.sub.bedfile)<-c("Chr","Start","End")
          
          GeneOverlap.r.bedfile<-rbind(GeneOverlap.r.bedfile,GeneOverlap.r.sub.bedfile)
          
          
          GeneOverlap.FilteredGeneUTRFrame<-rbind(GeneOverlap.FilteredGeneUTRFrame,chrfiltered.gene.UTR)
        }
      }
    }
  }
}

dev.off()

write.table(GeneOverlap.f.bedfile,"/Users/jmattick/Documents/GeneOverlap.f.bed",row.names = FALSE,col.names = F,sep = "\t",quote = FALSE)
write.table(GeneOverlap.r.bedfile,"/Users/jmattick/Documents/GeneOverlap.r.bed",row.names = FALSE,col.names = F,sep = "\t",quote = FALSE)





################################BUSCO MULTISPECIES TREE#######################################
Bp.busco<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/BPahangi.fulltable.tsv",header = FALSE,stringsAsFactors = FALSE))
Bm.busco<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/BMalayi.fulltable.tsv",header = FALSE,stringsAsFactors = FALSE))
Wb.busco<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/WBancrofti.fulltable.tsv",header = FALSE,stringsAsFactors = FALSE))
Bt.busco<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/BTimori.fulltable.tsv",header = FALSE,stringsAsFactors = FALSE))
Ov.busco<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/OVolvulus.fulltable.tsv",header = FALSE,stringsAsFactors = FALSE))

Bp.busco<-Bp.busco[Bp.busco$V2 == "Complete",]
Bp.busco<-Bp.busco[grep("Chr",Bp.busco$V3),]

Bm.busco<-Bm.busco[Bm.busco$V2 == "Complete",]
Bm.busco<-Bm.busco[grep("Chr",Bm.busco$V3),]

Wb.busco<-Wb.busco[Wb.busco$V2 == "Complete",]

Bt.busco<-Bt.busco[Bt.busco$V2 == "Complete",]

Ov.busco<-Ov.busco[Ov.busco$V2 == "Complete",]


Total.busco<-rbind(Bp.busco,Bm.busco,Wb.busco,Bt.busco,Ov.busco)
BuscoFreq<-as.data.frame(table(Total.busco$V1))

BMBP.busco<-rbind(Bp.busco,Bm.busco)
BMBPBuscoFreq<-as.data.frame(table(BMBP.busco$V1))


CommonBuscos<-BuscoFreq[BuscoFreq$Freq == 5,]
BMBPCommonBuscos<-BMBPBuscoFreq[BMBPBuscoFreq$Freq == 2,]


Bp.busco.all<-Bp.busco[Bp.busco$V1 %in% CommonBuscos$Var1,]
Bm.busco.all<-Bm.busco[Bm.busco$V1 %in% CommonBuscos$Var1,]
Wb.busco.all<-Wb.busco[Wb.busco$V1 %in% CommonBuscos$Var1,]
Bt.busco.all<-Bt.busco[Bt.busco$V1 %in% CommonBuscos$Var1,]
Ov.busco.all<-Ov.busco[Ov.busco$V1 %in% CommonBuscos$Var1,]

Bp.busco.BMBP<-Bp.busco[Bp.busco$V1 %in% BMBPCommonBuscos$Var1,]
Bm.busco.BMBP<-Bm.busco[Bm.busco$V1 %in% BMBPCommonBuscos$Var1,]



FinalFrame<-as.data.frame(cbind(Bp.busco.all$V1,Bp.busco.all$V3,Bm.busco.all$V3,Wb.busco.all$V3,Bt.busco.all$V3,Ov.busco.all$V3))
FinalFrame$V7<-gsub("_.","",gsub("BP_","",FinalFrame$V2))
FinalFrame$V8<-gsub("_scaffold_001|_contig_001","",gsub("Bm_v4_","",FinalFrame$V3))
FinalFrame$V9<-"Same"
FinalFrame[FinalFrame$V7 != FinalFrame$V8,]$V9<-"Shifted"

FinalFrame.BMBP<-as.data.frame(cbind(Bp.busco.BMBP$V1,Bp.busco.BMBP$V3,Bm.busco.BMBP$V3))
FinalFrame.BMBP$V7<-gsub("_.","",gsub("BP_","",FinalFrame.BMBP$V2))
FinalFrame.BMBP$V8<-gsub("_scaffold_001|_contig_001","",gsub("Bm_v4_","",FinalFrame.BMBP$V3))
FinalFrame.BMBP$V9<-"Same"
FinalFrame.BMBP[FinalFrame.BMBP$V7 != FinalFrame.BMBP$V8,]$V9<-"Shifted"


write.table(FinalFrame,"/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/ROutput.busco.multispecies.tsv",row.names = FALSE,col.names = F,sep = "\t",quote = FALSE)
write.table(FinalFrame.BMBP,"/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/ROutput.busco.BMBP.tsv",row.names = FALSE,col.names = F,sep = "\t",quote = FALSE)


###Plot Tree
XChrom.phy<-as.character(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/RAxML_bipartitionsBranchLabels.RaxML.XChrom",header = FALSE,stringsAsFactors = FALSE))
vert.tree.X<-read.tree(text=XChrom.phy)
plot(vert.tree.X,no.margin=TRUE,edge.width=2)
plot(unroot(vert.tree.X),type="unrooted",no.margin=TRUE,lab4ut="axial",edge.width=2)


Autosome.phy<-as.character(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/RAxML_bipartitionsBranchLabels.RaxML.autosome",header = FALSE,stringsAsFactors = FALSE))
vert.tree.auto<-read.tree(text=Autosome.phy)
plot(vert.tree.auto,no.margin=TRUE,edge.width=2)
plot(unroot(vert.tree.auto),type="unrooted",no.margin=TRUE,lab4ut="axial",edge.width=2)

pdf("/Users/jmattick/Documents/RaxMLTrees_FilarialNematodes.pdf",width = 20,height = 10,useDingbats=FALSE)
par(mfrow=c(1, 2))
plot(unroot(vert.tree.X),type="unrooted",no.margin=F,lab4ut="axial",edge.width=2,main="Busco Genes on the X Chromosome")
plot(unroot(vert.tree.auto),type="unrooted",no.margin=F,lab4ut="axial",edge.width=2,main="Busco Genes on the autosome Chromosome")
dev.off()



###B Timori Mito
BMmt_vs_BTmt <- as.data.frame(read.delim("/Users/jmattick/Documents/BTvsBMFinal.coords", sep="\t", header=FALSE,stringsAsFactors = F))
colnames(BMmt_vs_BTmt)<-c("BMStart","BMEnd","BTStart","BTEnd","BMLength","BTLength","PercID","BMSize","BTSize","BMContig","BTContig")
pdf("/Users/jmattick/Documents/BMvsBTMitochondria.pdf")

print(ggplot() + geom_segment(data=BMmt_vs_BTmt,mapping=aes(x=BMStart,xend=BMEnd,y=BTStart,yend=BTEnd,color=PercID)) + ggtitle("B. Timori mitochondria vs. B. malayi mitochondria") + theme_bw()+theme(plot.title = element_text(size=14))+labs(x="BM mitochondria",y="BT Mitochondria")+scale_color_gradientn(colours = rainbow(5)))
dev.off()



#####LD Plotting
LDPlot <- as.data.frame(read.delim("/Users/jmatticki/Documents/ld_window_200kb_filtered.hap.ld.nohead.chr", sep="\t", header=F,stringsAsFactors = F))
colnames(LDPlot)<-c("Chr","BP_A","BP_B","Samples","R2","D","Dprime")
chrs<-unique(LDPlot$Chr)
LDPlot$Diff<-LDPlot$BP_B-LDPlot$BP_A
TotalLDMeanFrame<-data.frame()
pdf("/Users/jmattick/Documents/Chromosome_LD.pdf",height=8,width=16)
increments<-c(1000,10000,20000,40000,60000,100000,200000)
for (j in 1:length(chrs)){
  if(chrs[j] == "Bm_v4_ChrX_scaffold_001"){
    PARThreshold<-c(20000000,24943668)
    PARstart<-1
    for (a in 1:length(PARThreshold))
    {
      subLDPlot<-LDPlot[LDPlot$Chr == chrs[j],]
      subLDPlot<-subLDPlot[subLDPlot$BP_A >= PARstart & subLDPlot$BP_A <= PARThreshold[a],]
      subLDPlot$Diff<-subLDPlot$BP_B-subLDPlot$BP_A
      
      #length<-as.numeric(Bm.chr.res[Bm.chr.res$V1 == chrs[j],]$V2)
      LDMeanFrame<-data.frame(matrix(ncol=4,nrow=length(increments)),stringsAsFactors = FALSE)
      colnames(LDMeanFrame)<-c("Range","AvgLD","Chr","SitePairings")
      ChrName<-strsplit(chrs[j],"_")[[1]][3]
      ChrName<-paste(ChrName,"_",format(PARThreshold[a], scientific=F),sep="")
      
      incrementFrame<-data.frame()
      start<-1
      for (d in 1:length(increments))
      {
        
        #end<-d*1000000
        #start<-((d-1)*1000000)+1
        subincrementLD<-subLDPlot[subLDPlot$Diff >= start & subLDPlot$Diff <= increments[d],]
        subincrementLD$Increment<-increments[d]
        incrementFrame<-rbind(incrementFrame,subincrementLD)
        subLDMean<-mean(subLDPlot[subLDPlot$Diff >= start & subLDPlot$Diff <= increments[d],]$Dprime,na.rm = T)
        LDMeanFrame[d,1]<-increments[d]
        LDMeanFrame[d,2]<-subLDMean
        LDMeanFrame[d,3]<-ChrName
        LDMeanFrame[d,4]<-length(subincrementLD$Chr)
        
        
        start<-increments[d]+1
      }
      my.order<-c("1000","10000","20000","40000","60000","1e+05","2e+05")
      incrementFrame<-na.omit(incrementFrame)
      incrementFrame$Increment<-as.character(incrementFrame$Increment)
      incrementFrame$Increment <- ordered(incrementFrame$Increment, levels = my.order)
      #LDMeanFrame$SitePairings<-LDMeanFrame$SitePairings/max(LDMeanFrame$SitePairings)
      print(ggplot(LDMeanFrame) + geom_point(aes(x=Range, y = AvgLD),size=1) + ggtitle(paste("Linkage Disequilibrium for ",ChrName,sep="")) + theme_bw() + ylim(0,0.8))
      g1<-ggplot(incrementFrame) + geom_density(aes(x=Dprime, color = Increment)) + theme_bw() + ggtitle(paste("Linkage Disequilibrium Distribution for ",ChrName,sep=""))
      LDMeanFramenew<-LDMeanFrame
      LDMeanFramenew$Range<-as.character(LDMeanFramenew$Range)
      LDMeanFramenew$Range <- ordered(LDMeanFramenew$Range, levels = my.order)
      g2<-ggplot(LDMeanFramenew) + geom_bar(aes(x=Range, y=SitePairings),stat="identity") + theme_bw() + ggtitle(paste("Linkage Disequilibrium Total Pairings for ",ChrName,sep=""))
      gl<-list(g1,g2)
      print(grid.arrange(grobs = gl, widths = c(2,1)))
      
      TotalLDMeanFrame<-rbind(TotalLDMeanFrame,LDMeanFrame)
      PARstart<-PARThreshold[a]
    }
  }else{
    subLDPlot<-LDPlot[LDPlot$Chr == chrs[j],]
    subLDPlot$Diff<-subLDPlot$BP_B-subLDPlot$BP_A
    
    #length<-as.numeric(Bm.chr.res[Bm.chr.res$V1 == chrs[j],]$V2)
    LDMeanFrame<-data.frame(matrix(ncol=4,nrow=length(increments)),stringsAsFactors = FALSE)
    colnames(LDMeanFrame)<-c("Range","AvgLD","Chr","SitePairings")
    ChrName<-strsplit(chrs[j],"_")[[1]][3]
  
    incrementFrame<-data.frame()
    start<-1
    for (d in 1:length(increments))
    {
      
      #end<-d*1000000
      #start<-((d-1)*1000000)+1
      subincrementLD<-subLDPlot[subLDPlot$Diff >= start & subLDPlot$Diff <= increments[d],]
      subincrementLD$Increment<-increments[d]
      incrementFrame<-rbind(incrementFrame,subincrementLD)
      subLDMean<-mean(subLDPlot[subLDPlot$Diff >= start & subLDPlot$Diff <= increments[d],]$Dprime,na.rm = T)
      LDMeanFrame[d,1]<-increments[d]
      LDMeanFrame[d,2]<-subLDMean
      LDMeanFrame[d,3]<-ChrName
      LDMeanFrame[d,4]<-length(subincrementLD$Chr)
      
      
      start<-increments[d]+1
    }
    my.order<-c("1000","10000","20000","40000","60000","1e+05","2e+05")
    incrementFrame<-na.omit(incrementFrame)
    incrementFrame$Increment<-as.character(incrementFrame$Increment)
    incrementFrame$Increment <- ordered(incrementFrame$Increment, levels = my.order)
    #LDMeanFrame$SitePairings<-LDMeanFrame$SitePairings/max(LDMeanFrame$SitePairings)
    print(ggplot(LDMeanFrame) + geom_point(aes(x=Range, y = AvgLD),size=1) + ggtitle(paste("Linkage Disequilibrium for ",ChrName,sep="")) + theme_bw() + ylim(0,0.8))
    g1<-ggplot(incrementFrame) + geom_density(aes(x=Dprime, color = Increment)) + theme_bw() + ggtitle(paste("Linkage Disequilibrium Distribution for ",ChrName,sep=""))
    LDMeanFramenew<-LDMeanFrame
    LDMeanFramenew$Range<-as.character(LDMeanFramenew$Range)
    LDMeanFramenew$Range <- ordered(LDMeanFramenew$Range, levels = my.order)
    g2<-ggplot(LDMeanFramenew) + geom_bar(aes(x=Range, y=SitePairings),stat="identity") + theme_bw() + ggtitle(paste("Linkage Disequilibrium Total Pairings for ",ChrName,sep=""))
    gl<-list(g1,g2)
    print(grid.arrange(grobs = gl, widths = c(2,1)))
    
    TotalLDMeanFrame<-rbind(TotalLDMeanFrame,LDMeanFrame)
  }
}

dev.off()

pdf("/Users/jmattick/Documents/Chromosome_LD_Combined.pdf")

print(ggplot(TotalLDMeanFrame) + geom_point(aes(x=Range, y = AvgLD,color = Chr),size=1) + ggtitle("Linkage Disequilibrium across B. malayi") + theme_bw())

dev.off()


####LD for BP
#ld_window_BP_200kb_hardfiltered.hap.ld
#LDPlot <- as.data.frame(read.delim("/Users/jmattick/Documents/ld_window_BP_200kb_hardfiltered.hap.ld.nohead.chr", sep="\t", header=F,stringsAsFactors = F))
#colnames(LDPlot)<-c("Chr","BP_A","BP_B","Samples","R2","D","Dprime")
chrs<-c("Chr1","Chr2","Chr3","Chr4","ChrX")
FileChrs<-c("chr1","chr2","chr3","chr4","chrX")

#LDPlot$Diff<-LDPlot$BP_B-LDPlot$BP_A
TotalLDMeanFrame<-data.frame()
pdf("/Users/jmattick/Documents/Chromosome_LD_BP.pdf",height=8,width=16)
increments<-c(1000,10000,20000,40000,60000,100000,200000)
for (j in 1:length(chrs)){
  #subLDPlot<-LDPlot[LDPlot$Chr == chrs[j],]
  subLDPlot <- as.data.frame(read.delim(paste("/Users/jmattick/Documents/ld_window_BP_200kb_hardfiltered.hap.ld.nohead.chr.real.",FileChrs[j],sep=""), sep="\t", header=F,stringsAsFactors = F))
  colnames(subLDPlot)<-c("Chr","BP_A","BP_B","Samples","R2","D","Dprime")
  subLDPlot$Diff<-subLDPlot$BP_B-subLDPlot$BP_A
  
  #length<-as.numeric(Bm.chr.res[Bm.chr.res$V1 == chrs[j],]$V2)
  LDMeanFrame<-data.frame(matrix(ncol=4,nrow=length(increments)),stringsAsFactors = FALSE)
  colnames(LDMeanFrame)<-c("Range","AvgLD","Chr","SitePairings")
  ChrName<-chrs[j]
  
  incrementFrame<-data.frame()
  start<-1
  for (d in 1:length(increments))
  {
    
    #end<-d*1000000
    #start<-((d-1)*1000000)+1
    subincrementLD<-subLDPlot[subLDPlot$Diff >= start & subLDPlot$Diff <= increments[d],]
    subincrementLD$Increment<-increments[d]
    incrementFrame<-rbind(incrementFrame,subincrementLD)
    subLDMean<-mean(subLDPlot[subLDPlot$Diff >= start & subLDPlot$Diff <= increments[d],]$Dprime,na.rm = T)
    LDMeanFrame[d,1]<-increments[d]
    LDMeanFrame[d,2]<-subLDMean
    LDMeanFrame[d,3]<-ChrName
    LDMeanFrame[d,4]<-length(subincrementLD$Chr)
    
    
    start<-increments[d]+1
  }
  my.order<-c("1000","10000","20000","40000","60000","1e+05","2e+05")
  incrementFrame<-na.omit(incrementFrame)
  incrementFrame$Increment<-as.character(incrementFrame$Increment)
  incrementFrame$Increment <- ordered(incrementFrame$Increment, levels = my.order)
  #LDMeanFrame$SitePairings<-LDMeanFrame$SitePairings/max(LDMeanFrame$SitePairings)
  print(ggplot(LDMeanFrame) + geom_point(aes(x=Range, y = AvgLD),size=1) + ggtitle(paste("Linkage Disequilibrium for B. pahangi ",chrs[j],sep="")) + theme_bw() + ylim(0,0.8))
  g1<-ggplot(incrementFrame) + geom_density(aes(x=Dprime, color = Increment)) + theme_bw() + ggtitle(paste("Linkage Disequilibrium Distribution for ",chrs[j],sep=""))
  LDMeanFramenew<-LDMeanFrame
  LDMeanFramenew$Range<-as.character(LDMeanFramenew$Range)
  LDMeanFramenew$Range <- ordered(LDMeanFramenew$Range, levels = my.order)
  g2<-ggplot(LDMeanFramenew) + geom_bar(aes(x=Range, y=SitePairings),stat="identity") + theme_bw() + ggtitle(paste("Linkage Disequilibrium Total Pairings for ",chrs[j],sep=""))
  gl<-list(g1,g2)
  print(grid.arrange(grobs = gl, widths = c(2,1)))
  
  TotalLDMeanFrame<-rbind(TotalLDMeanFrame,LDMeanFrame)
}

dev.off()
pdf("/Users/jmattick/Documents/Chromosome_LD_Combined_BP.pdf")

print(ggplot(TotalLDMeanFrame) + geom_point(aes(x=Range, y = AvgLD,color = Chr),size=1) + ggtitle("Linkage Disequilibrium across B. pahangi") + theme_bw())

dev.off()



###############TAJIMAS D#############################
library(ggpubr)
BMTajD <- as.data.frame(read.delim("/Users/jmattick/Documents/BM.TajD.tsv", sep="\t", header=F,stringsAsFactors = F))


#B malayi
colnames(BMTajD)<-c("DefaultChr","Pos","SNP_Number","TajD")
BMChrs<-unique(BMTajD$DefaultChr)



BMTajD$Chr<-gsub("_.*","",gsub("Bm_v4_","",BMTajD$DefaultChr))
BMTajD$Pos<-BMTajD$Pos/1000000
BMTajD$SweepRegion<-"Remaining Genome"
BMTajD[BMTajD$Chr == "ChrX" & BMTajD$Pos >= 4 & BMTajD$Pos <= 21,]$SweepRegion<-"Central X Chr"

g1<-ggplot(BMTajD[BMTajD$Chr == "Chr1",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + ggtitle(expression(paste("Tajima's D in ",italic("B. malayi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-4,4))
g2<-ggplot(BMTajD[BMTajD$Chr == "Chr2",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + ggtitle(expression(paste("Tajima's D in ",italic("B. malayi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-4,4))
g3<-ggplot(BMTajD[BMTajD$Chr == "Chr3",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + ggtitle(expression(paste("Tajima's D in ",italic("B. malayi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-4,4))
g4<-ggplot(BMTajD[BMTajD$Chr == "Chr4",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + ggtitle(expression(paste("Tajima's D in ",italic("B. malayi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-4,4))
g5<-ggplot(BMTajD[BMTajD$Chr == "ChrX",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + ggtitle(expression(paste("Tajima's D in ",italic("B. malayi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-4,4))


gl<-list(g1,g2,g3,g4,g5)
pdf("/Users/jmattick/Documents/BM.TajD.Chr.pdf",height=8,width=12,useDingbats = F)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
dev.off()

p1<-ggplot(BMTajD) + geom_point(aes(x=SNP_Number,y=TajD),size=1) + theme_bw() + ggtitle("Tajimas D vs SNP number in B. malayi") + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))
p2<-ggplot(BMTajD, aes(x=SweepRegion, y=TajD,color=SweepRegion)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Tajima's D in the X Sweep Region in ",italic("B. malayi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))
gl<-list(p1,p2)

pdf("/Users/jmattick/Documents/BM.TajD.Supp.pdf",height=8,width=16,useDingbats = F)

print(grid.arrange(grobs = gl, widths = c(2,1)))
dev.off()

BPTajD <- as.data.frame(read.delim("/Users/jmattick/Documents/BP.TajD.tsv", sep="\t", header=F,stringsAsFactors = F))

colnames(BPTajD)<-c("DefaultChr","Pos","SNP_Number","TajD")
BPChrs<-unique(BPTajD$DefaultChr)



BPTajD$Chr<-gsub("_.*","",gsub("BP_","",BPTajD$DefaultChr))
BPTajD$Pos<-BPTajD$Pos/1000000
BPTajD$SweepRegion<-"Remaining Genome"
BPTajD[BPTajD$DefaultChr == "BP_ChrX_b",]$SweepRegion<-"Central X Chr"


#B pahangi

g1<-ggplot(BPTajD[BPTajD$Chr == "Chr1",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + facet_grid(~DefaultChr,scales="free",space="free") + ggtitle(expression(paste("Tajima's D in ",italic("B. pahangi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-3,3))
g2<-ggplot(BPTajD[BPTajD$Chr == "Chr2",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + facet_grid(~DefaultChr,scales="free",space="free") + ggtitle(expression(paste("Tajima's D in ",italic("B. pahangi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-3,3))
g3<-ggplot(BPTajD[BPTajD$Chr == "Chr3",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + facet_grid(~DefaultChr,scales="free",space="free") + ggtitle(expression(paste("Tajima's D in ",italic("B. pahangi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-3,3))
g4<-ggplot(BPTajD[BPTajD$Chr == "Chr4",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + facet_grid(~DefaultChr,scales="free",space="free") + ggtitle(expression(paste("Tajima's D in ",italic("B. pahangi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-3,3))
g5<-ggplot(BPTajD[BPTajD$Chr == "ChrX",]) + geom_point(aes(x=Pos,y=TajD),size=1) + theme_bw() + facet_grid(~DefaultChr,scales="free",space="free") + ggtitle(expression(paste("Tajima's D in ",italic("B. pahangi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-3,3))

gl<-list(g1,g2,g3,g4,g5)
pdf("/Users/jmattick/Documents/BP.TajD.Chr.pdf",height=8,width=12,useDingbats = F)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
dev.off()

p1<-ggplot(BPTajD) + geom_point(aes(x=SNP_Number,y=TajD),size=1) + theme_bw() + ggtitle("Tajimas D vs SNP number in B. pahangi") + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))
p2<-ggplot(BPTajD, aes(x=SweepRegion, y=TajD,color=SweepRegion)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Tajima's D in the X Sweep Region in ",italic("B. pahangi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))
gl<-list(p1,p2)

pdf("/Users/jmattick/Documents/BP.TajD.Supp.pdf",height=8,width=16,useDingbats = F)

print(grid.arrange(grobs = gl, widths = c(2,1)))
dev.off()


p1<-ggplot(BMTajD, aes(x=SweepRegion, y=TajD,color=SweepRegion)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Tajima's D in the X Sweep Region in ",italic("B. malayi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) + stat_compare_means(aes(label = paste("p =", ..p.format..,sep="")), label.x = 1.5, label.y = 4.2)
p2<-ggplot(BPTajD, aes(x=SweepRegion, y=TajD,color=SweepRegion)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Tajima's D in the X Sweep Region in ",italic("B. pahangi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) + stat_compare_means(aes(label = paste("p =", ..p.format..,sep="")), label.x = 1.5, label.y = 2.5)
pdf("/Users/jmattick/Documents/BP.TajD.BoxPlots.pdf",height=8,width=12,useDingbats = F)
gl<-list(p1,p2)

print(grid.arrange(grobs = gl, widths = c(1,1)))
dev.off()

################BM BP Match#############################
Bm_vs_newBP <- as.data.frame(read.delim("/Users/jmattick/Documents/BmvsBp_Final.coords", sep="\t", header=FALSE,stringsAsFactors = F))

BMRelevants<-c("Bm_v4_Chr1_scaffold_001","Bm_v4_Chr2_contig_001","Bm_v4_Chr3_scaffold_001","Bm_v4_Chr4_scaffold_001","Bm_v4_ChrX_scaffold_001","Bm_006")
Matches<-data.frame()
pdf("/Users/jmattick/Documents/BPBMRef_FinalMatch.pdf")
OutputFrame<-data.frame()
for (a in 1:length(BMRelevants))
{
  
  ##Get BP Contigs Attached to BM via Promer
  sub_Bm_vs_newBP<-Bm_vs_newBP[Bm_vs_newBP$V9 == BMRelevants[a],]
  chrLength<-as.numeric(unique(sub_Bm_vs_newBP$V7))
  contigs<-unique(sub_Bm_vs_newBP$V10)
  HighLengthContigs<-data.frame(matrix(ncol=2,nrow=length(contigs)),stringsAsFactors = FALSE)
  colnames(HighLengthContigs)<-c("Contig","Length")
  for (b in 1:length(contigs))
  {
    totallen<-sum(sub_Bm_vs_newBP[sub_Bm_vs_newBP$V10 == contigs[b],]$V6)
    HighLengthContigs[b,1]<-contigs[b]
    HighLengthContigs[b,2]<-totallen
  }
  HighLengthContigs$Length<-HighLengthContigs$Length/chrLength
  HighLengthContigsFiltered<-HighLengthContigs[HighLengthContigs$Length > 0.04,]
  sumcovered<-round(sum(HighLengthContigsFiltered$Length)*100,digits=1)
  
  ##Subset Data Frames
  
  
  
  
  signif_sub_Bm_vs_newBP<-sub_Bm_vs_newBP[sub_Bm_vs_newBP$V10 %in% HighLengthContigsFiltered$Contig,]
  colnames(signif_sub_Bm_vs_newBP)<-c("BMStart","BMEnd","BPStart","BPEnd","BMLength","BPLength","BMSize","BPSize","BMContig","BPContig")
  if(length(signif_sub_Bm_vs_newBP$BMStart) > 0){
    print(ggplot() + geom_segment(data=signif_sub_Bm_vs_newBP,mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=BPContig)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Promer match for ",BMRelevants[a]," vs BP Contigs\nCovering ",sumcovered,"%",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[a],y="BPContigs"))
    OutputFrame<-rbind(OutputFrame,signif_sub_Bm_vs_newBP)
    MatchesSpecific<-HighLengthContigsFiltered
    MatchesSpecific$BMChr<-BMRelevants[a]
    Matches<-rbind(Matches,MatchesSpecific)
    #+ facet_grid(V10~.,scales="free",space="free")
  }
}

dev.off()

write.table(OutputFrame,"/Users/jmattick/Documents/BMBP_FilteredMatch.txt",row.names = FALSE,col.names = FALSE,sep = "\t",quote = FALSE)



###############Windowed FST#############################

BP_FST <- as.data.frame(read.delim("/Users/jmattick/Documents/FST_FR3_vs_Endemic.windowed.weir.fst", sep="\t", header=T,stringsAsFactors = F))
BP_FST<-BP_FST[grep("Chr",BP_FST$CHROM),]
BP_FST$Chr<-gsub("_.*","",gsub("BP_","",BP_FST$CHROM))

pdf("/Users/jmattick/Documents/BP.FST.Chr.pdf",height=8,width=12,useDingbats = F)

ggplot(BP_FST, aes(x=Chr, y=N_VARIANTS,color=Chr)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Comparable Variants per Chromosome in ",italic("B. pahangi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) + stat_compare_means(aes(label = paste("p ", ..p.format..,sep="")), label.x = 3, label.y = 3000)

g1<-ggplot(BP_FST[BP_FST$Chr == "Chr1",]) + geom_point(aes(x=BIN_START,y=N_VARIANTS),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Variant Number in ",italic("B. pahangi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,3600))
g2<-ggplot(BP_FST[BP_FST$Chr == "Chr2",]) + geom_point(aes(x=BIN_START,y=N_VARIANTS),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Variant Number in ",italic("B. pahangi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,3600))
g3<-ggplot(BP_FST[BP_FST$Chr == "Chr3",]) + geom_point(aes(x=BIN_START,y=N_VARIANTS),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Variant Number in ",italic("B. pahangi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,3600))
g4<-ggplot(BP_FST[BP_FST$Chr == "Chr4",]) + geom_point(aes(x=BIN_START,y=N_VARIANTS),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Variant Number in ",italic("B. pahangi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,3600))
g5<-ggplot(BP_FST[BP_FST$Chr == "ChrX",]) + geom_point(aes(x=BIN_START,y=N_VARIANTS),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Variant Number in ",italic("B. pahangi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,3600))

gl<-list(g1,g2,g3,g4,g5)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))

ggplot(BP_FST, aes(x=Chr, y=WEIGHTED_FST,color=Chr)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Fst Distribution (weighted) per Chromosome in ",italic("B. pahangi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) + stat_compare_means(aes(label = paste("p ", ..p.format..,sep="")), label.x = 3, label.y = 1)

BP_FST$BIN_START<-BP_FST$BIN_START/1000000
BP_FST$BIN_END<-BP_FST$BIN_END/1000000

g1<-ggplot(BP_FST[BP_FST$Chr == "Chr1",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g2<-ggplot(BP_FST[BP_FST$Chr == "Chr2",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g3<-ggplot(BP_FST[BP_FST$Chr == "Chr3",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g4<-ggplot(BP_FST[BP_FST$Chr == "Chr4",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g5<-ggplot(BP_FST[BP_FST$Chr == "ChrX",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))

gl<-list(g1,g2,g3,g4,g5)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))

ggplot(BP_FST, aes(x=Chr, y=MEAN_FST,color=Chr)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Fst Distribution (mean) per Chromosome in ",italic("B. pahangi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) + stat_compare_means(aes(label = paste("p ", ..p.format..,sep="")), label.x = 3, label.y = 1)


g1<-ggplot(BP_FST[BP_FST$Chr == "Chr1",]) + geom_point(aes(x=BIN_START,y=MEAN_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g2<-ggplot(BP_FST[BP_FST$Chr == "Chr2",]) + geom_point(aes(x=BIN_START,y=MEAN_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g3<-ggplot(BP_FST[BP_FST$Chr == "Chr3",]) + geom_point(aes(x=BIN_START,y=MEAN_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g4<-ggplot(BP_FST[BP_FST$Chr == "Chr4",]) + geom_point(aes(x=BIN_START,y=MEAN_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g5<-ggplot(BP_FST[BP_FST$Chr == "ChrX",]) + geom_point(aes(x=BIN_START,y=MEAN_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))

gl<-list(g1,g2,g3,g4,g5)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))

dev.off()

##BM FST##

pdf("/Users/jmattick/Documents/BM.FST.Chr.pdf",height=8,width=12,useDingbats = F)


##All pops seperated
BM_FST <- as.data.frame(read.delim("/Users/jmattick/Documents/FST_BMalayi_All.windowed.weir.fst", sep="\t", header=T,stringsAsFactors = F))
BM_FST<-BM_FST[grep("Chr",BM_FST$CHROM),]
BM_FST$Chr<-gsub("_.*","",gsub("Bm_v4_","",BM_FST$CHROM))


ggplot(BM_FST, aes(x=Chr, y=WEIGHTED_FST,color=Chr)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Fst Distribution (all seperated) per Chromosome in ",italic("B. malayi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) + stat_compare_means(aes(label = paste("p ", ..p.format..,sep="")), label.x = 3, label.y = 1)

g1<-ggplot(BM_FST[BM_FST$Chr == "Chr1",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (all seperated) in ",italic("B. malayi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g2<-ggplot(BM_FST[BM_FST$Chr == "Chr2",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (all seperated) in ",italic("B. malayi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g3<-ggplot(BM_FST[BM_FST$Chr == "Chr3",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (all seperated) in ",italic("B. malayi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g4<-ggplot(BM_FST[BM_FST$Chr == "Chr4",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (all seperated) in ",italic("B. malayi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g5<-ggplot(BM_FST[BM_FST$Chr == "ChrX",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (all seperated) in ",italic("B. malayi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))

gl<-list(g1,g2,g3,g4,g5)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))

##PCA pops seperated
BM_FST <- as.data.frame(read.delim("/Users/jmattick/Documents/FST_BMalayi_Lineage.windowed.weir.fst", sep="\t", header=T,stringsAsFactors = F))
BM_FST<-BM_FST[grep("Chr",BM_FST$CHROM),]
BM_FST$Chr<-gsub("_.*","",gsub("Bm_v4_","",BM_FST$CHROM))


ggplot(BM_FST, aes(x=Chr, y=WEIGHTED_FST,color=Chr)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Fst Distribution (PCA seperated) per Chromosome in ",italic("B. malayi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) + stat_compare_means(aes(label = paste("p ", ..p.format..,sep="")), label.x = 3, label.y = 1)

g1<-ggplot(BM_FST[BM_FST$Chr == "Chr1",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (PCA seperated) in ",italic("B. malayi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g2<-ggplot(BM_FST[BM_FST$Chr == "Chr2",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (PCA seperated) in ",italic("B. malayi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g3<-ggplot(BM_FST[BM_FST$Chr == "Chr3",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (PCA seperated) in ",italic("B. malayi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g4<-ggplot(BM_FST[BM_FST$Chr == "Chr4",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (PCA seperated) in ",italic("B. malayi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g5<-ggplot(BM_FST[BM_FST$Chr == "ChrX",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (PCA seperated) in ",italic("B. malayi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))

gl<-list(g1,g2,g3,g4,g5)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))



##Paired pops seperated
BM_FST <- as.data.frame(read.delim("/Users/jmattick/Documents/FST_BMalayi_Pair.windowed.weir.fst", sep="\t", header=T,stringsAsFactors = F))
BM_FST<-BM_FST[grep("Chr",BM_FST$CHROM),]
BM_FST$Chr<-gsub("_.*","",gsub("Bm_v4_","",BM_FST$CHROM))


ggplot(BM_FST, aes(x=Chr, y=WEIGHTED_FST,color=Chr)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("Fst Distribution (Lineage seperated) per Chromosome in ",italic("B. malayi"),sep=""))) + theme_bw() + theme(legend.position="none",plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm")) + stat_compare_means(aes(label = paste("p ", ..p.format..,sep="")), label.x = 3, label.y = 1)


BM_FST$BIN_START<-BM_FST$BIN_START/1000000
BM_FST$BIN_END<-BM_FST$BIN_END/1000000

g1<-ggplot(BM_FST[BM_FST$Chr == "Chr1",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g2<-ggplot(BM_FST[BM_FST$Chr == "Chr2",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g3<-ggplot(BM_FST[BM_FST$Chr == "Chr3",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g4<-ggplot(BM_FST[BM_FST$Chr == "Chr4",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g5<-ggplot(BM_FST[BM_FST$Chr == "ChrX",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))

gl<-list(g1,g2,g3,g4,g5)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))



dev.off()

###Comparative FST

BM_FST <- as.data.frame(read.delim("/Users/jmattick/Documents/FST_BMalayi_Pair.windowed.weir.fst", sep="\t", header=T,stringsAsFactors = F))
BM_FST<-BM_FST[grep("Chr",BM_FST$CHROM),]
BM_FST$Chr<-gsub("_.*","",gsub("Bm_v4_","",BM_FST$CHROM))

BP_FST <- as.data.frame(read.delim("/Users/jmattick/Documents/FST_FR3_vs_Endemic.windowed.weir.fst", sep="\t", header=T,stringsAsFactors = F))
BP_FST<-BP_FST[grep("Chr",BP_FST$CHROM),]
BP_FST$Chr<-gsub("_.*","",gsub("BP_","",BP_FST$CHROM))


Bm_vs_newBP <- as.data.frame(read.delim("/Users/jmattick/Documents/BMBP_FilteredMatch.txt", sep="\t", header=FALSE,stringsAsFactors = F))
colnames(Bm_vs_newBP)<-c("BMStart","BMEnd","BPStart","BPEnd","BMLength","BPLength","BMSize","BPSize","BMContig","BPContig")
FSTComparison<-data.frame()

for (j in 1:length(BM_FST$CHROM)){
  subBm_vs_newBP<-Bm_vs_newBP[Bm_vs_newBP$BMContig == BM_FST[j,]$CHROM,]
  ChrName<-gsub("_.*","",gsub("Bm_v4_","",BM_FST[j,]$CHROM))
  
  subBm_vs_newBP<-subBm_vs_newBP[subBm_vs_newBP$BMStart >= BM_FST[j,]$BIN_START & subBm_vs_newBP$BMEnd <= BM_FST[j,]$BIN_END,]
  if(length(subBm_vs_newBP$BMStart) > 0){
    subBm_vs_newBP<-subBm_vs_newBP[order(subBm_vs_newBP$BMLength,decreasing = T),]
    
    subBP_FST<-BP_FST[BP_FST$CHROM == subBm_vs_newBP[1,]$BPContig,]
    subBP_FST$ModStart<-subBP_FST$BIN_START-subBm_vs_newBP[1,]$BPStart
    subBP_FST<-subBP_FST[subBP_FST$ModStart < 0,]
    subBP_FST<-subBP_FST[order(subBP_FST$ModStart,decreasing=T),]
    
    subFSTComparison<-data.frame(BM_FST[j,]$WEIGHTED_FST,subBP_FST[1,]$WEIGHTED_FST,ChrName,BM_FST[j,]$CHROM,subBP_FST[1,]$Chr,subBP_FST[1,]$CHROM,BM_FST[j,]$BIN_START,subBP_FST[1,]$BIN_START)
    colnames(subFSTComparison)<-c("BM_FST","BP_FST","BM_Chr","BM_Chr_Full","BP_Chr","BP_Chr_Full","BMPos","BPPos")
    FSTComparison<-rbind(FSTComparison,subFSTComparison)
  }
}
pdf("/Users/jmattick/Documents/Nucmer.FSTComparison.pdf",height=8,width=12,useDingbats = F)
ggplot(FSTComparison, aes(x=BM_FST, y=BP_FST,col=BM_Chr)) + geom_point() + ggtitle("Fst Comparison between sequences in BM and BP") + theme_bw() + theme(plot.title = element_text(size=12,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))

dev.off()


######IDENTIFY GENES IN SHARED FST REGIONS######################




High_FSTRegions<-FSTComparison[FSTComparison$BM_FST > 0.5 & FSTComparison$BP_FST > 0.5,]
BM_Anno <- as.data.frame(read.delim("/Users/jmattick/Documents/b_malayi.PRJNA10729.WS259.canonical_geneset.nocomment.gtf", sep="\t", header=F,stringsAsFactors = F))
BM_Anno_genes<-BM_Anno[BM_Anno$V3 == "gene",]
GeneList<-data.frame()
for (j in 1:length(High_FSTRegions$BM_FST)){
  subBM_Anno_genes<-BM_Anno_genes[BM_Anno_genes$V1 == High_FSTRegions[j,]$BM_Chr_Full,]
  BM_start<-as.numeric(High_FSTRegions[j,]$BMPos)
  BM_end<-as.numeric(High_FSTRegions[j,]$BMPos+49999)
  windowed_genes<-subBM_Anno_genes[subBM_Anno_genes$V4 >= BM_start & subBM_Anno_genes$V5 <= BM_end,]
  if(length(windowed_genes$V1) > 0){
    sub_GeneList<-data.frame(gsub(";.*","",gsub("gene_id ","",windowed_genes$V9)),gsub(";","",gsub(".*gene_biotype ","",windowed_genes$V9)))
    colnames(sub_GeneList)<-c("GeneName","GeneType")
    GeneList<-rbind(GeneList,sub_GeneList)
  }
}

#         miRNA protein_coding           tRNA 
#         3            196              1 
coding_GeneList<-GeneList[GeneList$GeneType == "protein_coding",]
rownames(coding_GeneList)<-coding_GeneList$GeneName

High_FSTRegions_new<-High_FSTRegions
High_FSTRegions_new$BMPaste<-paste(High_FSTRegions_new$BM_Chr,"_",High_FSTRegions_new$BMPos,sep="")
High_FSTRegions_new$BPPaste<-paste(High_FSTRegions_new$BP_Chr_Full,"_",High_FSTRegions_new$BPPos,sep="")

pdf("/Users/jmattick/Documents/BM.FST.Chr.pdf",height=8,width=12,useDingbats = F)

BM_FST <- as.data.frame(read.delim("/Users/jmattick/Documents/FST_BMalayi_Pair.windowed.weir.fst", sep="\t", header=T,stringsAsFactors = F))
BM_FST<-BM_FST[grep("Chr",BM_FST$CHROM),]
BM_FST$Chr<-gsub("_.*","",gsub("Bm_v4_","",BM_FST$CHROM))
BM_FST$BMPaste<-paste(BM_FST$Chr,"_",BM_FST$BIN_START,sep="")
BM_FST$Shared<-"Unshared"
BM_FST[BM_FST$BMPaste %in% High_FSTRegions_new$BMPaste,]$Shared<-"Shared"

BM_FST$BIN_START<-BM_FST$BIN_START/1000000
BM_FST$BIN_END<-BM_FST$BIN_END/1000000

g1<-ggplot(BM_FST[BM_FST$Chr == "Chr1",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g2<-ggplot(BM_FST[BM_FST$Chr == "Chr2",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g3<-ggplot(BM_FST[BM_FST$Chr == "Chr3",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g4<-ggplot(BM_FST[BM_FST$Chr == "Chr4",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST,col=Shared),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1)) + scale_color_manual(values=c("red","black"))
g5<-ggplot(BM_FST[BM_FST$Chr == "ChrX",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST,col=Shared),size=1) + theme_bw() + ggtitle(expression(paste("Fst (Lineage seperated) in ",italic("B. malayi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1)) + scale_color_manual(values=c("red","black"))

gl<-list(g1,g2,g3,g4,g5)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))

dev.off()

pdf("/Users/jmattick/Documents/BP.FST.Chr.pdf",height=8,width=12,useDingbats = F)

BP_FST <- as.data.frame(read.delim("/Users/jmattick/Documents/FST_FR3_vs_Endemic.windowed.weir.fst", sep="\t", header=T,stringsAsFactors = F))
BP_FST<-BP_FST[grep("Chr",BP_FST$CHROM),]
BP_FST$Chr<-gsub("_.*","",gsub("BP_","",BP_FST$CHROM))

BP_FST$BPPaste<-paste(BP_FST$CHROM,"_",BP_FST$BIN_START,sep="")
BP_FST$Shared<-"Unshared"
BP_FST[BP_FST$BPPaste %in% High_FSTRegions_new$BPPaste,]$Shared<-"Shared"

BP_FST$BIN_START<-BP_FST$BIN_START/1000000
BP_FST$BIN_END<-BP_FST$BIN_END/1000000

g1<-ggplot(BP_FST[BP_FST$Chr == "Chr1",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr1",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g2<-ggplot(BP_FST[BP_FST$Chr == "Chr2",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr2",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g3<-ggplot(BP_FST[BP_FST$Chr == "Chr3",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr3",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1))
g4<-ggplot(BP_FST[BP_FST$Chr == "Chr4",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST,col=Shared),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": Chr4",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1)) + scale_color_manual(values=c("red","black"))
g5<-ggplot(BP_FST[BP_FST$Chr == "ChrX",]) + geom_point(aes(x=BIN_START,y=WEIGHTED_FST,col=Shared),size=1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle(expression(paste("Fst in ",italic("B. pahangi"),": ChrX",sep=""))) + theme(legend.position="none",plot.title = element_text(size=24,face = "bold"),text = element_text(size=20),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(-0.1,1)) + scale_color_manual(values=c("red","black"))

gl<-list(g1,g2,g3,g4,g5)

print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))

dev.off()

#####Functional Enrichment of Gene Content
modulefunctionaltermenrichment <- function(tpm.de.wgcnamodules, gene.info){
  #tpm.de.wgcnamodules.subset <- tpm.de.wgcnamodules[tpm.de.wgcnamodules$module == clusterofinterest,]
  #tpm.de.wgcnamodules.subset <- tpm.de.wgcnamodules.subset[tpm.de.wgcnamodules.subset$invert == inverted,]
  tpm.de.wgcnamodules.subset<-tpm.de.wgcnamodules
  #cluster.comparisons <- unique(paste0(tpm.de.wgcnamodules$module, tpm.de.wgcnamodules$invert))
  cluster.comparisons<-unique("High_Fst_Regions")
  clusterofinterest<-"Placeholder"
  inverted<-F
  gene.subset <- rownames(tpm.de.wgcnamodules.subset)
  
  for(i in 1:ncol(gene.info)){
    gene.info[,i] <- as.character(gene.info[,i])
  }
  gene.info$inteprodescription[which(is.na(gene.info$inteprodescription))] <- "No InterPro entry"
  gene.info$go_biologicalprocess[which(is.na(gene.info$go_biologicalprocess))] <- "No GO terms for biological process"
  gene.info$go_cellularcomponent[which(is.na(gene.info$go_cellularcomponent))] <- "No GO terms for cellular component"
  gene.info$go_molecularfunction[which(is.na(gene.info$go_molecularfunction))] <- "No GO terms for molecular function"
  
  functionalterms.list <- list(ipr = as.data.frame(table(unlist(strsplit(paste(gene.info$inteprodescription, collapse = "|"),  split = "[|]")))),
                               gobio = as.data.frame(table(unlist(strsplit(paste(gene.info$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                               gocell = as.data.frame(table(unlist(strsplit(paste(gene.info$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                               gomol = as.data.frame(table(unlist(strsplit(paste(gene.info$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  term <- c()
  clusteroccurences <- c()
  genomeoccurences <- c()
  pvalue <- c()
  correctedpvalue <- c()
  oddsratio <- c()
  cluster <- c()
  gene.info.subset <- gene.info[gene.info[,1] %in% gene.subset,]
  functionalterms.list.subset <- list(ipr = as.data.frame(table(unlist(strsplit(paste(gene.info.subset$inteprodescription, collapse = "|"),  split = "[|]")))),
                                      gobio = as.data.frame(table(unlist(strsplit(paste(gene.info.subset$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                                      gocell = as.data.frame(table(unlist(strsplit(paste(gene.info.subset$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                                      gomol = as.data.frame(table(unlist(strsplit(paste(gene.info.subset$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  for(j in 1:length(functionalterms.list)){
    for(k in 1:nrow(functionalterms.list[[j]])){
      freq.all <- functionalterms.list[[j]][k,2]
      if(length(which(functionalterms.list.subset[[j]][,1] == as.character(functionalterms.list[[j]][k,1]))) > 0){
        freq.subset <- functionalterms.list.subset[[j]][which(functionalterms.list.subset[[j]][,1] == as.character(functionalterms.list[[j]][k,1])),2]
      }else{
        freq.subset <- 0
      }
      clustergenes <- nrow(gene.info.subset)
      totalgenes <- nrow(gene.info)
      fisherexact.matrix <- matrix(c(freq.subset, freq.all - freq.subset,
                                     clustergenes - freq.subset, totalgenes - clustergenes - freq.all + freq.subset),
                                   nrow = 2,
                                   ncol = 2)
      fisher.test <- fisher.test(fisherexact.matrix)
      term[length(term) + 1] <- as.character(functionalterms.list[[j]][k,1])
      clusteroccurences[length(clusteroccurences) + 1] <- as.numeric(as.character(freq.subset))
      genomeoccurences[length(genomeoccurences) + 1] <- as.numeric(as.character(freq.all))
      pvalue[length(pvalue) + 1] <- as.numeric(as.character(fisher.test$p.value))
      correctedpvalue[length(correctedpvalue) + 1] <- p.adjust(as.numeric(as.character(fisher.test$p.value)), method = "bonferroni", n = nrow(functionalterms.list[[j]]) * length(cluster.comparisons))
      oddsratio[length(oddsratio) + 1] <- as.numeric(as.character(fisher.test$estimate))
      cluster[length(cluster) + 1] <- as.character(clusterofinterest)
    }
  }
  
  terms.df <- as.data.frame(cbind(term,
                                  clusteroccurences,
                                  genomeoccurences,
                                  pvalue,
                                  correctedpvalue,
                                  oddsratio,
                                  inverted,
                                  cluster))
  terms.df <- terms.df[order(as.numeric(as.character(terms.df$pvalue))),]
  return(terms.df)
}

gene.info <- read.delim("/Users/jmattick/Documents/miRNA_DE_Input_NewHTSeq/bmalayi_gene.info",
                        sep = "\t",
                        header = T)
gene.info <- gene.info[gene.info[,2] == "protein_coding",]
gene.info[,1] <- gsub("Gene:","",gene.info[,1])
gene.info[is.na(gene.info)] <- "Unknown"

terms.df <- as.data.frame(matrix(nrow = 0,ncol = 8))

terms.df <- as.data.frame(rbind(terms.df,modulefunctionaltermenrichment(coding_GeneList, gene.info)))

pvalue.cutoff <- 0.05
#write.table(terms.df,paste0(output_dir,"/clusterfunctionalterms.mRNAtargeted.tsv"),quote = F,row.names = F,col.names = T,sep = "\t")

#terms.df <- terms.df[as.numeric(as.character(terms.df$correctedpvalue)) < pvalue.cutoff,]
terms.df <- terms.df[as.numeric(as.character(terms.df$correctedpvalue)) < pvalue.cutoff,]

write.table(terms.df,"/Users/jmattick/Documents/significant_FST_enrichedterms.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

write.table(coding_GeneList,"/Users/jmattick/Documents/codingGenes_FST.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

CodingFST <- as.data.frame(read.delim("/Users/jmattick/Documents/codingGenes_FST.tsv", sep="\t", header=T,stringsAsFactors = F))

gene.info<-gene.info[gene.info$gene_id %in% CodingFST$GeneName,]
gene.info[is.na(gene.info)] <- "Unknown"
gene.info<-gene.info[order(match(gene.info[,1],CodingFST[,1])),]

write.table(gene.info,"/Users/jmattick/Documents/codingGenes_detailed_FST.tsv",quote = F,row.names = F,col.names = T,sep = "\t")


homologuetable<-as.data.frame(table(gene.info$protein))
homologuetable<-homologuetable[order(homologuetable$Freq,decreasing = T),]
homologuetable$ProteinID <- factor(homologuetable$Var1,levels = homologuetable$Var1)

pdf("/Users/jmattick/Documents/HighFstProteins.pdf",useDingbats=FALSE)

ggplot(data=homologuetable, aes(x=ProteinID, y=Freq)) + geom_bar(stat="identity") + theme_bw() + ggtitle(expression(paste("Proteins with high Fst separation in ",italic("B. malayi")," and ",italic("B. pahangi"),sep=""))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()


#######################INBREEDING#################################
if(!require(detectRUNS)) {
  install.packages("detectRUNS",repos = "http://cran.us.r-project.org")
}
library("detectRUNS")
if(!require(IRanges)) {
  install.packages("IRanges",repos = "http://cran.us.r-project.org")
}
library("IRanges")
if(!require(IRanges)) {
  install.packages("valr",repos = "http://cran.us.r-project.org")
}
library("valr")


genotypeFilePath <- paste("/Users/jmattick/Documents/PLINK_Brugia/",list.files("/Users/jmattick/Documents/PLINK_Brugia/",pattern="ped"),sep="")
mapFilePath <- paste("/Users/jmattick/Documents/PLINK_Brugia/",list.files("/Users/jmattick/Documents/PLINK_Brugia/",pattern="map"),sep="")

slidingRuns <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  windowSize = 15, 
  threshold = 0.05,
  minSNP = 20, 
  ROHet = FALSE, 
  maxOppWindow = 1, 
  maxMissWindow = 1,
  maxGap = 10^6, 
  minLengthBps = 250000, 
  minDensity = 1/10^3, # SNP/kbps
  maxOppRun = NULL,
  maxMissRun = NULL
) 

ChrList<-unique(slidingRuns$chrom)
RunsofHomozygosity<-data.frame()
for (a in 1:length(ChrList))
{
  RangeFrame<-data.frame(ChrList[a],slidingRuns[slidingRuns$chrom == ChrList[a],]$from,slidingRuns[slidingRuns$chrom == ChrList[a],]$to)
  colnames(RangeFrame)<-c("chrom","start","end")
  subRuns<-as.data.frame(bed_merge(RangeFrame, max_dist = 0))
  colnames(subRuns)<-c("Chr","Start","End")
  RunsofHomozygosity<-rbind(RunsofHomozygosity,subRuns)
}

write.table(RunsofHomozygosity,"/Users/jmattick/Documents/RunsofHomozygosity.tsv",row.names = FALSE,col.names = FALSE,sep = "\t",quote = FALSE)



inbreedingFilePath <- paste("/Users/jmattick/Documents/PLINK_Brugia/",list.files("/Users/jmattick/Documents/PLINK_Brugia",pattern="het"),sep="")
InbreedingFrame<-as.data.frame(read.delim(inbreedingFilePath,header = T,stringsAsFactors = FALSE,sep="\t"))
InbreedingFrame$Sample<-paste(InbreedingFrame$FID,InbreedingFrame$IID,sep="_")

HomozygosityTable<-data.frame(paste(slidingRuns$group,slidingRuns$id,sep="_"),slidingRuns$chrom,slidingRuns$from,slidingRuns$to)
colnames(HomozygosityTable)<-c("Sample","Chr","To","From")
Samples<-unique(InbreedingFrame$Sample)

InbreedingvsRoH<-data.frame(matrix(ncol=3,nrow=length(Samples)),stringsAsFactors = FALSE)
colnames(InbreedingvsRoH)<-c("Center","ROH_Count","F_ROH")
for (a in 1:length(Samples))
{
  CenterName<-strsplit(Samples[a],"_")[[1]][1]
  RoHNumber<-as.numeric(length(HomozygosityTable[HomozygosityTable$Sample == Samples[a],]$Sample))
  F_ROH<-as.numeric(InbreedingFrame[InbreedingFrame$Sample == Samples[a],]$F)
  
  if(grepl("Malaysia",Samples[a])){
    Center<-"Malaysia"
  }else if(grepl("TRS",Samples[a])){
    Center<-"TRS"
  }else if(grepl("F3",Samples[a])){
    Center<-"FR3"
  }else if(grepl("UN",Samples[a])){
    Center<-"Liverpool"
  }else if(grepl("W-male",Samples[a])){
    Center<-"WashU"
  }else if(grepl("AM",Samples[a])){
    Center<-"Lucknow"
  }
  InbreedingvsRoH[a,1]<-Center
  InbreedingvsRoH[a,2]<-RoHNumber
  InbreedingvsRoH[a,3]<-F_ROH
  
}
pdf("/Users/jmattick/Documents/BM.allSamples.inbreeding.pdf",height=4,width=8,useDingbats = F)
print(ggplot(InbreedingvsRoH,aes(x=ROH_Count,y=F_ROH,color=Center)) + geom_point(size=3) + theme_bw() + ggtitle("Inbreeding Across Centers")+ theme(plot.title = element_text(size=18),legend.text=element_text(size=12),legend.key.size = unit(12,"point"),axis.text=element_text(size=14),axis.title=element_text(size=16))+xlab("Number of RoH") + ylab("Inbreeding Coefficient"))
dev.off()
g1<-ggplot(InbreedingvsRoH,aes(x=ROH_Count,y=F_ROH,color=Center)) + geom_point(size=3) + theme_bw() + ggtitle("Inbreeding Across Centers")+ theme(plot.title = element_text(size=18),legend.text=element_text(size=12),legend.key.size = unit(12,"point"),axis.text=element_text(size=14),axis.title=element_text(size=16))+xlab("Number of RoH") + ylab("Inbreeding Coefficient")





########B PAHANGI VERSION
genotypeFilePath <- paste("/Users/jmattick/Documents/BP_ROHmaps/",list.files("/Users/jmattick/Documents/BP_ROHmaps/",pattern="ped"),sep="")
mapFilePath <- paste("/Users/jmattick/Documents/BP_ROHmaps/",list.files("/Users/jmattick/Documents/BP_ROHmaps/",pattern="map"),sep="")

slidingRuns <- slidingRUNS.run(
  genotypeFile = genotypeFilePath, 
  mapFile = mapFilePath, 
  windowSize = 15, 
  threshold = 0.05,
  minSNP = 20, 
  ROHet = FALSE, 
  maxOppWindow = 1, 
  maxMissWindow = 1,
  maxGap = 10^6, 
  minLengthBps = 250000, 
  minDensity = 1/10^3, # SNP/kbps
  maxOppRun = NULL,
  maxMissRun = NULL
) 

ChrList<-unique(slidingRuns$chrom)
RunsofHomozygosity<-data.frame()
for (a in 1:length(ChrList))
{
  RangeFrame<-data.frame(ChrList[a],slidingRuns[slidingRuns$chrom == ChrList[a],]$from,slidingRuns[slidingRuns$chrom == ChrList[a],]$to)
  colnames(RangeFrame)<-c("chrom","start","end")
  subRuns<-as.data.frame(bed_merge(RangeFrame, max_dist = 0))
  colnames(subRuns)<-c("Chr","Start","End")
  RunsofHomozygosity<-rbind(RunsofHomozygosity,subRuns)
}

write.table(RunsofHomozygosity,"/Users/jmattick/Documents/BP_RunsofHomozygosity.tsv",row.names = FALSE,col.names = FALSE,sep = "\t",quote = FALSE)




inbreedingFilePath <- paste("/Users/jmattick/Documents/BP_ROHmaps/",list.files("/Users/jmattick/Documents/BP_ROHmaps",pattern="het"),sep="")
InbreedingFrame<-as.data.frame(read.delim(inbreedingFilePath,header = T,stringsAsFactors = FALSE,sep="\t"))

HomozygosityTable<-data.frame(slidingRuns$group,slidingRuns$chrom,slidingRuns$from,slidingRuns$to)
colnames(HomozygosityTable)<-c("Sample","Chr","To","From")
Samples<-unique(InbreedingFrame$FID)
Samples<-Samples[grep("-FR3",Samples,invert = T)]

InbreedingvsRoH<-data.frame(matrix(ncol=3,nrow=length(Samples)),stringsAsFactors = FALSE)
colnames(InbreedingvsRoH)<-c("Center","ROH_Count","F_ROH")
for (a in 1:length(Samples))
{
  RoHNumber<-as.numeric(length(HomozygosityTable[HomozygosityTable$Sample == Samples[a],]$Sample))
  F_ROH<-as.numeric(InbreedingFrame[InbreedingFrame$FID == Samples[a],]$F)
  
  if(grepl("Bp1AM-",Samples[a])){
    Center<-"FR3"
  }else if(grepl("clinical",Samples[a])){
    Center<-"Clinical"
  }
  InbreedingvsRoH[a,1]<-Center
  InbreedingvsRoH[a,2]<-RoHNumber
  InbreedingvsRoH[a,3]<-F_ROH
  
}
pdf("/Users/jmattick/Documents/BP.allSamples.inbreeding.pdf",height=4,width=8)
print(ggplot(InbreedingvsRoH,aes(x=ROH_Count,y=F_ROH,color=Center)) + geom_point(size=3) + theme_bw() + ggtitle("Inbreeding Across Centers")+ theme(plot.title = element_text(size=18),legend.text=element_text(size=12),legend.key.size = unit(12,"point"),axis.text=element_text(size=14),axis.title=element_text(size=16))+xlab("Number of RoH") + ylab("Inbreeding Coefficient"))
dev.off()
g2<-ggplot(InbreedingvsRoH,aes(x=ROH_Count,y=F_ROH,color=Center)) + geom_point(size=3) + theme_bw() + ggtitle("Inbreeding Across Centers")+ theme(plot.title = element_text(size=18),legend.text=element_text(size=12),legend.key.size = unit(12,"point"),axis.text=element_text(size=14),axis.title=element_text(size=16))+xlab("Number of RoH") + ylab("Inbreeding Coefficient")
g3<-g1+ylim(-1.5,1)
g4<-g2+ylim(-1.5,1)
gl<-list(g4,g3)
pdf("/Users/jmattick/Documents/Both.inbreeding.pdf",height=12,width=10)

print(grid.arrange(grobs = gl))
dev.off()






#######################BP Haploid Calling

BPSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BPahangi_New/All.merged.ref.filteredonly.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))
BPSNPs<-BPSNPs[grepl("Chr",BPSNPs$CHROM),]
BPSNPs$Chromosome<-gsub("_.","",gsub("BP_","",BPSNPs$CHROM))

pi_ymax<-max(BPSNPs$PI)
BPSNPs$Position<-BPSNPs$BIN_START/1000000
g1<-ggplot(BPSNPs[BPSNPs$Chromosome == "Chr1",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g2<-ggplot(BPSNPs[BPSNPs$Chromosome == "Chr2",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g3<-ggplot(BPSNPs[BPSNPs$Chromosome == "Chr3",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g4<-ggplot(BPSNPs[BPSNPs$Chromosome == "Chr4",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g5<-ggplot(BPSNPs[BPSNPs$Chromosome == "ChrX",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
pdf("/Users/jmattick/Documents/BPahangiPi.pdf",width = 20,height = 10,useDingbats=FALSE)

  #print(ggarrange(g1,g2,g3,g4,g5,ncol=1,nrow=5,common.legend = TRUE,legend = "bottom"))
gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
BPahangiGraphing_vcftools<-BPSNPs
  

dev.off()







############Depth Plots###############################
pdf("/Users/jmattick/Documents/PahangiDepthHistograms.pdf")


BpFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Pahangi_Depth/",pattern=".depth")
for (b in 1:length(BpFileList))
{

  BPDepth <- read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Pahangi_Depth/",BpFileList[b],sep = ""), sep="\t", header=FALSE,stringsAsFactors=FALSE)
  
  BPDepth<-BPDepth[grep("Chr",BPDepth$V1),]
  #pdf("/Users/jmattick/Documents/PahangiDepthHistograms.pdf")
  
  print(ggplot() + geom_histogram(data=BPDepth[grep("ChrX",BPDepth$V1),], aes(BPDepth[grep("ChrX",BPDepth$V1),]$V3),alpha = .2,binwidth = 1) + geom_histogram(data=BPDepth[grep("ChrX",BPDepth$V1,invert=T),], aes(BPDepth[grep("ChrX",BPDepth$V1,invert=T),]$V3),alpha = .2,binwidth = 1) + theme_bw() + labs(title=paste("Depth Density across ",BpFileList[b],sep="")) + labs(x="Depth", y="Frequency")+xlim(0,1000)+ylim(0,500000))
    
  #dev.off()
  

}
dev.off()



Bp.res<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BPahangi.res.txt",header = FALSE,stringsAsFactors = FALSE))
Bm.res<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Masked_SNP_Tables_All/BMalayi.res.txt",header = FALSE,stringsAsFactors = FALSE))

Bp.chr.res<-Bp.res[grep("Chr",Bp.res$V1),1:2]
Bm.chr.res<-Bm.res[grep("Chr",Bm.res$V1),1:2]

BpFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Pahangi_Depth/",pattern=".depth")
ChrList<-c("Chr1","Chr2","Chr3","Chr4","ChrX")
pdf("/Users/jmattick/Documents/PahangiDepthPlots.pdf")

for (c in 1:length(BpFileList))
{
  samplename<-gsub("\\.sorted\\.dedup\\.depth","",BpFileList[c])
  WormHist <- read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Pahangi_Depth/",BpFileList[c],sep = ""), sep="\t", header=FALSE,stringsAsFactors=FALSE)
  
  
  WormHistRelevant<-WormHist[WormHist$V1 %in% Bp.chr.res$V1,]
  GenomeAvg<-mean(WormHist$V3)
  WormHistRelevant$V3<-(WormHistRelevant$V3/GenomeAvg)
  Graphing<-data.frame()
  
  for (a in 1:length(ChrList))
  {
    WormSubChr<-Bp.chr.res[grep(ChrList[a],Bp.chr.res$V1),]
    WormHistSub<-WormHistRelevant[WormHistRelevant$V1 %in% WormSubChr$V1,]
    subGraphing<-data.frame()
    for (d in 1:length(WormSubChr$V1))
    {
      WormHistSub<-WormHistRelevant[WormHistRelevant$V1 == WormSubChr[d,]$V1,]
      
      PartialGraphing<-data.frame(matrix(ncol=4,nrow=ceiling(length(WormHistSub$V1)/1000)),stringsAsFactors = FALSE)
      colnames(PartialGraphing)<-c("Position","Avg","SubChr","Chr")
      
      for (b in 1:ceiling(length(WormHistSub$V1)/1000))
      {
        if(b<ceiling(length(WormHistSub$V1)/1000)){
          avg<-mean(WormHistSub[((b*1000)-999):(b*1000),]$V3)*2
          PartialGraphing[b,]$Position<-b*1000
          PartialGraphing[b,]$Avg<-avg
          PartialGraphing[b,]$SubChr<-WormSubChr[d,]$V1
          PartialGraphing[b,]$Chr<-ChrList[a]
        } else {
          endpoint<-length(WormHistSub$V1)
          avg<-mean(WormHistSub[((b*1000)-999):endpoint,]$V3)*2
          PartialGraphing[b,]$Position<-endpoint
          PartialGraphing[b,]$Avg<-avg
          PartialGraphing[b,]$SubChr<-WormSubChr[d,]$V1
          PartialGraphing[b,]$Chr<-ChrList[a]
        }
      }
      subGraphing<-rbind(subGraphing,PartialGraphing)
    }
    Graphing<-rbind(Graphing,subGraphing)
    
  }
  
  colnames(Graphing)<-c("Position","N","SubChr","Chr")
  #pdf(paste("/home/jmattick/RScripts/NewPlots/",Species,"_",Sample,"_SequencingDepthPloidy.pdf",sep=""))
  Graphing$Position<-Graphing$Position/1000000
  g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": Chr1")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": Chr2")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": Chr3")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": Chr4")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": ChrX")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  gl<-list(g1,g2,g3,g4,g5)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
  
  write.table(Graphing,paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenomeDepth/",samplename,".",".Table.tsv",sep=""),row.names = FALSE,col.names = FALSE,sep = "\t")
  
}
dev.off()


BmFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Malayi_Depth/",pattern=".depth")
ChrList<-c("Chr1","Chr2","Chr3","Chr4","ChrX")
pdf("/Users/jmattick/Documents/MalayiDepthPlots.pdf")

for (c in 1:length(BmFileList))
{
  samplename<-gsub("\\.gz\\.sorted\\.dedup\\.depth","",BmFileList[c])
  WormHist <- read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Malayi_Depth/",BmFileList[c],sep = ""), sep="\t", header=FALSE,stringsAsFactors=FALSE)
  
  
  WormHistRelevant<-WormHist[WormHist$V1 %in% Bm.chr.res$V1,]
  GenomeAvg<-mean(WormHist$V3)
  WormHistRelevant$V3<-(WormHistRelevant$V3/GenomeAvg)
  Graphing<-data.frame()
  
  for (a in 1:length(ChrList))
  {
    WormSubChr<-Bm.chr.res[grep(ChrList[a],Bm.chr.res$V1),]
    WormHistSub<-WormHistRelevant[WormHistRelevant$V1 %in% WormSubChr$V1,]
    subGraphing<-data.frame()

    PartialGraphing<-data.frame(matrix(ncol=3,nrow=ceiling(length(WormHistSub$V1)/1000)),stringsAsFactors = FALSE)
    colnames(PartialGraphing)<-c("Position","Avg","Chr")
      
    for (b in 1:ceiling(length(WormHistSub$V1)/1000))
    {
      if(b<ceiling(length(WormHistSub$V1)/1000)){
        avg<-mean(WormHistSub[((b*1000)-999):(b*1000),]$V3)*2
        PartialGraphing[b,]$Position<-b*1000
        PartialGraphing[b,]$Avg<-avg
        PartialGraphing[b,]$Chr<-ChrList[a]
      } else {
        endpoint<-length(WormHistSub$V1)
        avg<-mean(WormHistSub[((b*1000)-999):endpoint,]$V3)*2
        PartialGraphing[b,]$Position<-endpoint
        PartialGraphing[b,]$Avg<-avg
        PartialGraphing[b,]$Chr<-ChrList[a]
      }
    }

    Graphing<-rbind(Graphing,PartialGraphing)
    
  }
  
  colnames(Graphing)<-c("Position","N","Chr")
  #pdf(paste("/home/jmattick/RScripts/NewPlots/",Species,"_",Sample,"_SequencingDepthPloidy.pdf",sep=""))
  Graphing$Position<-Graphing$Position/1000000
  g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": Chr1")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": Chr2")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": Chr3")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": Chr4")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": ChrX")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  gl<-list(g1,g2,g3,g4,g5)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
  
  write.table(Graphing,paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenomeDepth/",samplename,".",".Table.tsv",sep=""),row.names = FALSE,col.names = FALSE,sep = "\t")
  
}
dev.off()


######Replot from file lists

##Bp 


BpFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenomeDepth/",pattern="Bp")
ChrList<-c("Chr1","Chr2","Chr3","Chr4","ChrX")
pdf("/Users/jmattick/Documents/PahangiDepthPlots.pdf")

for (c in 1:length(BpFileList))
{
  samplename<-gsub("\\.\\.Table\\.tsv","",BpFileList[c])
  Graphing <- read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenomeDepth/",BpFileList[c],sep = ""), sep="\t", header=FALSE,stringsAsFactors=FALSE)
  colnames(Graphing)<-c("Position","N","SubChr","Chr")
  g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": Chr1",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": Chr2",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": Chr3",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": Chr4",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + facet_grid(~SubChr,scales="free",space="free") + ggtitle(paste(samplename,": ChrX",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,4))
  gl<-list(g1,g2,g3,g4,g5)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
  
}
dev.off()

pdf("/Users/jmattick/Documents/PahangiHistPlots.pdf")

for (c in 1:length(BpFileList))
{
  samplename<-gsub("\\.\\.Table\\.tsv","",BpFileList[c])
  Graphing <- read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenomeDepth/",BpFileList[c],sep = ""), sep="\t", header=FALSE,stringsAsFactors=FALSE)
  colnames(Graphing)<-c("Position","N","SubChr","Chr")
  Graphing$Chrom<-"Autosome"
  Graphing[Graphing$SubChr == "BP_ChrX_c",]$Chrom<-"PAR"
  Graphing[Graphing$Chr == "ChrX" & Graphing$SubChr != "BP_ChrX_c",]$Chrom<-"X-specific"
  Graphing<-Graphing[Graphing$N < 5,]
  
  print(ggplot(Graphing, aes(x=Chrom, y=N)) + geom_boxplot(notch=TRUE) + ggtitle(paste("Sequencing Depth of ",samplename,sep="")) + theme_bw())
  

}
dev.off()

##BM


BmFileList<-list.files("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenomeDepth/")
BmFileList<-BmFileList[grep("Bp",BmFileList,invert = T)]
ChrList<-c("Chr1","Chr2","Chr3","Chr4","ChrX")
pdf("/Users/jmattick/Documents/MalayiDepthPlots.pdf")

for (c in 1:length(BmFileList))
{
  samplename<-gsub("\\.\\.Table\\.tsv","",BmFileList[c])
  samplename<-gsub("\\.sorted\\.dedup\\.depth","",BmFileList[c])
  
  Graphing <- read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenomeDepth/",BmFileList[c],sep = ""), sep="\t", header=FALSE,stringsAsFactors=FALSE)
  colnames(Graphing)<-c("Position","N","Chr")
  g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": Chr1")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,6))
  g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": Chr2")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,6))
  g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": Chr3")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,6))
  g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": Chr4")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,6))
  g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=N),size=0.1) + theme_bw() + ggtitle(paste(samplename,": ChrX")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,6))
  gl<-list(g1,g2,g3,g4,g5)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))
  
}
dev.off()



pdf("/Users/jmattick/Documents/MalayiHistPlots.pdf")

for (c in 1:length(BmFileList))
{
  samplename<-gsub("\\.\\.Table\\.tsv","",BmFileList[c])
  samplename<-gsub("\\.sorted\\.dedup\\.depth","",BmFileList[c])
  Graphing <- read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/PopGenomeDepth/",BmFileList[c],sep = ""), sep="\t", header=FALSE,stringsAsFactors=FALSE)
  colnames(Graphing)<-c("Position","N","Chr")
  Graphing$Chrom<-"Autosome"
  Graphing[Graphing$Chr == "ChrX",]$Chrom<-"X-specific"
  Graphing[Graphing$Chr == "ChrX" & Graphing$Position <= 20,]$Chrom<-"PAR"
  
  Graphing<-Graphing[Graphing$N < 5,]
  
  print(ggplot(Graphing, aes(x=Chrom, y=N)) + geom_boxplot(notch=TRUE) + ggtitle(paste("Sequencing Depth of ",samplename,sep="")) + theme_bw())
  
  
}
dev.off()






##############Nucmer Plots: BM vs All Genomes
#BP
##Nucmer Run
#/usr/local/packages/mummer/show-coords -qlTH /local/scratch/jmattick/OldtoNew_Nucmer/NucmerBMtoNewEditedComparison.delta > /local/scratch/jmattick/OldtoNew_Nucmer/NucmerBMtoNewEditedComparison.coords

Bm_vs_newBP <- as.data.frame(read.delim("/Users/jmattick/Documents/NucmerBMtoFinalRotated.coords", sep="\t", header=FALSE,stringsAsFactors = F))

BMRelevants<-c("Bm_v4_Chr1_scaffold_001","Bm_v4_Chr2_contig_001","Bm_v4_Chr3_scaffold_001","Bm_v4_Chr4_scaffold_001","Bm_v4_ChrX_scaffold_001","Bm_006","gi|23307675|gb|AF538716.1|")
Matches<-data.frame()
pdf("/Users/jmattick/Documents/BP_V9Combined_BMRef_Final_Rotated_Nucmer.pdf")
Graphing<-data.frame()
AllHighLength<-data.frame()
for (a in 1:length(BMRelevants))
{
  
  ##Get BP Contigs Attached to BM via Promer
  sizefactor<-0.045
  sub_Bm_vs_newBP<-Bm_vs_newBP[Bm_vs_newBP$V10 == BMRelevants[a],]
  chrLength<-as.numeric(unique(sub_Bm_vs_newBP$V8))
  contigs<-unique(sub_Bm_vs_newBP$V11)
  HighLengthContigs<-data.frame(matrix(ncol=3,nrow=length(contigs)),stringsAsFactors = FALSE)
  colnames(HighLengthContigs)<-c("Contig","Length","IdentityLength")
  for (b in 1:length(contigs))
  {
    totallen<-sum(sub_Bm_vs_newBP[sub_Bm_vs_newBP$V11 == contigs[b],]$V6)
    totalidentity<-sum(sub_Bm_vs_newBP[sub_Bm_vs_newBP$V11 == contigs[b],]$V6*(sub_Bm_vs_newBP[sub_Bm_vs_newBP$V11 == contigs[b],]$V7/100))
    HighLengthContigs[b,1]<-contigs[b]
    HighLengthContigs[b,2]<-totallen
    HighLengthContigs[b,3]<-totalidentity
    
    
  }
  HighLengthContigs$Length<-HighLengthContigs$Length/chrLength
  HighLengthContigs$IdentityLength<-HighLengthContigs$IdentityLength/chrLength
  if(a == 7){
    sizefactor<-0.8
  }
  
  HighLengthContigsFiltered<-HighLengthContigs[HighLengthContigs$Length > sizefactor,]
  sumcovered<-round(sum(HighLengthContigsFiltered$Length)*100,digits=1)
  HighLengthContigsFiltered$Chr<-BMRelevants[a]
  AllHighLength<-rbind(AllHighLength,HighLengthContigsFiltered)
  ##Subset Data Frames
  
  
  
  
  signif_sub_Bm_vs_newBP<-sub_Bm_vs_newBP[sub_Bm_vs_newBP$V11 %in% HighLengthContigsFiltered$Contig,]
  colnames(signif_sub_Bm_vs_newBP)<-c("BMStart","BMEnd","BPStart","BPEnd","BMLength","BPLength","PercIdentity","BMSize","BPSize","BMContig","BPContig")
  if(length(signif_sub_Bm_vs_newBP$BMStart) > 0){
    #print(ggplot() + geom_segment(data=signif_sub_Bm_vs_newBP,mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[a]," vs BP Contigs\nCovering ",sumcovered,"%",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[a],y="BPContigs")+scale_color_viridis_c(option = "inferno",limits = c(58, 100)))
    print(ggplot() + geom_segment(data=signif_sub_Bm_vs_newBP,mapping=aes(x=BPStart,xend=BPEnd,y=BMStart,yend=BMEnd,color=PercIdentity)) + facet_grid(.~BPContig,scales="free",space="free") + ggtitle(paste("Nucmer match for BP Contigs vs ",BMRelevants[a],"\nCovering ",sumcovered,"%",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(y=BMRelevants[a],x="BPContigs")+scale_color_viridis_c(option = "inferno",limits = c(58, 100)))
    
    print(min(signif_sub_Bm_vs_newBP$PercIdentity))
    MatchesSpecific<-HighLengthContigsFiltered
    MatchesSpecific$BMChr<-BMRelevants[a]
    Matches<-rbind(Matches,MatchesSpecific)
    #+ facet_grid(V10~.,scales="free",space="free")
    #min=58
    Graphing<-rbind(Graphing,signif_sub_Bm_vs_newBP)
  }
}

dev.off()

TotalPercIdentity<-sum(AllHighLength$IdentityLength)/sum(AllHighLength$Length)

GraphingStorage<-Graphing
Graphing$BMStart<-Graphing$BMStart/1000000
Graphing$BMEnd<-Graphing$BMEnd/1000000
Graphing$BPStart<-Graphing$BPStart/1000000
Graphing$BPEnd<-Graphing$BPEnd/1000000
#Graphing<-GraphingStorage
GraphingStorage<-Graphing

var_width = 3
Graphing <- mutate(Graphing, pretty_varname = str_wrap(BPContig, width = var_width))


lowlimit<-80
p1<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[1]," vs BP Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[1],y="BPContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p2<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[2]," vs BP Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[2],y="BPContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p3<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[3]," vs BP Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[3],y="BPContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p4<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[4]," vs BP Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[4],y="BPContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p5<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[5]," vs BP Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[5],y="BPContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p6<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "gi|23307675|gb|AF538716.1|",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[7]," vs BP Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[7],y="BPContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})


Graphing<-GraphingStorage
Graphing$BPContig_f<-factor(Graphing$BPContig, levels=c('BP_Chr1_d','BP_Chr1_c','BP_Chr1_b','BP_Chr1_a'))
p1<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig_f~.,scales="free",space="free") + ggtitle("A") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[1]," (Mbp)",sep=""),y="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
Graphing$BPContig_f<-factor(Graphing$BPContig, levels=c('BP_Chr2_b','BP_Chr2_a'))
p2<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig_f~.,scales="free",space="free") + ggtitle("B") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[2]," (Mbp)",sep=""),y="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
Graphing$BPContig_f<-factor(Graphing$BPContig, levels=c('BP_Chr3_b','BP_Chr3_a'))
p3<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig_f~.,scales="free",space="free") + ggtitle("C") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[3]," (Mbp)",sep=""),y="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",]$BPContig<-gsub("BP_Chr4","BP_Chr4_a",Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",]$BPContig)
p4<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle("D") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[4]," (Mbp)",sep=""),y="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
Graphing$BPContig_f<-factor(Graphing$BPContig, levels=c('BP_ChrX_c','BP_ChrX_b','BP_ChrX_a'))
p5<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig_f~.,scales="free",space="free") + ggtitle("E") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[5]," (Mbp)",sep=""),y="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
p6<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "gi|23307675|gb|AF538716.1|",],mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity,size=5)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle("F") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[7]," (Mbp)",sep=""),y="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

###Flip Figures
Graphing<-GraphingStorage
Graphing$BPContig_f<-factor(Graphing$BPContig, levels=c('BP_Chr1_a','BP_Chr1_b','BP_Chr1_c','BP_Chr1_d'))
p1<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",],mapping=aes(x=BPStart,xend=BPEnd,y=BMStart,yend=BMEnd,color=PercIdentity,size=5)) + facet_grid(.~BPContig_f,scales="free",space="free") + ggtitle("A") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(y=paste(BMRelevants[1]," (Mbp)",sep=""),x="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
Graphing$BPContig_f<-factor(Graphing$BPContig, levels=c('BP_Chr2_a','BP_Chr2_b'))
p2<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",],mapping=aes(x=BPStart,xend=BPEnd,y=BMStart,yend=BMEnd,color=PercIdentity,size=5)) + facet_grid(.~BPContig_f,scales="free",space="free") + ggtitle("B") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(y=paste(BMRelevants[2]," (Mbp)",sep=""),x="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
Graphing$BPContig_f<-factor(Graphing$BPContig, levels=c('BP_Chr3_a','BP_Chr3_b'))
p3<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",],mapping=aes(x=BPStart,xend=BPEnd,y=BMStart,yend=BMEnd,color=PercIdentity,size=5)) + facet_grid(.~BPContig_f,scales="free",space="free") + ggtitle("C") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(y=paste(BMRelevants[3]," (Mbp)",sep=""),x="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",]$BPContig<-gsub("BP_Chr4","BP_Chr4_a",Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",]$BPContig)
p4<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",],mapping=aes(x=BPStart,xend=BPEnd,y=BMStart,yend=BMEnd,color=PercIdentity,size=5)) + facet_grid(.~BPContig,scales="free",space="free") + ggtitle("D") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(y=paste(BMRelevants[4]," (Mbp)",sep=""),x="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
Graphing$BPContig_f<-factor(Graphing$BPContig, levels=c('BP_ChrX_a','BP_ChrX_b','BP_ChrX_c'))
p5<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",],mapping=aes(x=BPStart,xend=BPEnd,y=BMStart,yend=BMEnd,color=PercIdentity,size=5)) + facet_grid(.~BPContig_f,scales="free",space="free") + ggtitle("E") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(y=paste(BMRelevants[5]," (Mbp)",sep=""),x="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
p6<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "gi|23307675|gb|AF538716.1|",],mapping=aes(x=BPStart,xend=BPEnd,y=BMStart,yend=BMEnd,color=PercIdentity,size=5)) + facet_grid(.~BPContig,scales="free",space="free") + ggtitle("F") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(y=paste(BMRelevants[7]," (Mbp)",sep=""),x="BPContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


mylegend<-g_legend(p1)
library(ggplot2)
library(gridExtra)
library(grid)

#png("/Users/jmattick/Documents/BP_NucmerMatch.png", width = 36, height = 16, units = 'in', res = 300)
pdf("/Users/jmattick/Documents/BP_NucmerMatch.V3.pdf", width = 24, height = 36, useDingbats=FALSE)

gl<-list(p1 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p2 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p3 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p4 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p5 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),mylegend)
grid.arrange(
  grobs = gl,
  widths = c(5,5,2),
  layout_matrix = rbind(c(1,2,6),
                        c(3,4,6),
                        c(5,5,6)))
#ggarrange(p1,p2,p3,p4,p5, ncol=2,nrow=3,common.legend = TRUE,legend = "bottom")

dev.off()


#OV
##Nucmer Run
#/usr/local/packages/mummer/show-coords -qlTH /local/scratch/jmattick/OldtoNew_Nucmer/NucmerBMtoNewEditedComparison.delta > /local/scratch/jmattick/OldtoNew_Nucmer/NucmerBMtoNewEditedComparison.coords

Bm_vs_OV <- as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Nucmer_WBOV/BMvsOV.coords", sep="\t", header=FALSE,stringsAsFactors = F))

BMRelevants<-c("Bm_v4_Chr1_scaffold_001","Bm_v4_Chr2_contig_001","Bm_v4_Chr3_scaffold_001","Bm_v4_Chr4_scaffold_001","Bm_v4_ChrX_scaffold_001")
Matches<-data.frame()
pdf("/Users/jmattick/Documents/OV_Nucmer.pdf")
Graphing<-data.frame()
AllHighLength<-data.frame()
for (a in 1:length(BMRelevants))
{
  
  ##Get BP Contigs Attached to BM via Promer
  sizefactor<-0.0045
  sub_Bm_vs_OV<-Bm_vs_OV[Bm_vs_OV$V10 == BMRelevants[a],]
  chrLength<-as.numeric(unique(sub_Bm_vs_OV$V8))
  contigs<-unique(sub_Bm_vs_OV$V11)
  HighLengthContigs<-data.frame(matrix(ncol=3,nrow=length(contigs)),stringsAsFactors = FALSE)
  colnames(HighLengthContigs)<-c("Contig","Length","IdentityLength")
  for (b in 1:length(contigs))
  {
    totallen<-sum(sub_Bm_vs_OV[sub_Bm_vs_OV$V11 == contigs[b],]$V6)
    totalidentity<-sum(sub_Bm_vs_OV[sub_Bm_vs_OV$V11 == contigs[b],]$V6*(sub_Bm_vs_OV[sub_Bm_vs_OV$V11 == contigs[b],]$V7/100))
    HighLengthContigs[b,1]<-contigs[b]
    HighLengthContigs[b,2]<-totallen
    HighLengthContigs[b,3]<-totalidentity
    
    
  }
  HighLengthContigs$Length<-HighLengthContigs$Length/chrLength
  HighLengthContigs$IdentityLength<-HighLengthContigs$IdentityLength/chrLength
  if(a == 7){
    sizefactor<-0.8
  }
  
  HighLengthContigsFiltered<-HighLengthContigs[HighLengthContigs$Length > sizefactor,]
  sumcovered<-round(sum(HighLengthContigsFiltered$Length)*100,digits=1)
  HighLengthContigsFiltered$Chr<-BMRelevants[a]
  AllHighLength<-rbind(AllHighLength,HighLengthContigsFiltered)
  ##Subset Data Frames
  
  
  
  
  signif_sub_Bm_vs_OV<-sub_Bm_vs_OV[sub_Bm_vs_OV$V11 %in% HighLengthContigsFiltered$Contig,]
  colnames(signif_sub_Bm_vs_OV)<-c("BMStart","BMEnd","OVStart","OVEnd","BMLength","OVLength","PercIdentity","BMSize","OVSize","BMContig","OVContig")
  if(length(signif_sub_Bm_vs_OV$BMStart) > 0){
    #print(ggplot() + geom_segment(data=signif_sub_Bm_vs_newBP,mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[a]," vs BP Contigs\nCovering ",sumcovered,"%",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[a],y="BPContigs")+scale_color_viridis_c(option = "inferno",limits = c(58, 100)))
    print(ggplot() + geom_segment(data=signif_sub_Bm_vs_OV,mapping=aes(x=OVStart,xend=OVEnd,y=BMStart,yend=BMEnd,color=PercIdentity)) + facet_grid(.~OVContig,scales="free",space="free") + ggtitle(paste("Nucmer match for OV Contigs vs ",BMRelevants[a],"\nCovering ",sumcovered,"%",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(y=BMRelevants[a],x="OVContigs")+scale_color_viridis_c(option = "inferno",limits = c(58, 100)))
    
    print(min(signif_sub_Bm_vs_OV$PercIdentity))
    MatchesSpecific<-HighLengthContigsFiltered
    MatchesSpecific$BMChr<-BMRelevants[a]
    Matches<-rbind(Matches,MatchesSpecific)
    #+ facet_grid(V10~.,scales="free",space="free")
    #min=58
    Graphing<-rbind(Graphing,signif_sub_Bm_vs_OV)
  }
}

dev.off()

TotalPercIdentity<-sum(AllHighLength$IdentityLength)/sum(AllHighLength$Length)

GraphingStorage<-Graphing
Graphing$BMStart<-Graphing$BMStart/1000000
Graphing$BMEnd<-Graphing$BMEnd/1000000
Graphing$OVStart<-Graphing$OVStart/1000000
Graphing$OVEnd<-Graphing$OVEnd/1000000
#Graphing<-GraphingStorage
GraphingStorage<-Graphing

var_width = 3
Graphing <- mutate(Graphing, pretty_varname = str_wrap(OVContig, width = var_width))


lowlimit<-80
p1<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[1]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[1],y="OVContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p2<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[2]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[2],y="OVContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p3<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[3]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[3],y="OVContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p4<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[4]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[4],y="OVContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p5<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[5]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[5],y="OVContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p6<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "gi|23307675|gb|AF538716.1|",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[7]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[7],y="OVContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})



Graphing<-GraphingStorage

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",]$OVContig)
Graphing$OVContig_f<-factor(Graphing$OVContig, levels=factors)
p1<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig_f~.,scales="free",space="free") + ggtitle("A") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[1]," (Mbp)",sep=""),y="OVContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",]$OVContig)
Graphing$OVContig_f<-factor(Graphing$OVContig, levels=factors)
p2<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig_f~.,scales="free",space="free") + ggtitle("B") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[2]," (Mbp)",sep=""),y="OVContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",]$OVContig)
Graphing$OVContig_f<-factor(Graphing$OVContig, levels=factors)
p3<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig_f~.,scales="free",space="free") + ggtitle("C") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[3]," (Mbp)",sep=""),y="OVContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",]$OVContig)
Graphing$OVContig_f<-factor(Graphing$OVContig, levels=factors)
p4<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig~.,scales="free",space="free") + ggtitle("D") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[4]," (Mbp)",sep=""),y="OVContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",]$OVContig)
Graphing$OVContig_f<-factor(Graphing$OVContig, levels=factors)
p5<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=OVStart,yend=OVEnd,color=PercIdentity,size=5)) + facet_grid(OVContig_f~.,scales="free",space="free") + ggtitle("E") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[5]," (Mbp)",sep=""),y="OVContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


mylegend<-g_legend(p1)
library(ggplot2)
library(gridExtra)
library(grid)

#png("/Users/jmattick/Documents/BP_NucmerMatch.png", width = 36, height = 16, units = 'in', res = 300)
pdf("/Users/jmattick/Documents/OV_NucmerMatch.V3.pdf", width = 24, height = 36, useDingbats=FALSE)

gl<-list(p1 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p2 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p3 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p4 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p5 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),mylegend)
grid.arrange(
  grobs = gl,
  widths = c(5,5,2),
  layout_matrix = rbind(c(1,2,6),
                        c(3,4,6),
                        c(5,5,6)))
#ggarrange(p1,p2,p3,p4,p5, ncol=2,nrow=3,common.legend = TRUE,legend = "bottom")

dev.off()


#WB
##Nucmer Run
#/usr/local/packages/mummer/show-coords -qlTH /local/scratch/jmattick/OldtoNew_Nucmer/NucmerBMtoNewEditedComparison.delta > /local/scratch/jmattick/OldtoNew_Nucmer/NucmerBMtoNewEditedComparison.coords

Bm_vs_WB <- as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Nucmer_WBOV/BMvsWB.coords", sep="\t", header=FALSE,stringsAsFactors = F))
#Bm_vs_WB <- as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/Nucmer_WBOV/BM_vs_WB.pseudo.coords", sep="\t", header=FALSE,stringsAsFactors = F))
#colnames(Bm_vs_WB)<-c("V1","V2","V3","V4","V5","V6","V7","V8","V9","skip","V10","V11")
BMRelevants<-c("Bm_v4_Chr1_scaffold_001","Bm_v4_Chr2_contig_001","Bm_v4_Chr3_scaffold_001","Bm_v4_Chr4_scaffold_001","Bm_v4_ChrX_scaffold_001")
Matches<-data.frame()
pdf("/Users/jmattick/Documents/WB_Nucmer.pdf")
Graphing<-data.frame()
fullcontigmap<-data.frame()
AllHighLength<-data.frame()
for (a in 1:length(BMRelevants))
{
  
  ##Get BP Contigs Attached to BM via Promer
  sizefactor<-0.0015
  sub_Bm_vs_WB<-Bm_vs_WB[Bm_vs_WB$V10 == BMRelevants[a],]
  chrLength<-as.numeric(unique(sub_Bm_vs_WB$V8))
  contigs<-unique(sub_Bm_vs_WB$V11)
  HighLengthContigs<-data.frame(matrix(ncol=3,nrow=length(contigs)),stringsAsFactors = FALSE)
  colnames(HighLengthContigs)<-c("Contig","Length","IdentityLength")
  for (b in 1:length(contigs))
  {
    totallen<-sum(sub_Bm_vs_WB[sub_Bm_vs_WB$V11 == contigs[b],]$V6)
    totalidentity<-sum(sub_Bm_vs_WB[sub_Bm_vs_WB$V11 == contigs[b],]$V6*(sub_Bm_vs_WB[sub_Bm_vs_WB$V11 == contigs[b],]$V7/100))
    HighLengthContigs[b,1]<-contigs[b]
    HighLengthContigs[b,2]<-totallen
    HighLengthContigs[b,3]<-totalidentity
    
    
  }
  HighLengthContigs$Length<-HighLengthContigs$Length/chrLength
  HighLengthContigs$IdentityLength<-HighLengthContigs$IdentityLength/chrLength
  if(a == 7){
    sizefactor<-0.8
  }
  
  HighLengthContigsFiltered<-HighLengthContigs[HighLengthContigs$Length > sizefactor,]
  sumcovered<-round(sum(HighLengthContigsFiltered$Length)*100,digits=1)
  HighLengthContigsFiltered$Chr<-BMRelevants[a]
  AllHighLength<-rbind(AllHighLength,HighLengthContigsFiltered)
  ##Subset Data Frames
  
  
  
  
  signif_sub_Bm_vs_WB<-sub_Bm_vs_WB[sub_Bm_vs_WB$V11 %in% HighLengthContigsFiltered$Contig,]
  colnames(signif_sub_Bm_vs_WB)<-c("BMStart","BMEnd","WBStart","WBEnd","BMLength","WBLength","PercIdentity","BMSize","WBSize","BMContig","WBContig")
  if(length(signif_sub_Bm_vs_WB$BMStart) > 0){
    #print(ggplot() + geom_segment(data=signif_sub_Bm_vs_newBP,mapping=aes(x=BMStart,xend=BMEnd,y=BPStart,yend=BPEnd,color=PercIdentity)) + facet_grid(BPContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[a]," vs BP Contigs\nCovering ",sumcovered,"%",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[a],y="BPContigs")+scale_color_viridis_c(option = "inferno",limits = c(58, 100)))
    ###Order, Orient, Pseudomoleculize
    contigmap<-data.frame()
    for (c in 1:length(HighLengthContigsFiltered$Contig))
    {
      Orientation<-sum(signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == HighLengthContigsFiltered[c,]$Contig,]$WBEnd - signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == HighLengthContigsFiltered[c,]$Contig,]$WBStart)
      if(Orientation > 0){
        OrientationResult<-"Forward"
      }else{
        OrientationResult<-"Reverse"
      }
      StartPosition<-weighted.mean(signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == HighLengthContigsFiltered[c,]$Contig,]$BMStart,(signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == HighLengthContigsFiltered[c,]$Contig,]$WBLength/signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == HighLengthContigsFiltered[c,]$Contig,]$WBSize))
      subcontigmap<-data.frame(BMRelevants[a],HighLengthContigsFiltered[c,]$Contig,StartPosition,OrientationResult)
      colnames(subcontigmap)<-c("BMChr","Contig","Start","Orientation")
      contigmap<-rbind(contigmap,subcontigmap)
    }
    contigmap<-contigmap[order(contigmap$Start),]
    fullcontigmap<-rbind(fullcontigmap,contigmap)
    WBPseudoName<-paste("WB_Pseudo_",gsub("_.*","",gsub("Bm_v4_","",BMRelevants[a])),sep="")
    totallength<-0
    for (d in 1:length(contigmap$Contig))
    {
      if(contigmap[d,]$Orientation == "Reverse"){
        tempstarts<-signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBStart
        tempend<-signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBEnd
        signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBStart<-tempend
        signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBEnd<-tempstarts
      }
      signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBStart<-signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBStart+totallength
      signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBEnd<-signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBEnd+totallength
      
      totallength<-totallength+as.numeric(unique(signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBSize))
      signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == contigmap[d,]$Contig,]$WBContig<-WBPseudoName
    }
    signif_sub_Bm_vs_WB[signif_sub_Bm_vs_WB$WBContig == WBPseudoName,]$WBSize<-totallength
    print(ggplot() + geom_segment(data=signif_sub_Bm_vs_WB,mapping=aes(x=WBStart,xend=WBEnd,y=BMStart,yend=BMEnd,color=PercIdentity)) + ggtitle(paste("Nucmer match for WB Contigs vs ",BMRelevants[a],"\nCovering ",sumcovered,"%",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(y=BMRelevants[a],x=WBPseudoName)+scale_color_viridis_c(option = "inferno",limits = c(58, 100)))
    
    print(min(signif_sub_Bm_vs_WB$PercIdentity))
    MatchesSpecific<-HighLengthContigsFiltered
    MatchesSpecific$BMChr<-BMRelevants[a]
    Matches<-rbind(Matches,MatchesSpecific)
    #+ facet_grid(V10~.,scales="free",space="free")
    #min=58
    Graphing<-rbind(Graphing,signif_sub_Bm_vs_WB)
  }
}

dev.off()
write.table(fullcontigmap,"/Users/jmattick/Documents/WB_Contig_OrderOrient.tsv",row.names = FALSE,col.names = FALSE,sep = "\t",quote = F)

TotalPercIdentity<-sum(AllHighLength$IdentityLength)/sum(AllHighLength$Length)

GraphingStorage<-Graphing
Graphing<-Graphing[order(Graphing$BMStart),]
Graphing$BMStart<-Graphing$BMStart/1000000
Graphing$BMEnd<-Graphing$BMEnd/1000000
Graphing$WBStart<-Graphing$WBStart/1000000
Graphing$WBEnd<-Graphing$WBEnd/1000000
#Graphing<-GraphingStorage
GraphingStorage<-Graphing

var_width = 3
Graphing <- mutate(Graphing, pretty_varname = str_wrap(WBContig, width = var_width))


lowlimit<-80
p1<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[1]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[1],y="WBContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p2<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[2]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[2],y="WBContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p3<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[3]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[3],y="WBContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p4<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[4]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[4],y="WBContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p5<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[5]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[5],y="WBContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})
p6<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "gi|23307675|gb|AF538716.1|",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig~.,scales="free",space="free") + ggtitle(paste("Nucmer match for ",BMRelevants[7]," vs OV Contigs",sep="")) + theme_bw()+theme(plot.title = element_text(size=14))+labs(x=BMRelevants[7],y="WBContigs")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>80,scales::rescale(x,to = to,from = c(80, max(x, na.rm = TRUE))),0)})



Graphing<-GraphingStorage

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",]$WBContig)
Graphing$WBContig_f<-factor(Graphing$WBContig, levels=factors)
p1<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + ggtitle("A") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[1]," (Mbp)",sep=""),y=paste(unique(Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",]$WBContig)," (Mbp)",sep=""))+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
#p1<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig_f~.,scales="free",space="free") + ggtitle("A") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[1]," (Mbp)",sep=""),y="WBContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",]$WBContig)
Graphing$WBContig_f<-factor(Graphing$WBContig, levels=factors)
p2<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + ggtitle("B") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[2]," (Mbp)",sep=""),y=paste(unique(Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",]$WBContig)," (Mbp)",sep=""))+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
#p2<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig_f~.,scales="free",space="free") + ggtitle("B") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[2]," (Mbp)",sep=""),y="WBContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",]$WBContig)
Graphing$WBContig_f<-factor(Graphing$WBContig, levels=factors)
p3<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + ggtitle("C") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[3]," (Mbp)",sep=""),y=paste(unique(Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",]$WBContig)," (Mbp)",sep=""))+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
#p3<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig_f~.,scales="free",space="free") + ggtitle("C") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[3]," (Mbp)",sep=""),y="WBContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",]$WBContig)
Graphing$WBContig_f<-factor(Graphing$WBContig, levels=factors)
p4<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + ggtitle("D") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[4]," (Mbp)",sep=""),y=paste(unique(Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",]$WBContig)," (Mbp)",sep=""))+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
#p4<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig~.,scales="free",space="free") + ggtitle("D") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[4]," (Mbp)",sep=""),y="WBContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))

factors<-unique(Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",]$WBContig)
Graphing$WBContig_f<-factor(Graphing$WBContig, levels=factors)
p5<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + ggtitle("E") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[5]," (Mbp)",sep=""),y=paste(unique(Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",]$WBContig)," (Mbp)",sep=""))+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))
#p5<-ggplot() + geom_segment(data=Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",],mapping=aes(x=BMStart,xend=BMEnd,y=WBStart,yend=WBEnd,color=PercIdentity,size=5)) + facet_grid(WBContig_f~.,scales="free",space="free") + ggtitle("E") + theme_bw()+theme(legend.text = element_text(size=15,face="bold"),strip.text.y = element_text(angle = 0),text = element_text(size=36),plot.title = element_text(size=60))+labs(x=paste(BMRelevants[5]," (Mbp)",sep=""),y="WBContigs (Mbp)")+scale_color_viridis_c(option = "inferno",rescaler = function(x, to = c(0, 1), from = NULL) {ifelse(x>65,scales::rescale(x,to = to,from = c(65, max(x, na.rm = TRUE))),0)}) + scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


mylegend<-g_legend(p1)
library(ggplot2)
library(gridExtra)
library(grid)

#png("/Users/jmattick/Documents/BP_NucmerMatch.png", width = 36, height = 16, units = 'in', res = 300)
#pdf("/Users/jmattick/Documents/WB_NucmerMatch.V3.pdf", width = 24, height = 86, useDingbats=FALSE)
pdf("/Users/jmattick/Documents/WB_NucmerMatch.pseudo.V3.pdf", width = 24, height = 86, useDingbats=FALSE)

gl<-list(p1 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p2 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p3 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p4 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),p5 + theme(legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"),plot.title = element_text(face = "bold")),mylegend)
grid.arrange(
  grobs = gl,
  widths = c(5,5,2),
  layout_matrix = rbind(c(1,2,6),
                        c(3,4,6),
                        c(5,5,6)))
ggarrange(p1,p2,p3,p4,p5, ncol=2,nrow=3,common.legend = TRUE,legend = "bottom")

dev.off()



#pdf("/Users/jmattick/Documents/All_SNVHist.pdf", useDingbats=FALSE)
pdf("/Users/jmattick/Documents/All_SNVHist.filtered.pdf", useDingbats=FALSE)

##WB
WuchContigsX<-Graphing[Graphing$BMContig == "Bm_v4_ChrX_scaffold_001",]$WBContig
WuchContigs1<-Graphing[Graphing$BMContig == "Bm_v4_Chr1_scaffold_001",]$WBContig
WuchContigs2<-Graphing[Graphing$BMContig == "Bm_v4_Chr2_contig_001",]$WBContig
WuchContigs3<-Graphing[Graphing$BMContig == "Bm_v4_Chr3_scaffold_001",]$WBContig
WuchContigs4<-Graphing[Graphing$BMContig == "Bm_v4_Chr4_scaffold_001",]$WBContig

#WBSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/W_bancrofti_jointcall.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))
WBSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/W_bancrofti_jointcall.filteredonly.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))

WBSNPs$Chromosome<-"Other"
WBSNPs[WBSNPs$CHROM %in% WuchContigsX,]$Chromosome<-"ChrX"
WBSNPs[WBSNPs$CHROM %in% WuchContigs1,]$Chromosome<-"Chr1"
WBSNPs[WBSNPs$CHROM %in% WuchContigs2,]$Chromosome<-"Chr2"
WBSNPs[WBSNPs$CHROM %in% WuchContigs3,]$Chromosome<-"Chr3"
WBSNPs[WBSNPs$CHROM %in% WuchContigs4,]$Chromosome<-"Chr4"
WBSNPs<-WBSNPs[WBSNPs$Chromosome != "Other",]

print(ggplot(WBSNPs, aes(x=Chromosome, y=PI)) + geom_boxplot(notch=TRUE) + ggtitle("WB SNV Distribution (10kb windows)") + theme_bw() + ylim(0,0.002))
print(mean(WBSNPs[WBSNPs$Chromosome == "ChrX",]$PI)/mean(WBSNPs[WBSNPs$Chromosome != "ChrX",]$PI))
print(median(WBSNPs[WBSNPs$Chromosome == "ChrX",]$PI)/median(WBSNPs[WBSNPs$Chromosome != "ChrX",]$PI))

#dev.off()


##BM
###BM Haploid SNPs
#HaploidSNPs <- readData("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/CombinedVCFs/Hap/", format = "VCF", include.unknown = TRUE, FAST = TRUE)
HaploidSNPs <- readData("/Users/jmattick/Documents/Hap/", format = "VCF", include.unknown = TRUE, FAST = TRUE)
HaploidSNPs <- diversity.stats(HaploidSNPs, pi = TRUE)
HaploidSNPs_sw <- sliding.window.transform(HaploidSNPs, width = 10000, jump = 10000, type = 2)
HaploidSNPs_sw <- diversity.stats(HaploidSNPs_sw, pi = TRUE)
HaploidSNPs_sw_div_sw <- HaploidSNPs_sw@nuc.diversity.within
HaploidSNPs_sw_div_sw <- HaploidSNPs_sw_div_sw/10000
position <- seq(from = 10001, to = 24943668, by = 10000)

#BMSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BM.final.jointcall.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))
BMSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BM.final.jointcall.filteredonly.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))

BMSNPs<-BMSNPs[grepl("Chr",BMSNPs$CHROM),]
BMSNPs$Chromosome<-gsub("_.*","",gsub("Bm_v4_","",BMSNPs$CHROM))
BMSNPs<-BMSNPs[BMSNPs$Chromosome != "ChrX",]
BMSNPs_ChrX<-data.frame("Bm_v4_ChrX_scaffold_001",position-10000,position,1,HaploidSNPs_sw_div_sw,"ChrX")
colnames(BMSNPs_ChrX)<-colnames(BMSNPs)
BMSNPs<-rbind(BMSNPs,BMSNPs_ChrX)
#pdf("/Users/jmattick/Documents/BP_SNVHist.pdf", useDingbats=FALSE)

print(ggplot(BMSNPs, aes(x=Chromosome, y=PI)) + geom_boxplot(notch=TRUE) + ggtitle("BM SNV Distribution (10kb windows)") + theme_bw())
print(mean(BMSNPs[BMSNPs$Chromosome == "ChrX",]$PI)/mean(BMSNPs[BMSNPs$Chromosome != "ChrX",]$PI))
print(median(BMSNPs[BMSNPs$Chromosome == "ChrX",]$PI)/median(BMSNPs[BMSNPs$Chromosome != "ChrX",]$PI))


##BP
#BPSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BP.final.jointcall.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))
HaploidSNPs <- readData("/Users/jmattick/Documents/Hap/", format = "VCF", include.unknown = TRUE, FAST = TRUE)
HaploidSNPs <- diversity.stats(HaploidSNPs, pi = TRUE)
HaploidSNPs_sw <- sliding.window.transform(HaploidSNPs, width = 10000, jump = 10000, type = 2)
HaploidSNPs_sw <- diversity.stats(HaploidSNPs_sw, pi = TRUE)
HaploidSNPs_sw_div_sw <- HaploidSNPs_sw@nuc.diversity.within
HaploidSNPs_sw_div_sw <- HaploidSNPs_sw_div_sw/10000
position <- seq(from = 10001, to = 24943668, by = 10000)

BPSNPs<-BPSNPs[grepl("Chr",BPSNPs$CHROM),]
BPSNPs$Chromosome<-gsub("_.","",gsub("BP_","",BPSNPs$CHROM))
#pdf("/Users/jmattick/Documents/BP_SNVHist.pdf", useDingbats=FALSE)
BPSNPs<-BPSNPs[BPSNPs$Chromosome != "ChrX",]
BPSNPs_ChrX<-data.frame("Bm_ChrX_c",position-10000,position,1,HaploidSNPs_sw_div_sw,"ChrX")
colnames(BPSNPs_ChrX)<-colnames(BPSNPs)
BPSNPs_ChrX<-rbind(BPSNPs,BPSNPs_ChrX)

print(ggplot(BPSNPs, aes(x=Chromosome, y=PI)) + geom_boxplot(notch=TRUE) + ggtitle("BP SNV Distribution (10kb windows)") + theme_bw())
print(mean(BPSNPs[BPSNPs$Chromosome == "ChrX",]$PI)/mean(BPSNPs[BPSNPs$Chromosome != "ChrX",]$PI))
print(median(BPSNPs[BPSNPs$Chromosome == "ChrX",]$PI)/median(BPSNPs[BPSNPs$Chromosome != "ChrX",]$PI))

#OV
#OVSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OV.final.jointcall.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))
OVSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OV.final.jointcall.filteredonly.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))

OVSNPs<-OVSNPs[grepl("OM",OVSNPs$CHROM),]
OVSNPs$Chromosome<-"Chr1"
OVSNPs[OVSNPs$CHROM == "OVOC_OM2",]$Chromosome<-"ChrX"
OVSNPs[OVSNPs$CHROM == "OVOC_OM5",]$Chromosome<-"ChrX"
OVSNPs[OVSNPs$CHROM == "OVOC_OM3",]$Chromosome<-"Chr2"
OVSNPs[OVSNPs$CHROM == "OVOC_OM4",]$Chromosome<-"Chr3"
HaploidSNPs <- readData("/Users/jmattick/Documents/Hap/", format = "VCF", include.unknown = TRUE, FAST = TRUE)
HaploidSNPs <- diversity.stats(HaploidSNPs, pi = TRUE)
HaploidSNPs_sw <- sliding.window.transform(HaploidSNPs, width = 10000, jump = 10000, type = 2)
HaploidSNPs_sw <- diversity.stats(HaploidSNPs_sw, pi = TRUE)
HaploidSNPs_sw_div_sw <- HaploidSNPs_sw@nuc.diversity.within
HaploidSNPs_sw_div_sw <- HaploidSNPs_sw_div_sw/10000
position <- seq(from = 10001, to = 24943668, by = 10000)

#pdf("/Users/jmattick/Documents/BP_SNVHist.pdf", useDingbats=FALSE)
OVSNPs<-OVSNPs[OVSNPs$Chromosome != "ChrX",]
OVSNPs_ChrX<-data.frame("OVOC_OM2",position-10000,position,1,HaploidSNPs_sw_div_sw,"ChrX")
colnames(OVSNPs_ChrX)<-colnames(OVSNPs)
OVSNPs_ChrX<-rbind(OVSNPs,OVSNPs_ChrX)
#pdf("/Users/jmattick/Documents/BP_SNVHist.pdf", useDingbats=FALSE)

print(ggplot(OVSNPs, aes(x=Chromosome, y=PI)) + geom_boxplot(notch=TRUE) + ggtitle("OV SNV Distribution (10kb windows)") + theme_bw())
print(mean(OVSNPs[OVSNPs$Chromosome == "ChrX",]$PI)/mean(OVSNPs[OVSNPs$Chromosome != "ChrX",]$PI))
print(median(OVSNPs[OVSNPs$Chromosome == "ChrX",]$PI)/median(OVSNPs[OVSNPs$Chromosome != "ChrX",]$PI))


dev.off()



####Distribution Plots
#pdf("/Users/jmattick/Documents/All_Pi.pdf",width = 20,height = 10,useDingbats=FALSE)
pdf("/Users/jmattick/Documents/All_Pi.filtered.pdf",width = 20,height = 10,useDingbats=FALSE)

pi_ymax<-max(c(BMSNPs$PI,BPSNPs$PI,OVSNPs$PI))
BMSNPs$Position<-BMSNPs$BIN_START/1000000
g1<-ggplot(BMSNPs[BMSNPs$Chromosome == "Chr1",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BM:A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g2<-ggplot(BMSNPs[BMSNPs$Chromosome == "Chr2",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BM:B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g3<-ggplot(BMSNPs[BMSNPs$Chromosome == "Chr3",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BM:C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g4<-ggplot(BMSNPs[BMSNPs$Chromosome == "Chr4",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BM:D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g5<-ggplot(BMSNPs[BMSNPs$Chromosome == "ChrX",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BM:E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))

gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))


BPSNPs$Position<-BPSNPs$BIN_START/1000000
g1<-ggplot(BPSNPs[BPSNPs$Chromosome == "Chr1",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BP:A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g2<-ggplot(BPSNPs[BPSNPs$Chromosome == "Chr2",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BP:B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g3<-ggplot(BPSNPs[BPSNPs$Chromosome == "Chr3",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BP:C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g4<-ggplot(BPSNPs[BPSNPs$Chromosome == "Chr4",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BP:D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g5<-ggplot(BPSNPs[BPSNPs$Chromosome == "ChrX",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("BP:E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))

gl<-list(g1,g2,g3,g4,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4),c(5,5))))

OVSNPs$Position<-OVSNPs$BIN_START/1000000
g1<-ggplot(OVSNPs[OVSNPs$Chromosome == "Chr1",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("OV:A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g2<-ggplot(OVSNPs[OVSNPs$Chromosome == "Chr2",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("OV:B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g3<-ggplot(OVSNPs[OVSNPs$Chromosome == "Chr3",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("OV:C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))
g5<-ggplot(OVSNPs[OVSNPs$Chromosome == "ChrX",]) + geom_point(aes(x=Position,y=PI),size=0.1) + theme_bw() + facet_grid(~CHROM,scales="free",space="free") + ggtitle("OV:E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,pi_ymax))

gl<-list(g1,g2,g3,g5)
print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,3),c(4,4))))


dev.off()

###OV Nigon Elements
Species<-"O_Volvulus"
n<-1
OVSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/OV.final.jointcall.filteredonly.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))
OVSNPs[OVSNPs$CHROM == 'OVOC_OM1a',]$BIN_START<-OVSNPs[OVSNPs$CHROM == 'OVOC_OM1a',]$BIN_START+(28345163*(n*45/100))
OVSNPs$CHROM<-gsub('OVOC_OM1a', 'Nigon-I', OVSNPs$CHROM)
OVSNPs[OVSNPs$BIN_START <= (28345163*(n*45/100)),]$CHROM <- gsub('OVOC_OM1b', 'Nigon-I', OVSNPs[OVSNPs$BIN_START <= (28345163*(n*45/100)),]$CHROM)
OVSNPs$CHROM <- gsub('OVOC_OM4', 'Nigon-II', OVSNPs$CHROM)
OVSNPs$CHROM <- gsub('OVOC_OM3', 'Nigon-III', OVSNPs$CHROM)
OVSNPs[OVSNPs$BIN_START <= (25485961*(n*55/100)),]$CHROM <- gsub('OVOC_OM2', 'Nigon-IV', OVSNPs[OVSNPs$BIN_START <= (25485961*(n*55/100)),]$CHROM)
OVSNPs[OVSNPs$BIN_START > (25485961*(n*55/100)),]$CHROM <- gsub('OVOC_OM2', 'Nigon-V', OVSNPs[OVSNPs$BIN_START > (25485961*(n*55/100)),]$CHROM)
OVSNPs[OVSNPs$CHROM == 'Nigon-V',]$BIN_START<-OVSNPs[OVSNPs$CHROM == 'Nigon-V',]$BIN_START-(25485961*(n*55/100))
OVSNPs[OVSNPs$CHROM == 'OVOC_OM5',]$BIN_START<-OVSNPs[OVSNPs$CHROM == 'OVOC_OM5',]$BIN_START+(25485961*(n*45/100))
OVSNPs$CHROM <- gsub('OVOC_OM5', 'Nigon-V', OVSNPs$CHROM)
OVSNPs[OVSNPs$BIN_START > (28345163*(n*45/100)),]$CHROM <- gsub('OVOC_OM1b', 'Nigon-X', OVSNPs[OVSNPs$BIN_START > (28345163*(n*45/100)),]$CHROM)
OVSNPs<-OVSNPs[grep("Nigon",OVSNPs$CHROM),]
colnames(OVSNPs)<-c("Nigon","Position","PosEnd","NVar","PI")
OVSNPs$Position<-OVSNPs$Position/1000000
pdf("/Users/jmattick/Documents/OV_Nigon_Pi.filtered.pdf")
ggplot(OVSNPs) + geom_point(aes(x=Position,y=PI,colour=Nigon),size=0.1) + facet_grid(~Nigon,scales="free",space="free") + ggtitle("O. volvulus PI by Nigon Element") + theme(axis.text.x=element_blank()) + theme(legend.position="none",plot.title = element_text(size=14))
dev.off()
write.table(OVSNPs,paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""),row.names = FALSE,col.names = FALSE,sep = "\t")

Species<-"B_malayi"

BMSNPs<-as.data.frame(read.delim("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BM.final.jointcall.filteredonly.pi.windowed.pi", sep="\t", header=T,stringsAsFactors = FALSE))
BMSNPs[BMSNPs$BIN_START <= (24943668*(n*50/100)),]$CHROM <- gsub('Bm_v4_ChrX_scaffold_001', 'Nigon-IV', BMSNPs[BMSNPs$BIN_START <= (24943668*(n*50/100)),]$CHROM)
BMSNPs[BMSNPs$BIN_START > (24943668*(n*50/100)),]$CHROM <- gsub('Bm_v4_ChrX_scaffold_001', 'Nigon-X', BMSNPs[BMSNPs$BIN_START > (24943668*(n*50/100)),]$CHROM)
BMSNPs$CHROM <- gsub('Bm_v4_Chr1_scaffold_001', 'Nigon-III', BMSNPs$CHROM)
BMSNPs$CHROM <- gsub('Bm_v4_Chr2_contig_001', 'Nigon-II', BMSNPs$CHROM)
BMSNPs$CHROM <- gsub('Bm_v4_Chr3_scaffold_001', 'Nigon-I', BMSNPs$CHROM)
BMSNPs$CHROM <- gsub('Bm_v4_Chr4_scaffold_001', 'Nigon-V', BMSNPs$CHROM)
BMSNPs<-BMSNPs[grep("Nigon",BMSNPs$CHROM),]
colnames(BMSNPs)<-c("Nigon","Position","PosEnd","NVar","PI")
BMSNPs$Position<-BMSNPs$Position/1000000
pdf("/Users/jmattick/Documents/BM_Nigon_Pi.filtered.pdf")
ggplot(BMSNPs) + geom_point(aes(x=Position,y=PI,colour=Nigon),size=0.1) + facet_grid(~Nigon,scales="free",space="free") + ggtitle("B. malayi PI by Nigon Element") + theme(axis.text.x=element_blank()) + theme(legend.position="none",plot.title = element_text(size=14))
dev.off()
write.table(BMSNPs,paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""),row.names = FALSE,col.names = FALSE,sep = "\t")

Species<-"Celegans"
PiArray<-list.files(path = paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",sep=""),pattern = ".pi")[1]
CESNPs<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",PiArray,sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
CESNPs$CHROM<-gsub('\\<I\\>', 'Nigon-I', CESNPs$CHROM)
CESNPs$CHROM<-gsub('\\<II\\>', 'Nigon-II', CESNPs$CHROM)
CESNPs$CHROM<-gsub('\\<III\\>', 'Nigon-III', CESNPs$CHROM)
CESNPs$CHROM<-gsub('\\<IV\\>', 'Nigon-IV', CESNPs$CHROM)
CESNPs$CHROM<-gsub('\\<V\\>', 'Nigon-V', CESNPs$CHROM)
CESNPs$CHROM<-gsub('\\<X\\>', 'Nigon-X', CESNPs$CHROM)
CESNPs<-CESNPs[grep("MtDNA",CESNPs$CHROM,invert=T),]
colnames(CESNPs)<-c("Nigon","Position","PosEnd","NVar","PI")
CESNPs$Position<-CESNPs$Position/1000000
pdf("/Users/jmattick/Documents/CE_Nigon_Pi.filtered.pdf")
ggplot(CESNPs) + geom_point(aes(x=Position,y=PI,colour=Nigon),size=0.1) + facet_grid(~Nigon,scales="free",space="free") + ggtitle("C. elegans PI by Nigon Element") + theme(axis.text.x=element_blank()) + theme(legend.position="none",plot.title = element_text(size=14))
dev.off()
write.table(CESNPs,paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""),row.names = FALSE,col.names = FALSE,sep = "\t")


Species<-"Dmelanogaster"
PiArray<-list.files(path = paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",sep=""),pattern = ".pi")[1]
DMSNPs<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",PiArray,sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
DMSNPs$CHROM<-gsub('NC_004354.4', 'ChrX', DMSNPs$CHROM)
DMSNPs$CHROM<-gsub('NT_033779.5', 'Chr2L', DMSNPs$CHROM)
DMSNPs$CHROM<-gsub('NT_033778.4', 'Chr2R', DMSNPs$CHROM)
DMSNPs$CHROM<-gsub('NT_037436.4', 'Chr3L', DMSNPs$CHROM)
DMSNPs$CHROM<-gsub('NT_033777.3', 'Chr3R', DMSNPs$CHROM)
DMSNPs$CHROM<-gsub('NC_004353.4', 'Chr4', DMSNPs$CHROM)
DMSNPs$CHROM<-gsub('NC_024512.1', 'ChrY', DMSNPs$CHROM)
DMSNPs<-DMSNPs[grep("Chr",DMSNPs$CHROM),]
colnames(DMSNPs)<-c("Chromosome","Position","PosEnd","NVar","PI")
DMSNPs$Position<-DMSNPs$Position/1000000
pdf("/Users/jmattick/Documents/DM_Chromosome_Pi.filtered.pdf")
ggplot(DMSNPs) + geom_point(aes(x=Position,y=PI,colour=Chromosome),size=0.1) + facet_grid(~Chromosome,scales="free",space="free") + ggtitle("D. melanogaster PI by Chromosome") + theme(axis.text.x=element_blank()) + theme(legend.position="none",plot.title = element_text(size=14))
dev.off()
write.table(DMSNPs,paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""),row.names = FALSE,col.names = FALSE,sep = "\t")


####Box plot and length plots
pdf("/Users/jmattick/Documents/All_Pi_Boxplots.filtered.pdf")

my_comparisons <- list( c("Nigon-D", "Nigon-A"), c("Nigon-D", "Nigon-B"), c("Nigon-D", "Nigon-C") )
my_comparisons2 <- list( c("ChrX", "Chr2L"), c("ChrX", "Chr2R"), c("ChrX", "Chr3L"), c("ChrX", "Chr3R"), c("ChrX", "Chr4") )
my_comparisons3 <- list( c("Nigon-X", "Nigon-A"), c("Nigon-X", "Nigon-B"), c("Nigon-X", "Nigon-C") )

NigonLengthFrame<-data.frame()
NigonElements<-c("Nigon-A","Nigon-B","Nigon-C","Nigon-D","Nigon-E","Nigon-X")
Species<-"O_Volvulus"
Graphing<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
colnames(Graphing)<-c("Nigon","Position","PosEnd","NVar","PI")
Graphing$Nigon<-gsub("Nigon-III","Nigon-C",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-II","Nigon-B",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-IV","Nigon-D",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-I","Nigon-A",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-V","Nigon-E",Graphing$Nigon)
Graphing$PI<-log(Graphing$PI)
ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("O. volvulus"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#OVplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("O. volvulus"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.005) + stat_compare_means(comparisons = my_comparisons,label.y = c(0.004,0.0045,0.005))
OVplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("O. volvulus"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none")
subNigonLength<-data.frame()
for (a in 1:length(NigonElements))
{
  NigonLength<-length(Graphing[Graphing$Nigon == NigonElements[a],]$Position)*10000
  subsubNigonLength<-data.frame(Species,NigonElements[a],NigonLength)
  subNigonLength<-rbind(subNigonLength,subsubNigonLength)
}
NigonLengthFrame<-rbind(NigonLengthFrame,subNigonLength)

Species<-"B_malayi"
Graphing<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
colnames(Graphing)<-c("Nigon","Position","PosEnd","NVar","PI")
Graphing$Nigon<-gsub("Nigon-III","Nigon-C",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-II","Nigon-B",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-IV","Nigon-D",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-I","Nigon-A",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-V","Nigon-E",Graphing$Nigon)
Graphing$PI<-log(Graphing$PI)
ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("O. volvulus"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#BMplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("B. malayi"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.005) + stat_compare_means(comparisons = my_comparisons,label.y = c(0.004,0.0045,0.005))
BMplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("B. malayi"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none")

subNigonLength<-data.frame()
for (a in 1:length(NigonElements))
{
  NigonLength<-length(Graphing[Graphing$Nigon == NigonElements[a],]$Position)*10000
  subsubNigonLength<-data.frame(Species,NigonElements[a],NigonLength)
  subNigonLength<-rbind(subNigonLength,subsubNigonLength)
}
NigonLengthFrame<-rbind(NigonLengthFrame,subNigonLength)

Species<-"Celegans"
Graphing<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
colnames(Graphing)<-c("Nigon","Position","PosEnd","NVar","PI")
Graphing$Nigon<-gsub("Nigon-III","Nigon-C",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-II","Nigon-B",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-IV","Nigon-D",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-I","Nigon-A",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-V","Nigon-E",Graphing$Nigon)
if(1==2){
  p1<-ggplot() + geom_histogram(data=Graphing[Graphing$Nigon == "Nigon-D",], aes(PI),binwidth=0.0001,alpha = .5) + theme_bw() + ggtitle("Nigon D for C. elegans")+xlim(0,0.005)+geom_vline(xintercept=mean(Graphing[Graphing$Nigon == "Nigon-D",]$PI))
  p2<-ggplot() + geom_histogram(data=Graphing[Graphing$Nigon == "Nigon-C",], aes(PI),binwidth=0.0001,alpha = .5) + theme_bw() + ggtitle("Nigon C for C. elegans")+xlim(0,0.005)+geom_vline(xintercept=mean(Graphing[Graphing$Nigon == "Nigon-C",]$PI))
  p3<-ggplot() + geom_histogram(data=Graphing[Graphing$Nigon == "Nigon-B",], aes(PI),binwidth=0.0001,alpha = .5) + theme_bw() + ggtitle("Nigon B for C. elegans")+xlim(0,0.005)+geom_vline(xintercept=mean(Graphing[Graphing$Nigon == "Nigon-B",]$PI))
  p4<-ggplot() + geom_histogram(data=Graphing[Graphing$Nigon == "Nigon-A",], aes(PI),binwidth=0.0001,alpha = .5) + theme_bw() + ggtitle("Nigon A for C. elegans")+xlim(0,0.005)+geom_vline(xintercept=mean(Graphing[Graphing$Nigon == "Nigon-A",]$PI))
  pdf("/Users/jmattick/Documents/CelegansDistributions.pdf",width=10,height=15)
  gl<-list(p1,p2,p3,p4)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4))))
  dev.off()
}
ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("C. elegans"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
Graphing$PI<-log(Graphing$PI)
#CEplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("C. elegans"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.005) + stat_compare_means(comparisons = my_comparisons3,label.y = c(0.004,0.0045,0.005))
CEplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("C. elegans"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + stat_compare_means(comparisons = my_comparisons3,label.y = c(max(Graphing$PI)-0.5,max(Graphing$PI),max(Graphing$PI)+0.5))

if(1==2){
###Outlier Removal
  CEplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("C. elegans"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.005) + stat_compare_means(comparisons = my_comparisons,label.y = c(0.004,0.0045,0.005))
  

  outliers <- boxplot(Graphing[Graphing$Nigon == "Nigon-A",]$PI, plot=FALSE)$out
  Graphing2<-Graphing[!(Graphing$PI %in% outliers) & Graphing$Nigon == "Nigon-A",]
  
  outliers <- boxplot(Graphing[Graphing$Nigon == "Nigon-B",]$PI, plot=FALSE)$out
  Graphing2<-rbind(Graphing2,Graphing[!(Graphing$PI %in% outliers) & Graphing$Nigon == "Nigon-B",])
  
  outliers <- boxplot(Graphing[Graphing$Nigon == "Nigon-C",]$PI, plot=FALSE)$out
  Graphing2<-rbind(Graphing2,Graphing[!(Graphing$PI %in% outliers) & Graphing$Nigon == "Nigon-C",])
  
  outliers <- boxplot(Graphing[Graphing$Nigon == "Nigon-D",]$PI, plot=FALSE)$out
  Graphing2<-rbind(Graphing2,Graphing[!(Graphing$PI %in% outliers) & Graphing$Nigon == "Nigon-D",])
  
  outliers <- boxplot(Graphing[Graphing$Nigon == "Nigon-E",]$PI, plot=FALSE)$out
  Graphing2<-rbind(Graphing2,Graphing[!(Graphing$PI %in% outliers) & Graphing$Nigon == "Nigon-E",])
  
  outliers <- boxplot(Graphing[Graphing$Nigon == "Nigon-X",]$PI, plot=FALSE)$out
  Graphing2<-rbind(Graphing2,Graphing[!(Graphing$PI %in% outliers) & Graphing$Nigon == "Nigon-X",])
  
  p1<-ggplot() + geom_histogram(data=Graphing2[Graphing2$Nigon == "Nigon-D",], aes(PI),binwidth=0.0001,alpha = .5) + theme_bw() + ggtitle("Nigon D for C. elegans")+xlim(0,0.005)+geom_vline(xintercept=mean(Graphing2[Graphing2$Nigon == "Nigon-D",]$PI))
  p2<-ggplot() + geom_histogram(data=Graphing2[Graphing2$Nigon == "Nigon-C",], aes(PI),binwidth=0.0001,alpha = .5) + theme_bw() + ggtitle("Nigon C for C. elegans")+xlim(0,0.005)+geom_vline(xintercept=mean(Graphing2[Graphing2$Nigon == "Nigon-C",]$PI))
  p3<-ggplot() + geom_histogram(data=Graphing2[Graphing2$Nigon == "Nigon-B",], aes(PI),binwidth=0.0001,alpha = .5) + theme_bw() + ggtitle("Nigon B for C. elegans")+xlim(0,0.005)+geom_vline(xintercept=mean(Graphing2[Graphing2$Nigon == "Nigon-B",]$PI))
  p4<-ggplot() + geom_histogram(data=Graphing2[Graphing2$Nigon == "Nigon-A",], aes(PI),binwidth=0.0001,alpha = .5) + theme_bw() + ggtitle("Nigon A for C. elegans")+xlim(0,0.005)+geom_vline(xintercept=mean(Graphing2[Graphing2$Nigon == "Nigon-A",]$PI))
  pdf("/Users/jmattick/Documents/CelegansDistributionsOutlierRemoval.pdf",width=10,height=15)
  gl<-list(p1,p2,p3,p4)
  print(grid.arrange(grobs = gl, widths = c(5,5),layout_matrix = rbind(c(1,2),c(3,4))))
  dev.off()
  
  CEplot<-ggplot(Graphing2, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("C. elegans"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.005) + stat_compare_means(comparisons = my_comparisons,label.y = c(0.004,0.0045,0.005))
  print(CEplot)
}
dev.off()
subNigonLength<-data.frame()
for (a in 1:length(NigonElements))
{
  NigonLength<-length(Graphing[Graphing$Nigon == NigonElements[a],]$Position)*10000
  subsubNigonLength<-data.frame(Species,NigonElements[a],NigonLength)
  subNigonLength<-rbind(subNigonLength,subsubNigonLength)
}
NigonLengthFrame<-rbind(NigonLengthFrame,subNigonLength)

Species<-"Dmelanogaster"
Graphing<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
colnames(Graphing)<-c("Chromosome","Position","PosEnd","NVar","PI")
Graphing$PI<-log(Graphing$PI)
ggplot(Graphing, aes(x=Chromosome, y=PI,color=Chromosome)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("D. melanogaster"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#Dmelplot<-ggplot(Graphing, aes(x=Chromosome, y=PI,color=Chromosome)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("D. melanogaster"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.025) + stat_compare_means(comparisons = my_comparisons2,label.y = c(0.021,0.022,0.023,0.024,0.025))
Dmelplot<-ggplot(Graphing, aes(x=Chromosome, y=PI,color=Chromosome)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("D. melanogaster"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none")

Species<-"Wbancrofti"
Graphing<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
colnames(Graphing)<-c("Position","PI","Nigon")
Graphing$Nigon<-gsub("Nigon-III","Nigon-C",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-II","Nigon-B",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-IV","Nigon-D",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-I","Nigon-A",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-V","Nigon-E",Graphing$Nigon)
Graphing$PI<-log(Graphing$PI)
ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("W. bancrofti"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#WBplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("W. bancrofti"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.0025) + stat_compare_means(comparisons = my_comparisons,label.y = c(0.002,0.00225,0.0025))
WBplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("W. bancrofti"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none")
subNigonLength<-data.frame()
for (a in 1:length(NigonElements))
{
  NigonLength<-length(Graphing[Graphing$Nigon == NigonElements[a],]$Position)*10000
  subsubNigonLength<-data.frame(Species,NigonElements[a],NigonLength)
  subNigonLength<-rbind(subNigonLength,subsubNigonLength)
}
NigonLengthFrame<-rbind(NigonLengthFrame,subNigonLength)

Species<-"Dirofilariaimmitis"
Graphing<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
colnames(Graphing)<-c("Position","PI","Nigon")
Graphing$Nigon<-gsub("Nigon-III","Nigon-C",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-II","Nigon-B",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-IV","Nigon-D",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-I","Nigon-A",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-V","Nigon-E",Graphing$Nigon)
Graphing$PI<-log(Graphing$PI)
ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("D. immitis"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#DIplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("D. immitis"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.0025) + stat_compare_means(comparisons = my_comparisons,label.y = c(0.002,0.00225,0.0025))
DIplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("D. immitis"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none")
subNigonLength<-data.frame()
for (a in 1:length(NigonElements))
{
  NigonLength<-length(Graphing[Graphing$Nigon == NigonElements[a],]$Position)*10000
  subsubNigonLength<-data.frame(Species,NigonElements[a],NigonLength)
  subNigonLength<-rbind(subNigonLength,subsubNigonLength)
}
NigonLengthFrame<-rbind(NigonLengthFrame,subNigonLength)

Species<-"LoaLoa"
Graphing<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
colnames(Graphing)<-c("Position","PI","Nigon")
Graphing$Nigon<-gsub("Nigon-III","Nigon-C",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-II","Nigon-B",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-IV","Nigon-D",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-I","Nigon-A",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-V","Nigon-E",Graphing$Nigon)
Graphing$PI<-log(Graphing$PI)
ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("L. loa"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#LLplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("L. loa"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.002) + stat_compare_means(comparisons = my_comparisons,label.y = c(0.0015,0.00175,0.002))
LLplot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("L. loa"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none")
subNigonLength<-data.frame()
for (a in 1:length(NigonElements))
{
  NigonLength<-length(Graphing[Graphing$Nigon == NigonElements[a],]$Position)*10000
  subsubNigonLength<-data.frame(Species,NigonElements[a],NigonLength)
  subNigonLength<-rbind(subNigonLength,subsubNigonLength)
}
NigonLengthFrame<-rbind(NigonLengthFrame,subNigonLength)

Species<-"B_pahangi"
Graphing<-as.data.frame(read.delim(paste("/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/NigonElement_Pi/",Species,"/",Species,".PIGraphing.tsv",sep=""), sep="\t", header=T,stringsAsFactors = FALSE))
colnames(Graphing)<-c("Position","PI","Nigon")
Graphing$Nigon<-gsub("Nigon-III","Nigon-C",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-II","Nigon-B",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-IV","Nigon-D",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-I","Nigon-A",Graphing$Nigon)
Graphing$Nigon<-gsub("Nigon-V","Nigon-E",Graphing$Nigon)
Graphing$PI<-log(Graphing$PI)
ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("B. pahangi"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
#BPPlot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("B. pahangi"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none") + ylim(0,0.02) + stat_compare_means(comparisons = my_comparisons,label.y = c(0.015,0.0175,0.02))
BPPlot<-ggplot(Graphing, aes(x=Nigon, y=PI,color=Nigon)) + geom_boxplot(notch=TRUE) + ggtitle(expression(paste("PI of ",italic("B. pahangi"),sep=""))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+theme(legend.position = "none")
subNigonLength<-data.frame()
for (a in 1:length(NigonElements))
{
  NigonLength<-length(Graphing[Graphing$Nigon == NigonElements[a],]$Position)*10000
  subsubNigonLength<-data.frame(Species,NigonElements[a],NigonLength)
  subNigonLength<-rbind(subNigonLength,subsubNigonLength)
}
NigonLengthFrame<-rbind(NigonLengthFrame,subNigonLength)

dev.off()
pdf("/Users/jmattick/Documents/All_Pi_Boxplots.combined.filtered.pdf",width=10,height=15)
gl<-list(BMplot,BPPlot,WBplot,LLplot,DIplot,OVplot,CEplot,Dmelplot)
print(grid.arrange(grobs = gl, widths = c(5,5,5),layout_matrix = rbind(c(1,2,3),c(4,5,6),c(7,8,9))))
dev.off()


pdf("/Users/jmattick/Documents/All_Nigon_Lengths.pdf")
colnames(NigonLengthFrame)<-c("Species","Nigon","Length")
ggplot(NigonLengthFrame, aes(fill=Nigon, y=Length, x=Species)) + geom_bar(position="dodge", stat="identity") + theme_bw() + ggtitle("Nigon element lengths for all species")
dev.off()