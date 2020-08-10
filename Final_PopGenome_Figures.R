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
      colnames(PartialGraphing)<-c("Position","Pi","Chr","subChr")
      for (d in 1:ceiling(ChrSize/10000))
      {
        if(d<ceiling(ChrSize/10000)){
          pi<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/(20000*SampleFactor)
          PartialGraphing[d,]$Position<-(d*10000)
          PartialGraphing[d,]$Pi<-pi
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-subChromosomes[c]
        } else {
          pi<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= ChrSize & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/((ChrSize-((d*10000)-10000))*SampleFactor)
          PartialGraphing[d,]$Position<-ChrSize
          PartialGraphing[d,]$Pi<-pi
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-subChromosomes[c]
        }
      }
      Graphing<-rbind(Graphing,PartialGraphing)
    }
  }
  pi_ymax<-max(Graphing$Pi)
  TotalMax<-c(TotalMax,pi_ymax)
  Graphing$Position<-Graphing$Position/1000000
  g1<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g2<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g3<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g4<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g5<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  
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
TotalMaxPreCalc<-0.01


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
      colnames(PartialGraphing)<-c("Position","Pi","Chr","subChr")
      for (d in 1:ceiling(ChrSize/10000))
      {
        if(d<ceiling(ChrSize/10000)){
          pi<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/(20000*SampleFactor)
          PartialGraphing[d,]$Position<-(d*10000)
          PartialGraphing[d,]$Pi<-pi
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-subChromosomes[c]
        } else {
          pi<-(length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= ChrSize & subSNPs$Zyg == "het",]$Pos)+(2*length(subSNPs[subSNPs$Pos >= ((d*10000)-9999) & subSNPs$Pos <= (d*10000) & subSNPs$Zyg == "hom",]$Pos)))/((ChrSize-((d*10000)-10000))*SampleFactor)
          PartialGraphing[d,]$Position<-ChrSize
          PartialGraphing[d,]$Pi<-pi
          PartialGraphing[d,]$Chr<-ChrName
          PartialGraphing[d,]$subChr<-subChromosomes[c]
        }
      }
      Graphing<-rbind(Graphing,PartialGraphing)
    }
  }
  pi_ymax<-max(Graphing$Pi)
  TotalMax<-c(TotalMax,pi_ymax)
  Graphing$Position<-Graphing$Position/1000000
  #g6<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("A") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  #g7<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("B") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  #g8<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("C") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  #g9<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("D") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  #g10<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle("E") + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  
  
  g6<-ggplot(Graphing[Graphing$Chr == "Chr1",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":A",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g7<-ggplot(Graphing[Graphing$Chr == "Chr2",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":B",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g8<-ggplot(Graphing[Graphing$Chr == "Chr3",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":C",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g9<-ggplot(Graphing[Graphing$Chr == "Chr4",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":D",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  g10<-ggplot(Graphing[Graphing$Chr == "ChrX",]) + geom_point(aes(x=Position,y=Pi),size=0.1) + theme_bw() + facet_grid(~subChr,scales="free",space="free") + ggtitle(paste(BMGroups[a],":E",sep="")) + theme(legend.position="none",plot.title = element_text(size=14,face = "bold"),plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))+xlab("Position (Mb)")+scale_x_continuous(breaks = trans_breaks(identity, identity, n = 3),expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits=c(0,TotalMaxPreCalc))
  
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
Bp.vec.subset$SampleGroup<-c("MultiAF",rep("AM",3),rep("Clinical",3),rep("FR3",4))
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

CommonBuscos<-BuscoFreq[BuscoFreq$Freq == 5,]

Bp.busco<-Bp.busco[Bp.busco$V1 %in% CommonBuscos$Var1,]
Bm.busco<-Bm.busco[Bm.busco$V1 %in% CommonBuscos$Var1,]
Wb.busco<-Wb.busco[Wb.busco$V1 %in% CommonBuscos$Var1,]
Bt.busco<-Bt.busco[Bt.busco$V1 %in% CommonBuscos$Var1,]
Ov.busco<-Ov.busco[Ov.busco$V1 %in% CommonBuscos$Var1,]


FinalFrame<-as.data.frame(cbind(Bp.busco$V1,Bp.busco$V3,Bm.busco$V3,Wb.busco$V3,Bt.busco$V3,Ov.busco$V3))
FinalFrame$V7<-gsub("_.","",gsub("BP_","",FinalFrame$V2))
FinalFrame$V8<-gsub("_scaffold_001|_contig_001","",gsub("Bm_v4_","",FinalFrame$V3))
FinalFrame$V9<-"Same"
FinalFrame[FinalFrame$V7 != FinalFrame$V8,]$V9<-"Shifted"

write.table(FinalFrame,"/Volumes/projects-t3/EBMAL/SNP_MALE_GATK_Best_Practices/BUSCO_MultiSpecies/ROutput.busco.multispecies.tsv",row.names = FALSE,col.names = F,sep = "\t",quote = FALSE)


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


