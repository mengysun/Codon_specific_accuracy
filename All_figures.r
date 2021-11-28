#load packages
library(Biostrings)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(extrafont)
library(reshape2)
library(ggtree)
library(phytools)
library(R.utils)

#load tables
Relative_error_tab<-read.table(file="Table1",sep="\t",header=TRUE)
Bacteria_meta<-read.table(file="Table2",sep="\t",header=TRUE)
Bacteria_OR<-read.table(file="Table3",sep="\t",header=TRUE)
phy<-read.tree(file="tree.nwk")
Translation_evolution_SimResults<-read.table(file="Table4",sep="\t",header=TRUE)
bacteria_codon_counts<-read.table("Table5",sep="\t",header=TRUE)
Archaea_meta<-read.table(file="Table6",sep="\t",header=TRUE)
Archaea_OR<-read.table(file="Table7",sep="\t",header=TRUE)
Eukaryotic_OR<-read.table(file="Table8",sep="\t",header=TRUE)

#Figure plots
#Fig 1a
aa_list<-unique(as.character(Relative_error_tab$aa))
p = list()
for(i in 1:length(aa_list)){
  codon_tab_error<-Relative_error_tab%>%
    filter(aa==aa_list[i])
  codon_tab_error_ordered<-codon_tab_error[order(codon_tab_error$RSCU),]
  codon_tab_error$codons<-factor(codon_tab_error$Origin_codon,levels=unique(codon_tab_error_ordered$Origin_codon))
  amino_acid<-AMINO_ACID_CODE[aa_list[i]]
  g<-ggplot(codon_tab_error,aes(codons))+geom_bar(aes(y=RSCU),stat="identity",fill="gray",alpha=1,color="black")
  g<-g+geom_point(aes(y=RMR),size=5,color="#D56C55")
  g<-g+scale_y_continuous(limits = c(-0.1,3.5))
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(family="Lucida Sans",size=20,hjust=0.5,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+geom_errorbar(aes(ymin=RMR-sd_of_RMR,ymax=RMR+sd_of_RMR),width=0.02,position = position_dodge(0.9))
  g<-g+theme(axis.text.x =element_text(family="Helvectica",size=15,face="bold"))
  g<-g+theme(axis.text.y =element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
             panel.background=element_rect(color="black",size=2,fill="white"),axis.line=element_line(colour="black"))
  p[[i]]<-g
}
png(filename="Fig1a.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

#Fig 1b
png(file="./Fig1b.png",units="in",width=5,height=5,res=600)
ggplot(Relative_error_tab,aes(x=RSCU,y=RMR))+
  geom_point(color="#738FC1",alpha=1,size=3)+
  geom_errorbar(aes(ymin=RMR-sd_of_RMR,ymax=RMR+sd_of_RMR),size=0.1,width=0.001,position = position_dodge(0.9))+
  geom_smooth(method="lm",se=FALSE,color="#D56C55")+
  xlab("RSCU")+
  ylab("RMR")+
  theme(axis.text = element_text(family="Helvectica",size=18))+
  theme(axis.title=element_text(family="Lucida Sans",size=15,face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="#FFEDCE"),axis.line=element_line(colour="black"))
dev.off()

#Fig 2b
png(file="./Fig2b.png",units="in",width=5,height=5,res=600)
ggplot(Relative_error_tab,aes(x=OR,y=RMR))+
  geom_point(color="#738FC1",alpha=1,size=3)+
  geom_smooth(method="lm",se=FALSE,color="#D56C55")+
  geom_errorbar(aes(ymin=RMR-sd_of_RMR,ymax=RMR+sd_of_RMR),size=0.1,width=0.001,position = position_dodge(0.9))+
  geom_errorbar(aes(xmin=OR-sd_of_OR,xmax=OR+sd_of_OR),size=0.1,width=0.001,position = position_dodge(0.9))+
  xlab("OR")+
  ylab("RMR")+
  theme(axis.text = element_text(family="Helvectica",size=15))+
  theme(axis.title=element_text(family="Lucida Sans",size=15,face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="#FFEDCE"),axis.line=element_line(colour="black"))
dev.off()

#Fig 2c
png(filename="./Fig2c.png",width=8,height=5,unit="in",res=600)
ggplot(Bacteria_meta,aes(x=cor_OR_RSCU))+
  theme_bw()+
  geom_histogram(colour="black",fill="#738FC1",aes(y=..count../sum(..count..)),alpha=0.8)+
  xlab(label="Correlations between RSCU and OR")+
  ylab(label="Frequency")+
  theme(axis.text.x=element_text(size=15,family="Helvetica",color="black"))+
  theme(axis.text.y=element_text(size=15,family="Helvetica",color="black"))+
  theme(axis.title.x=element_text(size=15,family="Lucida Sans",color="black",face="bold"))+
  theme(axis.title.y=element_text(size=15,family="Lucida Sans",color="black",face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#FFEDCE"), axis.line = element_line(colour = "black"))
dev.off()

#Fig 2d
unique_order<-unique(as.character(Bacteria_meta$order))
order_unique_idx<-match(unique_order,Bacteria_meta$order)
keep_genome<-as.character(Bacteria_meta$taxa[order_unique_idx])
keep_genome_idx<-match(keep_genome,names(Bacteria_OR))
Bacteria_OR_CAT<-Bacteria_OR[21,keep_genome_idx]
species<-names(Bacteria_OR_CAT)
pruned.phy<-drop.tip(phy,phy$tip.label[-match(species,phy$tip.label)])
species_idx<-match(pruned.phy$tip.label,species)
traits<-data.frame(log(c(as.numeric(as.character(Bacteria_OR_CAT[1,]))[species_idx])))
names(traits)<-"colors_plot"
row.names(traits)<-pruned.phy$tip.label
png(filename="./Fig2d.png",width=5,height=5,unit="in",res=600)
p<-ggtree(pruned.phy,layout='circular')
gheatmap(p, traits, offset=0.2, width=0.2,font.size=2,colnames = FALSE,low="red",high="green",legend_title = "ln(OR)")
dev.off()

#Fig 2e
aa_list<-unique(as.character(Bacteria_OR$aa))
p = list()
for(i in 1:9){
  odds_ratio_tmp<-Bacteria_OR%>%
    filter(aa==aa_list[i])
  amino_acid<-AMINO_ACID_CODE[aa_list[i]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],1197))
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:1199)])))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,log(OR)))+geom_violin(fill="gray",color="black")
  g<-g+ylim(-2,2)
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=13,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  p[[i]]<-g
}
png(filename="Fig2e1.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

p = list()
for(i in 1:9){
  odds_ratio_tmp<-Bacteria_OR%>%
    filter(aa==aa_list[i+9])
  amino_acid<-AMINO_ACID_CODE[aa_list[i+9]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],1197))
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:1199)])))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,log(OR)))+geom_violin(fill="gray",color="black")
  g<-g+ylim(-2,2)
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=13,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  p[[i]]<-g
}
png(filename="Fig2e2.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

#Fig 3a
png(file="./Fig3a.png",units="in",width=6,height=6,res=600)
ggplot(Relative_error_tab,aes(x=RRcnc_fromExp,y=RMR))+
  geom_point(color="#738FC1",alpha=1,size=3)+
  geom_smooth(method="lm",se=FALSE,color="#D56C55")+
  geom_errorbar(aes(ymin=RMR-sd_of_RMR,ymax=RMR+sd_of_RMR),size=0.1,width=0.001,position = position_dodge(0.9))+
  xlab("RR[c/nc]")+
  ylab("RMR")+
  theme(axis.text = element_text(family="Helvectica",size=15))+
  theme(axis.title=element_text(family="Lucida Sans",size=15,face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="#FFEDCE"),axis.line=element_line(colour="black"))
dev.off()

#Fig 3b
png(file="./Fig3b.png",units="in",width=6,height=6,res=600)
ggplot(Relative_error_tab,aes(x=RRcnc_fromExp,y=OR))+
  geom_point(color="#738FC1",alpha=1,size=3)+
  geom_smooth(method="lm",se=FALSE,color="#D56C55")+
  xlab("RR[c/nc]")+
  ylab("OR")+
  geom_errorbar(aes(ymin=OR-sd_of_OR,ymax=OR+sd_of_OR),size=0.1,width=0.001,position = position_dodge(0.9))+
  theme(axis.text = element_text(family="Helvectica",size=15))+
  theme(axis.title=element_text(family="Lucida Sans",size=15,face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="#FFEDCE"),axis.line=element_line(colour="black"))
dev.off()

#Fig 3c
png(file="./Fig3c.png",units="in",width=6,height=6,res=600)
ggplot(Relative_error_tab,aes(x=RRcnc_fromCopy,y=RMR))+
  geom_point(color="#738FC1",alpha=1,size=3)+
  geom_smooth(method="lm",se=FALSE,color="#D56C55")+
  geom_errorbar(aes(ymin=RMR-sd_of_RMR,ymax=RMR+sd_of_RMR),size=0.1,width=0.001,position = position_dodge(0.9))+
  xlab("RR[c/nc]")+
  ylab("RMR")+
  theme(axis.text = element_text(family="Helvectica",size=15))+
  theme(axis.title=element_text(family="Lucida Sans",size=15,face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="#FFEDCE"),axis.line=element_line(colour="black"))
dev.off()

#Fig 3d
png(filename="./Fig3d.png",width=6,height=6,unit="in",res=600)
Bacteria_meta_filtered_by_tRNA<-Bacteria_meta%>%
  filter(total_tRNA_counts>80)
ggplot(Bacteria_meta_filtered_by_tRNA,aes(x=cor_OR_RRcncfromCopy))+
  theme_bw()+
  geom_histogram(colour="black",fill="#738FC1",aes(y=..count../sum(..count..)),alpha=0.8)+
  xlab(bquote('Correlations'~'between'~RR[c/nc]~and~OR))+
  ylab(label="Frequency")+
  theme(axis.text.x=element_text(size=15,family="Helvetica",color="black"))+
  theme(axis.text.y=element_text(size=15,family="Helvetica",color="black"))+
  theme(axis.title.x=element_text(size=15,family="Lucida Sans",color="black",face="bold"))+
  theme(axis.title.y=element_text(size=15,family="Lucida Sans",color="black",face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#FFEDCE"), axis.line = element_line(colour = "black"))
dev.off()

#Fig 4c
png(filename="Fig4c.png",units="in",width=10,height=6,res=600)
Translation_evolution_SimResults$pop_size<-as.character(Translation_evolution_SimResults$pop_size)
ggplot(Translation_evolution_SimResults,aes(fill=with_mistranslation,y=translation_efficiency,x=pop_size))+
  geom_boxplot(position="dodge", alpha=0.5)+
  theme(axis.text = element_text(family="Helvectica",size=18))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white"),axis.line=element_line(colour="black"))+
  scale_fill_manual(values=c("white","pink"))
dev.off()

#Fig 4d
png(filename="Fig4d.png",units="in",width=10,height=6,res=600)
ggplot(Translation_evolution_SimResults,aes(fill=with_mistranslation,y=RSCU_diff,x=pop_size))+
  geom_boxplot(position="dodge", alpha=0.5)+
  theme(axis.text = element_text(family="Helvectica",size=18))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="white"),axis.line=element_line(colour="black"))+
  scale_fill_manual(values=c("white","pink"))
dev.off()

#Fig S1a
aa_list<-unique(as.character(Bacteria_OR$aa))
p = list()
for(i in 1:9){
  odds_ratio_tmp<-Bacteria_OR%>%
    filter(aa==aa_list[i])
  amino_acid<-AMINO_ACID_CODE[aa_list[i]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
   codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],1197))
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:1199)])))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab<-or_for_drawTab%>%
    group_by(codons)%>%
    dplyr::summarise(f=sum(OR>1)/length(OR))
  print(as.character(or_for_drawTab$f))
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,f))+geom_bar(stat="identity",fill="gray",color="black")
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=12,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=15))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  g<-g+scale_y_continuous(limits=c(0,1))
  p[[i]]<-g
}
png(filename="FigS1a1.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

p = list()
for(i in 1:9){
  odds_ratio_tmp<-Bacteria_OR%>%
    filter(aa==aa_list[i+9])
  amino_acid<-AMINO_ACID_CODE[aa_list[i+9]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],1197))
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:1199)])))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab<-or_for_drawTab%>%
    group_by(codons)%>%
    dplyr::summarise(f=sum(OR>1)/length(OR))
  print(as.character(or_for_drawTab$f))
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,f))+geom_bar(stat="identity",fill="gray",color="black")
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=12,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=15))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  g<-g+scale_y_continuous(limits=c(0,1))
  p[[i]]<-g
}
png(filename="FigS1a2.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

#Fig S2b
aa_list<-unique(as.character(Bacteria_OR$aa))
p = list()
for(i in 1:9){
  odds_ratio_tmp<-Bacteria_OR%>%
    filter(aa==aa_list[i])
  dat_codonC_tmp<-bacteria_codon_counts%>%
    filter(aa==aa_list[i])
  amino_acid<-AMINO_ACID_CODE[aa_list[i]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  #get the cols to keep.
  s_keep_dat<-data.frame(matrix(ncol=1197,nrow=length(codon_tmp)))
  for(k in 1:length(codon_tmp)){
    s_keep_dat[k,]<-(as.numeric(as.character(dat_codonC_tmp[k,3:1199]))>1000)
  }
  s_keep<-(colSums(s_keep_dat)>=length(codon_tmp))
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:1199)]))[s_keep])
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],sum(s_keep)))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,log(OR)))+geom_violin(fill="gray",color="black")
  g<-g+ylim(-2,2)
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=13,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  p[[i]]<-g
}
png(filename="FigS1b1.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

p = list()
for(i in 1:9){
  odds_ratio_tmp<-Bacteria_OR%>%
    filter(aa==aa_list[i+9])
  dat_codonC_tmp<-bacteria_codon_counts%>%
    filter(aa==aa_list[i+9])
  amino_acid<-AMINO_ACID_CODE[aa_list[i+9]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  s_keep_dat<-data.frame(matrix(ncol=1197,nrow=length(codon_tmp)))
  for(k in 1:length(codon_tmp)){
    s_keep_dat[k,]<-(as.numeric(as.character(dat_codonC_tmp[k,3:1199]))>1000)
  }
  s_keep<-(colSums(s_keep_dat)>=length(codon_tmp))
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:1199)]))[s_keep])
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],sum(s_keep)))
    
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,log(OR)))+geom_violin(fill="gray",color="black")
  g<-g+ylim(-2,2)
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=13,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  p[[i]]<-g
}
png(filename="FigS1b2.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

#Fig S1c
Bacteria_OR_StrongSelection<-Bacteria_OR[,c(1:2,c(3:1199)[Bacteria_meta$cor_OR_RSCU>0.5])]
aa_list<-unique(as.character(Bacteria_OR_StrongSelection$aa))
p = list()
for(i in 1:9){
  odds_ratio_tmp<-Bacteria_OR_StrongSelection%>%
    filter(aa==aa_list[i])
  amino_acid<-AMINO_ACID_CODE[aa_list[i]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],sum(Bacteria_meta$cor_OR_RSCU>0.5)))
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:701)])))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,log(OR)))+geom_violin(fill="gray",color="black")
  g<-g+ylim(-2,2)
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=13,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  p[[i]]<-g
}
png(filename="FigS1c1.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

p = list()
for(i in 1:9){
  odds_ratio_tmp<-Bacteria_OR_StrongSelection%>%
    filter(aa==aa_list[i+9])
  amino_acid<-AMINO_ACID_CODE[aa_list[i+9]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],sum(Bacteria_meta$cor_OR_RSCU>0.5)))
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:701)])))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,log(OR)))+geom_violin(fill="gray",color="black")
  g<-g+ylim(-2,2)
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=13,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  p[[i]]<-g
}
png(filename="FigS1c2.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

#Fig S2a
png(filename="./FigS2a.png",width=8,height=5,unit="in",res=600)
ggplot(Archaea_meta,aes(x=cor_oddratio_RSCU))+
  theme_bw()+
  geom_histogram(colour="black",fill="#738FC1",aes(y=..count../sum(..count..)),alpha=0.8)+
  xlab(label="Correlations between RSCU and Odds-ratios across species")+
  ylab(label="Frequency")+
  theme(axis.text.x=element_text(size=20,family="Helvetica",color="black"))+
  theme(axis.text.y=element_text(size=20,family="Helvetica",color="black"))+
  theme(axis.title.x=element_text(size=16,family="Lucida Sans",color="black",face="bold"))+
  theme(axis.title.y=element_text(size=20,family="Lucida Sans",color="black",face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#FFEDCE"), axis.line = element_line(colour = "black"))
dev.off()

#Fig S2b
species<-c(names(Archaea_OR)[3:65])
pruned.phy<-drop.tip(phy,phy$tip.label[-match(species,phy$tip.label)])
species_idx<-match(pruned.phy$tip.label,species)
traits<-data.frame(log(c(as.numeric(as.character(Archaea_OR[21,3:65])))[species_idx]))
names(traits)<-"colors_plot"
row.names(traits)<-pruned.phy$tip.label
png(filename="./FigS2b.png",width=5,height=5,unit="in",res=600)
p<-ggtree(pruned.phy,layout='circular')
gheatmap(p, traits, offset=0.2, width=0.2,font.size=2,colnames = FALSE,low="red",high="green",legend_title = "ln (OR)")
dev.off()

#Fig S2c
aa_list<-unique(as.character(Archaea_OR$aa))
p = list()
for(i in 1:9){
  odds_ratio_tmp<-Archaea_OR%>%
    filter(aa==aa_list[i])
  amino_acid<-AMINO_ACID_CODE[aa_list[i]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],63))
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:65)])))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,log(OR)))+geom_violin(fill="gray",color="black")
  g<-g+ylim(-2,2)
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=13,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  p[[i]]<-g
}
png(filename="FigS2c1.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

p = list()
for(i in 1:9){
  odds_ratio_tmp<-Archaea_OR%>%
    filter(aa==aa_list[i+9])
  amino_acid<-AMINO_ACID_CODE[aa_list[i+9]]
  codon_tmp<-as.character(odds_ratio_tmp$codon)
  codon_to_draw<-character(0)
  or_to_draw<-numeric(0)
  for(j in 1:length(codon_tmp)){
    codon_to_draw<-c(codon_to_draw,rep(codon_tmp[j],63))
    or_to_draw<-c(or_to_draw,as.numeric(as.character(odds_ratio_tmp[j,c(3:65)])))
  }
  or_for_drawTab<-data.frame(codon_to_draw,or_to_draw)
  names(or_for_drawTab)<-c("codons","OR")
  or_for_drawTab$codons<-factor(or_for_drawTab$codons,levels=unique(or_for_drawTab$codons))
  g<-ggplot(or_for_drawTab,aes(codons,log(OR)))+geom_violin(fill="gray",color="black")
  g<-g+ylim(-2,2)
  g<-g+ggtitle(amino_acid)+theme(plot.title = element_text(hjust=0.5,family="Lucida Sans",size=20,face="bold"))
  g<-g+theme(axis.title=element_blank())
  g<-g+theme(axis.text.x=element_text(family="Helvectica",size=13,face="bold"))
  g<-g+theme(axis.text.y=element_text(family="Helvectica",size=18,face="bold"))
  g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white",color="black",size=2), axis.line = element_line(colour = "black"))
  p[[i]]<-g
}
png(filename="FigS2c2.png",units="in",width=10,height=10,res=600)
do.call(grid.arrange,c(p,nrow=3))
dev.off()

#Fig S3
png(filename="./FigS3.png",width=5,height=5,unit="in",res=600)
odds_ratios_tab<-Eukaryotic_OR[,c(3:7)]
pairs(odds_ratios_tab,
      pch=20,
      col="#738FC1")
dev.off()

#Fig S4a
png(file="./FigS4a.png",units="in",width=6,height=6,res=600)
ggplot(Relative_error_tab,aes(x=relative_cog_tRNA,y=RMR))+
  geom_point(color="#738FC1",alpha=1,size=3)+
  geom_smooth(method="lm",se=FALSE,color="#D56C55")+
  geom_errorbar(aes(ymin=RMR-sd_of_RMR,ymax=RMR+sd_of_RMR),size=0.1,width=0.001,position = position_dodge(0.9))+
  xlab("Relative cognate tRNA concentration")+
  ylab("RMR")+
  theme(axis.text = element_text(family="Helvectica",size=15))+
  theme(axis.title=element_text(family="Lucida Sans",size=15,face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="#FFEDCE"),axis.line=element_line(colour="black"))

dev.off()

#Fig S4b
png(file="./FigS4b.png",units="in",width=6,height=6,res=600)
ggplot(Relative_error_tab,aes(x=relative_cog_tRNA,y=RSCU))+
  geom_point(color="#738FC1",alpha=1,size=3)+
  geom_smooth(method="lm",se=FALSE,color="#D56C55")+
  xlab("Relative cognate tRNA concentration")+
  ylab("RSCU")+
  theme(axis.text = element_text(family="Helvectica",size=15))+
  theme(axis.title=element_text(family="Lucida Sans",size=15,face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_rect(fill="#FFEDCE"),axis.line=element_line(colour="black"))
dev.off()
