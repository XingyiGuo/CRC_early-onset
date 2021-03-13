## Figure 1: barplot #### top genes ######
out.m<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-ALL.txt",head=TRUE,sep="\t")
out.m<-as.data.frame(out.m)
out.m <- out.m[out.m$P < 0.05,]
out.mbar<-out.m[,8:10]
rownames(out.mbar)<-out.m$Gene
#out.mbar<-out.mbar[,-(1:2)]
#out.mbar$V2 <- as.numeric(as.character(out.mbar$V2))
#out.mbar$V3 <- as.numeric(as.character(out.mbar$V3))
#out.mbar$V4 <- as.numeric(as.character(out.mbar$V4))
out.mbar<-out.mbar[rev(order(out.mbar$Freq)),]
#out.mbar<- out.mbar[rownames(out.mbar)%in% CD$Gene,]
library(reshape)
library(RColorBrewer)
coul <- brewer.pal(2, "Set2")[1:2] 
#out.mbar<-melt(out.m,measure.vars = c("V2","V3","V4"))
pdf("Figure1.barplot_early_old_topMF.pdf")
barplot(t(out.mbar[,1]),beside=TRUE,ylim=c(0,0.8),border=NA)
barplot(t(out.mbar[,2:3]),beside=TRUE,col=coul, ylim=c(0,0.8),border=NA)
dev.off()
pdf("Figure1.Sig.barplot_early_old_topMF.pdf")
#barplot(t(out.mbar[,1]),beside=TRUE,ylim=c(0,0.8),border=NA)
barplot(t(out.mbar[,2:3]),beside=TRUE,col=coul, ylim=c(0,0.8),border=NA)
dev.off()

#### mutation freq plot by race ####### 
sigList1<-c("PIK3CA","APC","FAT1")
sigList2<-c("TGFBR2","CREBBP")
sigList3<-c("LRP1B","TP53","TCF7L2","KDR","DOCK8","FLT4","SMAD2","SMAD3")
sigList<-sigList2
Htest1<-c(sigList1,sigList2,sigList3)
out.m1<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Asian1.txt",head=TRUE,sep="\t")
out.m2<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Black1.txt",head=TRUE,sep="\t")
out.m3<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-White1.txt",head=TRUE,sep="\t")

out.m1 <- out.m1[out.m1$Gene %in% sigList,]
out.m2 <- out.m2[out.m2$Gene %in% sigList,]
out.m3 <- out.m3[out.m3$Gene %in% sigList,]
rownames(out.m1)<-out.m1$Gene
rownames(out.m2)<-out.m2$Gene
rownames(out.m3)<-out.m3$Gene
out.m1<-out.m1[sigList,]
out.m2<-out.m2[sigList,]
out.m3<-out.m3[sigList,]


out.m<-rbind(out.m1,out.m2,out.m3)
out.m<-as.data.frame(out.m)
label.1<-rep(1:2,3)
label.2<-c(rep("a",2),rep("b",2),  rep("c",2))
label.1<-paste(label.1,label.2,sep="")
out.m<-cbind(out.m,label.1)
out.m<-out.m[order(out.m$label.1),]
out.mbar<-out.m[,8:10]
#rownames(out.mbar)<-out.m$Gene
#out.mbar<-out.mbar[,-(1:2)]
#out.mbar$V2 <- as.numeric(as.character(out.mbar$V2))
#out.mbar$V3 <- as.numeric(as.character(out.mbar$V3))
#out.mbar$V4 <- as.numeric(as.character(out.mbar$V4))
#out.mbar<-out.mbar[rev(order(out.mbar$Freq)),]
#out.mbar<- out.mbar[rownames(out.mbar)%in% CD$Gene,]
#library(reshape)
#library(RColorBrewer)
#coul <- brewer.pal(2, "Set2")[1:2] 
#out.mbar<-melt(out.m,measure.vars = c("V2","V3","V4"))
den1<-rep(c(10,10,20,20,30,30),2)
ang1<-rep(c(45,45,0,0,135,135),2)
pdf("Figure1.SigBlack.barplot_early_old_topMF.pdf")
#barplot(t(out.mbar[,1]),beside=TRUE,ylim=c(0,0.8),border=NA)
barplot(t(out.mbar[,2:3]),beside=TRUE,col=coul,density= den1, angle=ang1, ylim=c(0,0.8))
dev.off()

### foreast plot #####
# An example for foreast plot searched from the google. 
### "/Users/xingyi/Dropbox/project/GENIE/CRC/newanalysis03082021/R" ###
### figure: forest plot ####
library(forestplot)
library(metafor)
genes_df <- read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-ALL.txt", header=T, sep="\t")
genes_df <- genes_df[genes_df$P > 0 & genes_df$P < 0.05,c(1:7)]

labs <- genes_df$Gene
yi   <- genes_df$Beta
sei  <- genes_df$SE
res  <- rma(yi=yi, sei=sei, method="FE")
data <- structure(list(OR  = c(NA,genes_df$OR),
                       low = c(NA,genes_df$X95CI1),
                       high = c(NA,genes_df$X95CI2)),
                       .Names = c("OR", "low", "high"),
                 #      row.names = c(NA,-11L), 
                       class = "data.frame")

labels <- cbind(c("Gene_ID","Gene_1","Gene_2","Gene_3","Gene_4","Gene_5","Gene_6"),
                   c("HR","0.83","0.61","0.85","0.77","0.75","0.81"),
                   c("low","0.78","0.51","0.8","0.7","0.68","0.76"),
                   c("high","0.89","0.74","0.9","0.84","0.83","0.87"))

print("....Creating the plot....")
#jpeg(filename="Hazard_ratio_plot.jpg",units="cm",width=20,height=17, res=800)
pdf("MSS-all.pdf")
forestplot(labels,
           data,new_page = TRUE,
           boxsize = .25,
           zero = 0.707,
           ci.vertices = TRUE,
           ci.vertices.height = 0.25,
           xlog=TRUE,
           cex = 0.1,
           graph.pos = 2,
           lwd.zero = gpar(lty=1, alpha = 1),
           lineheight = "auto",
           title = " ",
           txt_gp = fpTxtGp(label=gpar(fontfamily="Calibri")),
           col = fpColors(box="blue",line="black",zero = "black"),
           xlab="Odd ratio")
dev.off()
