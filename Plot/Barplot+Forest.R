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

"/Users/xingyi/Dropbox/project/GENIE/CRC/newanalysis03082021/R"
library(metafor)
 
### copy BCG vaccine meta-analysis data into 'dat'
dat <- dat.bcg
 
### calculate log risk ratios and corresponding sampling variances (and use
### the 'slab' argument to store study labels as part of the data frame)
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat,
              slab=paste(author, year, sep=", "))
 
### fit random-effects model
res <- rma(yi, vi, data=dat)
 
### forest plot with extra annotations
forest(res, atransf=exp, at=log(c(.05, .25, 1, 4)), xlim=c(-16,6),
       ilab=cbind(dat.bcg$tpos, dat.bcg$tneg, dat.bcg$cpos, dat.bcg$cneg),
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75, header="Author(s) and Year",
       mlab="")
op <- par(cex=.75, font=2)
text(c(-9.5,-8,-6,-4.5), 15, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75,-5.25),     16, c("Vaccinated", "Control"))
par(op)
 
### add text with Q-value, dfs, p-value, and I^2 statistic
text(-16, -1, pos=4, cex=0.75, bquote(paste("RE Model (Q = ",
     .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
     ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
     .(formatC(res$I2, digits=1, format="f")), "%)")))

