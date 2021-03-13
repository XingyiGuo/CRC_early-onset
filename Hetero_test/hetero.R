#### Heterogeneity race  ######
out.m1<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Asian1.txt",head=TRUE,sep="\t")
out.m2<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Black1.txt",head=TRUE,sep="\t")
out.m3<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-White1.txt",head=TRUE,sep="\t")
sigList1<-c("PIK3CA","APC","FAT1")
sigList2<-c("TGFBR2","CREBBP")
sigList3<-c("LRP1B","TP53","TCF7L2","KDR","DOCK8","FLT4","SMAD2","SMAD3")
Htest1<-c(sigList1,sigList2,sigList3)
out.m1 <- out.m1[out.m1$Gene %in% Htest1,]
out.m2 <- out.m2[out.m2$Gene %in% Htest1,]
out.m3 <- out.m3[out.m3$Gene %in% Htest1,]
rownames(out.m1)<-out.m1$Gene
rownames(out.m2)<-out.m2$Gene
rownames(out.m3)<-out.m3$Gene
out.m1<-out.m1[Htest1,]
out.m2<-out.m2[Htest1,]
out.m3<-out.m3[Htest1,]

out.m<-cbind(out.m1,out.m2,out.m3)
out.m<-as.data.frame(out.m)


out.hetro<-matrix("na",13,51)
for (i in 1:dim(out.m)[1])
{
beta<- as.numeric(out.m[i,c(2,18,34)]) 
se <- as.numeric(out.m[i,c(3,19,35)])  

Y<-beta
V<-se^2
W<-1/V
WY<-W*Y
WY2<-W*Y^2
SumW<-sum(W)
SumWY<-sum(WY)
SumWY2<-sum(WY2)
MeanW<-mean(W)
VarW<-var(W)
Fixed_effect<-SumWY/SumW
Fixed_se<-1/sqrt(SumW)
Q<-SumWY2-SumWY^2/SumW
df<-length(Y)-1
P_hetero<-1-pchisq(Q,df)


   a1 <- as.numeric(out.m[i,c(11:12,27:28,43:44)])
        a1 <- c(a1[1],a1[2]-a1[1],a1[1+2],a1[2+2]-a1[1+2],a1[1+4],a1[2+4]-a1[1+4])
        p1<-chisq.test(t(matrix(a1,2,3)))$p.value
       

out.hetro[i,] <- c(as.character(out.m$Gene[i]),as.character(out.m[i,-1]),Q, P_hetero,p1)
}
write.table(out.hetro,"Genes_HeterobyRace.results",row.names=F,quote=F,sep="\t")



#### Heterogeneity test by gender ######
out.m1<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Male1.txt",head=TRUE,sep="\t")
out.m2<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Female1.txt",head=TRUE,sep="\t")
Htest1<-c(as.character(out.m1[out.m1$P < 0.05,]$Gene), as.character(out.m2[out.m2$P < 0.05,]$Gene))
out.m1 <- out.m1[out.m1$Gene %in% Htest1,]
out.m2 <- out.m2[out.m2$Gene %in% Htest1,]
rownames(out.m1)<-out.m1$Gene
rownames(out.m2)<-out.m2$Gene
out.m1<-out.m1[Htest1,]
out.m2<-out.m2[Htest1,]

out.m<-cbind(out.m1,out.m2)
out.m<-as.data.frame(out.m)


out.hetro<-matrix("na",12,35)
for (i in 1:dim(out.m)[1])
{
beta<- as.numeric(out.m[i,c(2,18)]) 
se <- as.numeric(out.m[i,c(3,19)])  

Y<-beta
V<-se^2
W<-1/V
WY<-W*Y
WY2<-W*Y^2
SumW<-sum(W)
SumWY<-sum(WY)
SumWY2<-sum(WY2)
MeanW<-mean(W)
VarW<-var(W)
Fixed_effect<-SumWY/SumW
Fixed_se<-1/sqrt(SumW)
Q<-SumWY2-SumWY^2/SumW
df<-length(Y)-1
P_hetero<-1-pchisq(Q,df)
   
	a1 <- as.numeric(out.m[i,c(11:12,27:28)])
        a1 <- c(a1[1],a1[2]-a1[1],a1[1+2],a1[2+2]-a1[1+2])
        p1<-chisq.test(t(matrix(a1,2,2)))$p.value
	
	out.hetro[i,] <- c(as.character(out.m$Gene[i]),as.character(out.m[i,-1]),Q, P_hetero,p1)
}
write.table(out.hetro,"Genes_HeterobySex.results",row.names=F,quote=F,sep="\t")


