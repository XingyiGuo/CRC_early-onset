##### main Function ######
##### section I. GENIE data analysis ######
#3'UTR
#5'UTR
#Frame_Shift_Del
#Frame_Shift_Ins
#In_Frame_Del
#In_Frame_Ins
#Intron
#Missense_Mutation
#Nonsense_Mutation
#Nonstop_Mutation
#Nonstop_Mutation
#RNA
#Silent
#Splice_Region
#Splice_Site
#Targeted_Region
#Translation_Start_Site
#Variant_Classification
library(gdata)
library(ggplot2)
require(scales)

library(data.table)
mutation<-c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site","Translation_Start_Site")

mut.dat<-fread(".../path/...Project-earlyonset/Endometrial/GENIE_9/data_mutations_extended.txt",head=TRUE,sep = "\t")

mut.dat$SAMPLE_ID <-mut.dat$Tumor_Sample_Barcode
cliP<- read.table(".../path/...Project-earlyonset/Endometrial/GENIE_9/data_clinical_patient.txt",head=TRUE,sep="\t")
cliS<- fread(".../path/...Project-earlyonset/Endometrial/GENIE_9/data_clinical_sample.txt",sep="\t",head=TRUE)
cliS<-as.data.frame(cliS)

cli<-merge(cliP,cliS,by = "PATIENT_ID",all = FALSE)
cli<-cli[!duplicated(cli$PATIENT_ID),]

cli$AGE_AT_SEQ_REPORT <- gsub(">89","Unknown",cli$AGE_AT_SEQ_REPORT,perl =TRUE)
cli$AGE_AT_SEQ_REPORT <- gsub("<18","Unknown",cli$AGE_AT_SEQ_REPORT,perl =TRUE)

cli$AGE <- cli$AGE_AT_SEQ_REPORT
#AGEM<-cli[!(cli$AGE %in% 'Unknown'),]$AGE
#AGEM <- as.numeric(as.character(AGEM))   
#AGEM<-median(AGEM)
#cli[cli$AGE %in% 'Unknown',]$AGE <- AGEM
cli<-cli[!(cli$AGE %in% 'Unknown'),]
cli$AGE <- as.numeric(as.character(cli$AGE))


PANEL<-read.table(".../path/...Project-earlyonset/Endometrial/GENIE_9/genomic_information_noFalse.txt",head=TRUE,sep="\t")
TMB_ref<-read.table(".../path/...Project-earlyonset/CRC/GENIE_PANEL.COV",head=TRUE)

PANEL$IDu<-paste(PANEL$Hugo_Symbol,PANEL$SEQ_ASSAY_ID,sep="")
PANEL<-PANEL[!duplicated(PANEL$IDu),]

Ctype<- as.data.frame(table(cli$CANCER_TYPE)) 
Ctype <- Ctype[Ctype$Freq > 200,]

        type <- 'Colorectal Cancer' 
        x1<-cli[cli$CANCER_TYPE %in% type,]
        x1<-x1[!duplicated(x1$SAMPLE_ID),]
        x1<-x1[!duplicated(x1$PATIENT_ID),]
	QC1<- c('Colon Adenocarcinoma In Situ','Medullary Carcinoma of the Colon','Signet Ring Cell Adenocarcinoma of the Colon and Rectum') 
	x1 <- x1[!(x1$CANCER_TYPE_DETAILED %in% QC1),]
	
	
	#histology two groups:
	hg1<-c("Colon Adenocarcinoma","Colorectal Adenocarcinoma","Rectal Adenocarcinoma")                                
	hg2<-"Mucinous Adenocarcinoma of the Colon and Rectum"
                                
  sub.mut<- mut.dat[mut.dat$SAMPLE_ID %in% x1$SAMPLE_ID,]
	sub.mut <- as.data.frame(sub.mut)
       # sub.mut<-sub.mut[sub.mut$Variant_Classification %in% mutation,]
        A.matrix.1 <- merge(sub.mut,x1,by.x = "SAMPLE_ID")
	
	### add samples without mutations dected: start  ####
	x2 <- x1[!(x1$SAMPLE_ID %in% sub.mut$SAMPLE_ID),] ## to add samples without mutations 
	sub.mut.1 <- sub.mut[1:dim(x2)[1],] ### total of samples: 385 -338 
	sub.mut.1$SAMPLE_ID  <- x2$SAMPLE_ID  	
	sub.mut.1$Hugo_Symbol <-  "ADD"
	A.matrix.2 <- merge(sub.mut.1,x2,by.x = "SAMPLE_ID")    
	### end #######
	#### merge A.matrix + A.matrix.1 ####
	A.matrix <- rbind(A.matrix.2,A.matrix.1)
       
	A.matrix$AGE<-cut(as.numeric(A.matrix$AGE), breaks = c(0,49,100),labels=c(0,1))
#	A.matrix$AGE<-cut(A.matrix$AGE, breaks = c(0,49,59,69,100),labels=c(0,1,2,3))
	######################################################
  ######################################################
	######################################################
	### main analysis I: ##### CRC / baselines ##########
	A.base <- A.matrix
	A.base$y<-ifelse((A.base$Variant_Classification %in% mutation),1,0) 		
	A.base1 <- merge(A.base,TMB_ref,by.x = 'SEQ_ASSAY_ID')
	assay<- A.base1[!duplicated(A.base1$PATIENT_ID),]
	tmb<- as.data.frame(table(as.character(A.base1$PATIENT_ID)))
	names(tmb)<- c("PATIENT_ID","Freq")
	tmb<-merge(tmb,assay,by.x = 'PATIENT_ID')
	tmb$tmb<- tmb$Freq/tmb$ts
	### test assays to rem those with poor QCs ####
	
	tmb<-tmb[tmb$ts > 50000 & !(tmb$SEQ_ASSAY_ID %in% c('YALE-OCP-V2','UCHI-ONCOHEME55-V1','NKI-PATH-NGS')),]
	pdf("CRC_tmb_assay.pdf")
	boxplot(tmb$tmb ~ tmb$SEQ_ASSAY_ID)
	dev.off()
	### main analysis I: end #####
	######################################################
  ######################################################
	######################################################
  #### only focusing on Asian/Pacific Islander + Non-Spanish/non-Hispanic Black + Non-Spanish/non-Hispanic White
  tmb1 <- tmb
	table(tmb1$PRIMARY_RACE,tmb1$ETHNICITY)
	tmb1$PRIMARY_RACE <- gsub("Pacific Islander","Asian",tmb1$PRIMARY_RACE,perl=TRUE)
	#tmb1$PRIMARY_RACE  <- ifelse(tmb1$PRIMARY_RACE %in% "Asian","Asian",ifelse(tmb1$PRIMARY_RACE %in% "Black","Black",ifelse(tmb1$PRIMARY_RACE %in% "White","White","Others")))
	tmb1<-tmb1[tmb1$PRIMARY_RACE %in% 'Asian' | (tmb1$PRIMARY_RACE %in% c("Black","White") & tmb1$ETHNICITY %in% 'Non-Spanish/non-Hispanic'),]	
	
  write.table(tmb1,"GENIE_CRC_Table1_populations",row.names=F,sep="\t",quote=F)

	######################################################
  ######################################################
	######################################################
	### main analysis II: ##### overall analysis in TMB ##########
	tmb1$lmtmb<-resid(lm(tmb ~SEQ_ASSAY_ID, data=tmb1))
	tmb1<-tmb1[!(tmb1$PRIMARY_RACE %in% 'Others'),]
	tmb1$tmb <- tmb1$tmb * 1e6
#	tmb1$CANCER_TYPE_DETAILED <- ifelse(!(tmb1$CANCER_TYPE_DETAILE %in% 'Mucinous Adenocarcinoma of the Colon and Rectum'),'Colon','non-Colon')
	st1<- summary(lm(tmb ~ relevel(AGE, ref = "1") +relevel(as.factor(PRIMARY_RACE), ref = "White") +SEX+  CANCER_TYPE_DETAILED + SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID,data=tmb1[tmb1$tmb > 10^1.25,]))
	st1<- summary(lm(tmb ~ relevel(AGE, ref = "1") +relevel(as.factor(PRIMARY_RACE), ref = "White") +SEX+  CANCER_TYPE_DETAILED + SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID,data=tmb1[tmb1$tmb <= 10^1.25,]))
	# by age for non-hyper group
	st1<- summary(lm(tmb ~ relevel(as.factor(PRIMARY_RACE), ref = "White") +SEX+  CANCER_TYPE_DETAILED + SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID,data=tmb1[tmb1$tmb <= 10^1.25 & tmb1$AGE==0,]))
	st1<- summary(lm(tmb ~ relevel(as.factor(PRIMARY_RACE), ref = "White") +SEX+  CANCER_TYPE_DETAILED + SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID,data=tmb1[tmb1$tmb <= 10^1.25 & tmb1$AGE==1,]))
	#st1<- summary(lm(lmtmb ~ relevel(AGE, ref = "1") + relevel(as.factor(PRIMARY_RACE), ref = "White") +  CANCER_TYPE_DETAILED + SAMPLE_TYPE_DETAILED,data=tmb1))
	#st1<- summary(lm(lmtmb ~ relevel(AGE, ref = "1") ,data=tmb1))
  
  	### main analysis II: end #####
	######################################################
  ######################################################
	######################################################
	
  ######################################################
  ######################################################
	######################################################	
  #### Table 1 ######
tmb1$AGE_AT_SEQ_REPORT <- as.numeric(tmb1$AGE_AT_SEQ_REPORT)
dim(tmb1[tmb1$AGE_AT_SEQ_REPORT < 30,])[1]
dim(tmb1[tmb1$AGE_AT_SEQ_REPORT >= 30 & tmb1$AGE_AT_SEQ_REPORT <40,])[1]
dim(tmb1[tmb1$AGE_AT_SEQ_REPORT >= 40 & tmb1$AGE_AT_SEQ_REPORT <50,])[1]
dim(tmb1[tmb1$AGE_AT_SEQ_REPORT >= 50 & tmb1$AGE_AT_SEQ_REPORT <60,])[1]
dim(tmb1[tmb1$AGE_AT_SEQ_REPORT >= 60 & tmb1$AGE_AT_SEQ_REPORT <70,])[1]
dim(tmb1[tmb1$AGE_AT_SEQ_REPORT >= 70 & tmb1$AGE_AT_SEQ_REPORT <80,])[1]
dim(tmb1[tmb1$AGE_AT_SEQ_REPORT >= 80 & tmb1$AGE_AT_SEQ_REPORT <100,])[1]

E1 <- tmb1[tmb1$AGE == 0,]
E1 <- tmb1[tmb1$AGE == 1,]

E1.1<-(E1[E1$AGE_AT_SEQ_REPORT < 30,])
E1.2<-(E1[E1$AGE_AT_SEQ_REPORT >= 30 & E1$AGE_AT_SEQ_REPORT <40,])
E1.3<-(E1[E1$AGE_AT_SEQ_REPORT >= 40 & E1$AGE_AT_SEQ_REPORT <50,])
E1.4<-(E1[E1$AGE_AT_SEQ_REPORT >= 50 & E1$AGE_AT_SEQ_REPORT <60,])
E1.5<-(E1[E1$AGE_AT_SEQ_REPORT >= 60 & E1$AGE_AT_SEQ_REPORT <70,])
E1.6<-(E1[E1$AGE_AT_SEQ_REPORT >= 70 & E1$AGE_AT_SEQ_REPORT <80,])
E1.7<-(E1[E1$AGE_AT_SEQ_REPORT >= 80 & E1$AGE_AT_SEQ_REPORT <100,])

MSI.1<-tmb1[tmb1$MSI %in% 'MSS' & tmb1$AGE == 0,]

table(MSI.1$CANCER_TYPE_DETAILED)

table(E1$PRIMARY_RACE, E1$SEX)

  ####  Table 1: END ###### 
  #####################################################
  ######################################################
  ######################################################

	
	#### end the analysis checking by feb17,2021 #### 
	#### end the tmb analysis I #### 
	tmb1a<-tmb1[tmb1$tmb <= 10^1.25,]
	pdf("Non-hyper_tmb_groups.pdf")
	boxplot(tmb1a$lmtmb~tmb1a$AGE,ylim=c(-1.369e-05,0.00001))
	boxplot(lmtmb~PRIMARY_RACE,data=tmb1a[tmb1a$AGE==0,],ylim=c(-1.369e-05,0.00001))
	boxplot(lmtmb~PRIMARY_RACE,data=tmb1a[tmb1a$AGE==1,],ylim=c(-1.369e-05,0.00001))
	boxplot(lmtmb~SEX,data=tmb1[tmb1a$AGE==0,],ylim=c(-1.369e-05,0.00001))
	boxplot(lmtmb~SEX,data=tmb1[tmb1a$AGE==1,],ylim=c(-1.369e-05,0.00001))
	dev.off()
	
	#### plot baseline Early-onset vs Old-onset ####
	x1<-seq(1:length(rev(sort(tmb1$tmb))))
	y1<-rev(sort(tmb1$tmb))
	fig1a<-as.data.frame(cbind(x1,y1))
	
	jpeg(file="Figure1_GENIE_tmb.jpeg",width = 480, height = 480,
     	pointsize = 12, quality = 100, bg = "white", res = NA)
	p2 <- ggplot(fig1a, aes(x = x1, y = y1)) + geom_point() +
     	scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
    	 theme_bw() 
	p2	
	dev.off()




######################################################
 ######################################################
######################################################		      
#### main analysis II: evaluating gene mutations assoicated with early-onset/old-onset.
### to select genes based on frequecny #### 
	A.g<-A.matrix[(A.matrix$PATIENT_ID %in% tmb1$PATIENT_ID),]
        gene<- as.data.frame(table(A.g$Hugo_Symbol))
        gene <- gene[gene$Freq > 100,] ### check all posible frequently mutated genes
	gene <- gene[!(as.character(gene$Var1) %in% 'ADD'),]
	names(gene)<-c("Gene","Freq")
	gene <- gene[!(gene$Gene %in% c("GNAS-AS1","INPPL1","MEF2BNB-MEF2B","ZFHX3")),]
	write.table(gene,"216gene.list",quote=F,row.names=F) ### need to modifiy based on the gene assocation testing (i.e. removing genes with low numbers of mutations but high frquency###
#### end : using main fucntion to select #####
CD <- read.table("cancer_driver.list",head=T) ## add addtional cancer discovery genes from TCGA using a mut frq at 1.5%.
Sigg<-read.table("216gene.list",head=TRUE,stringsAsFactor=F)
Sigg1<-Sigg[((Sigg$Gene %in% CD$Gene) & Sigg$Freq > 0.015)  | Sigg$Freq >0.02,]
j <-0
jjj <-0
#GENIETdata<-matrix("na",1000,23)
#out.m<-matrix("na",120,7)
out.m<-matrix("na",87,14) ### total 87 investigated genes ### 
	
	tmb1$MSI<- ifelse(tmb1$tmb < 10^1.25,"MSS","MSI")
	tmb1a<-tmb1
	rownames(tmb1a)<-tmb1$PATIENT_ID
	for(g1 in Sigg1$Gene)
	#for(g1 in gene$Gene)
	{
	#	g1 <- 'TP53'
	#	g1 <- 'APC'
	#	g1 <- 'RNF43'
	#	g1 <- 'MSH2'
	#	g1 <- 'SMAD2'
	#	g1<-'PMS2'
		panel<-as.character(PANEL[PANEL$Hugo_Symbol %in% g1,]$SEQ_ASSAY_ID)
		Asub <- A.matrix[(A.matrix$PATIENT_ID %in% tmb1[tmb1$tmb > 10^1.25,]$PATIENT_ID) & (A.matrix$SEQ_ASSAY_ID %in% panel),]
#		Asub$CANCER_TYPE_DETAILED <- ifelse(!(Asub$CANCER_TYPE_DETAILE %in% 'Mucinous Adenocarcinoma of the Colon and Rectum'),'Colon','non-Colon')
#		Asub <- A.matrix[(A.matrix$PATIENT_ID %in% tmb1[tmb1$tmb < 10^1.25,]$PATIENT_ID),]
#		Asub.w <- Asub[Asub$PRIMARY_RACE %in% 'White',]
#		Asub.b <- Asub[Asub$PRIMARY_RACE %in% 'Black',]
#		Asub.a <- Asub[Asub$PRIMARY_RACE %in% 'Asian',]
#		write.table(Asub.w,"Non-hyper-white.mut.matrix",quote=F,sep="\t")
#		write.table(Asub.b,"Non-hyper-black.mut.matrix",quote=F,sep="\t")
#		write.table(Asub.a,"Non-hyper-asian.mut.matrix",quote=F,sep="\t")
	#	Asub <- A.matrix[(A.matrix$PATIENT_ID %in% tmb1$PATIENT_ID) & (A.matrix$SEQ_ASSAY_ID %in% panel),]
	#	Asub<-Asub[Asub$PRIMARY_RACE %in% 'White',]
	#	Asub<-Asub[Asub$PRIMARY_RACE %in% 'Asian',]
	#	Asub<-Asub[Asub$PRIMARY_RACE %in% 'Black',]
	#	Asub<-Asub[Asub$SEX %in% 'Female',]
	#	Asub<-Asub[(Asub$CANCER_TYPE_DETAILED %in% c('Colon Adenocarcinoma',"Colorectal Adenocarcinoma", "Rectal Adenocarcinoma")),]
	#	Asub<-Asub[(Asub$CANCER_TYPE_DETAILED %in% "Mucinous Adenocarcinoma of the Colon and Rectum"),]
		### remove panel genes due to insufficent cover #####
		###### 
		Asub$y<-ifelse((Asub$Hugo_Symbol %in% g1) & (Asub$Variant_Classification %in% mutation),1,0)
		Asub.1 <- Asub[Asub$y == 1,]
		Asub.2 <- Asub[!(Asub$SAMPLE_ID %in% Asub.1$SAMPLE_ID),]
		Asub.1<-Asub.1[!duplicated(Asub.1$SAMPLE_ID),]
		Asub.2<-Asub.2[!duplicated(Asub.2$SAMPLE_ID),]
		Bsub<-rbind(Asub.1,Asub.2)
		
		frq<- dim(Asub.1)[1]/dim(Bsub)[1]
		frq1<- dim(Asub.1[Asub.1$AGE == 0,])[1]/dim(Bsub[Bsub$AGE==0,])[1]
		frq2<- dim(Asub.1[Asub.1$AGE ==1,])[1]/dim(Bsub[Bsub$AGE==1,])[1]
		My1<-dim(Asub.1[Asub.1$AGE == 0,])[1]
		My0<-dim(Bsub[Bsub$AGE == 0,])[1]
		Mo1<-dim(Asub.1[Asub.1$AGE == 1,])[1]
		Mo0<-dim(Bsub[Bsub$AGE == 1,])[1]
	#	j <- j+1
	#	out.m[j,]<- c(g1,frq,frq1,frq2,My1,My0,Mo1,Mo0)
	#	if(frq > 0.02)
	#	{
	#	table(Bsub$AGE,Bsub$y)    
		sample.list<- as.character(Bsub$PATIENT_ID)
		MSI <- tmb1a[sample.list,]$tmb
		
	#	Aout <- summary(glm(y ~relevel(AGE,ref='1')  + SEX + relevel(PRIMARY_RACE,ref='White') + as.factor(CANCER_TYPE_DETAILED)+ SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID + MSI ,data = Bsub,family="binomial"))
		Aout <- summary(glm(y ~relevel(AGE,ref='1')  + SEX + as.factor(CANCER_TYPE_DETAILED)+ SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID + MSI ,data = Bsub,family="binomial"))
	#	Aout <- summary(glm(y ~relevel(AGE,ref='1')  + SEX+relevel(PRIMARY_RACE,ref='White') +  SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID + MSI ,data = Bsub,family="binomial"))
	#	Aout <- summary(glm(y ~relevel(AGE,ref='1') + relevel(PRIMARY_RACE,ref='White') +as.factor(CANCER_TYPE_DETAILED) + SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID + MSI ,data = Bsub,family="binomial"))
		#Aout <- summary(glm(y ~relevel(AGE,ref='1') + relevel(PRIMARY_RACE,ref='White')+as.factor(CANCER_TYPE_DETAILED) + SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID + MSI ,data = Bsub,family="binomial"))
		#Aout <- summary(glm(y ~relevel(AGE,ref='1') +SEX+ relevel(PRIMARY_RACE,ref='White') + SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID + MSI ,data = Bsub,family="binomial"))
		#Aout <- summary(glm(y ~relevel(AGE,ref='1')  + relevel(as.factor(CANCER_TYPE_DETAILED),ref='Colon') + SAMPLE_TYPE_DETAILED + SEQ_ASSAY_ID + MSI ,data = Bsub,family="binomial"))
		
		j<- j+1

		sum.coef<- Aout$coef

		est<-exp(sum.coef[2,1])
		upper.ci<-exp(sum.coef[2,1]+1.96*sum.coef[2,2])
		lower.ci<-exp(sum.coef[2,1]-1.96*sum.coef[2,2])
		out.m[j,]<- c(g1,Aout$coeff[2,c(1,2,4)],est,lower.ci,upper.ci,frq,frq1,frq2,My1,My0,Mo1,Mo0)
		}
#	}

j
out.m1<-out.m[1:87,]
out.m1<-as.data.frame(out.m1)
names(out.m1) <- c("Gene","Beta","SE","P","OR","95CI1","95CI2","Freq","Freq_E","Freq_O","#_Carrie_E","#_EAll","#_Carrie_O","#_OAll")
out.m1$P<- as.numeric(as.character(out.m1$P))
#out.m1$Freq<- as.numeric(as.character(out.m1$Freq))

out.m1$adj<-p.adjust(out.m1$P, method="BH")
out.m1$ID<- paste(out.m1$Cancertype,out.m1$Gene)

#write.table(out.m1,"GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-ALLa.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
#write.table(out.m1,"GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-White.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
#write.table(out.m1,"GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Asian.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
write.table(out.m1,"GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Black.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
#write.table(out.m1,"GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Male.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
#write.table(out.m1,"GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Female.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
write.table(out.m1,"GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Colon+Rectal_Adenocarcinoma.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
#write.table(out.m1,"GENIE_DiffGENE_EarlyvsTypical_CRC_02192021-MSS-Mucinous.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
#write.table(tmb1,"GENIE_CRC_Table1_populations",quote=F,sep="\t",row.names=F) # add p.adj with BH 

