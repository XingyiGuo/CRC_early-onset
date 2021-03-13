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
	######################################################
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




write.table(rownames(out.mbar),"Figure1.genes.list",quote=F,row.names=F)
CD <- read.table("cancer_driver.list",head=T)
        #dim(x1[x1$PRIMARY_RACE %in% 'Black' & x1$ETHNICITY %in% 'Non-Spanish/non-Hispanic',])[1]
        #tt<-A.matrix[A.matrix$PRIMARY_RACE %in% 'Black' & A.matrix$ETHNICITY %in% 'Non-Spanish/non-Hispanic',]
#       dim(A.matrix[!duplicated(A.matrix$PATIENT_ID),])

#        write.table(gene,"CRC/gene_freq.txt",quote=F,sep="\t",row.names=F)
#        A.matrix <- A.matrix[A.matrix$PRIMARY_RACE %in% 'Black' & A.matrix$ETHNICITY %in% 'Non-Spanish/non-Hispanic',] # group 1
#        gene1<- as.data.frame(table(A.matrix$Hugo_Symbol))
#        gene1<-gene1[as.character(gene1$Var1) %in% as.character(gene$Var1),]
#        write.table(gene1,"CRC/gene_g3.2.freq.txt",quote=F,sep="\t",row.names=F)
#        A.matrix <- A.matrix[A.matrix$PRIMARY_RACE %in% 'White'  & A.matrix$ETHNICITY %in% 'Non-Spanish/non-Hispanic',] # group 2
#        A.matrix <- A.matrix[A.matrix$ETHNICITY %in% 'Spanish/Hispanic',] # group 3 
#        A.matrix <- A.matrix[A.matrix$PRIMARY_RACE %in% 'Asian',]
#       A.matrix <- rbind(A.matrix1,A.matrix2,A.matrix3,A.matrix4)
        
        # MAIN model - single gene analysis II:  #### 
        # 1: age ####
        # 2: ethich group
        # 3: histological gruoup
        # 4: sex
 #       A.matrix<- A.matrix[A.matrix$SEX %in% 'Female',] # group SEXf
 #       A.matrix<- A.matrix[A.matrix$SEX %in% 'Male',] # group SEXm
 #       A.matrix<- A.matrix[A.matrix$CANCER_TYPE_DETAILED %in% 'Colon Adenocarcinoma',] # group 4x
 #       A.matrix<- A.matrix[A.matrix$CANCER_TYPE_DETAILED %in% 'Rectal Adenocarcinoma',] # group 5x 
#Sigg<-read.table("Full_sig.genes",head=TRUE,stringsAsFactor=F)
### to select genes based on frequecny #### 
	A.g<-A.matrix[(A.matrix$PATIENT_ID %in% tmb1$PATIENT_ID),]
        gene<- as.data.frame(table(A.g$Hugo_Symbol))
        gene <- gene[gene$Freq > 100,]
	gene <- gene[!(as.character(gene$Var1) %in% 'ADD'),]
	names(gene)<-c("Gene","Freq")
	gene <- gene[!(gene$Gene %in% c("GNAS-AS1","INPPL1","MEF2BNB-MEF2B","ZFHX3")),]
	write.table(gene,"xxgene.list",quote=F,row.names=F)
#### end : using main fucntion to select #####

Sigg<-read.table("216gene.list",head=TRUE,stringsAsFactor=F)
Sigg1<-Sigg[((Sigg$Gene %in% CD$Gene) & Sigg$Freq > 0.015)  | Sigg$Freq >0.02,]
j <-0
jjj <-0
#GENIETdata<-matrix("na",1000,23)
#out.m<-matrix("na",120,7)
out.m<-matrix("na",87,14)
	
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
		g1<-'PMS2'
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

#### figure 2 plot #####
geneList<-c("TP53","APC","KRAS","PIK3CA","SMAD4","FBXW7","LRP1B","SOX9","COL7A1","BRAF","TCF7L2","ATM","KMT2D","ARID1A",'ADD','TGFBR2',"PIK3CA","APC","FAT1","KDR","SMAD2","FLT4")
panel<-as.character(PANEL[PANEL$Hugo_Symbol %in% geneList,]$SEQ_ASSAY_ID)
Asub <- A.matrix[(A.matrix$PATIENT_ID %in% tmb1[tmb1$tmb < 10^1.25,]$PATIENT_ID) & (A.matrix$SEQ_ASSAY_ID %in% panel),]                                                                                 
Asub.e <-Asub[Asub$Hugo_Symbol %in% geneList & Asub$AGE ==0 ,]
Asub.o <-Asub[Asub$Hugo_Symbol %in% geneList & Asub$AGE ==1,]
write.table(Asub.e,"Top0.05_Early_genes_mut.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 
write.table(Asub.o,"Top0.05_Old_genes_mut.txt",quote=F,sep="\t",row.names=F) # add p.adj with BH 



geneList <- read.table("GENIE_DiffGENE_EarlyvsTypical_byCancertype2020-.txt",sep="\t",head=TRUE)
geneList<- geneList[geneList$adj < 0.05,]
geneList1<-as.character(geneList$Gene) 


0.0787, 0.0566, 0.0643)
se <- c(0.0097, 0.0049, 0.0066)

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
Q
P_hetero


### 
xh1<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_112220-MSS-AsianF.txt",head=TRUE,sep="\t")
xh2<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_112220-MSS-BlackF.txt",head=TRUE)
xh3<-read.table("GENIE_DiffGENE_EarlyvsTypical_CRC_112220-MSS-WhiteF.txt",head=TRUE)
rownames(xh1)<- xh1$Gene
rownames(xh2)<- xh2$Gene
rownames(xh3)<- xh3$Gene
tmb1y<-tmb1[tmb1$tmb < 10^1.25,]	
s1<-as.data.frame(table(tmb1$PRIMARY_RACE))
tt<-0
mf.mat<-matrix('NA',19,2)
for(xg in as.character(xh1$Gene))
{
	#AF1 <- xh1[as.character(xg),9]*s1[s1$Var1 %in% 'Asian',2] 
	a1 <- as.numeric(c(xh1[xg,c(11:12)],xh2[xg,c(11:12)],xh3[xg,c(11:12)]))
	a1 <- c(a1[1],a1[2]-a1[1],a1[1+2],a1[2+2]-a1[1+2],a1[1+4],a1[2+4]-a1[1+4])
	p1<-chisq.test(t(matrix(a1,2,3)))$p.value
	tt<- tt+1
	mf.mat[tt,]<-c(xg,p1)

}
tt

#4'UTR
#5'UTR
#Frame_Shift_Del
#Frame_Shift_Ins
#In_Frame_Del
#In_Frame_Ins
#Intron
#Missense_Mutation
#Nonsense_Mutation
#Nonstop_Mutation
#RNA
#Silent
#Splice_Region
#Splice_Site
#Targeted_Region
#Translation_Start_Site
#Variant_Classification

##### section II. GENIE data analysis ######
#### TCGA data analysis ##########
#cli<- read.table("TCGA/coad_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli1<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/acc_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli2<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/blca_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli3<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/brca_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli4<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/cesc_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli5<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/chol_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli6<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/coadread_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
#cli<- read.table("TCGA/coad_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli1<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/acc_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli2<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/blca_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli3<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/brca_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli4<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/cesc_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli5<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/chol_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli6<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/coadread_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
#### TCGA data analysis ##########
#cli<- read.table("TCGA/coad_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli1<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/acc_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli2<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/blca_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli3<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/brca_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli4<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/cesc_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli5<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/chol_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli6<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/coadread_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli7<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/coad_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli8<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/esca_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli9<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/gbm_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli10<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli11<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/kich_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli12<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/kirc_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli13<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/kirp_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli14<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/laml_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli15<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/lgg_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli16<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/lihc_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli17<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/luad_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli18<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/lusc_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
#cli19<- read.table("TCGA/meso_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli20<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/ov_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli21<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/paad_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli22<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/pcpg_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli23<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/prad_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli24<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/read_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli25<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/sarc_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli26<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/skcm_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli27<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/stad_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
cli28<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/tgct_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)
#cli29<- read.table("TCGA/thym_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli30<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/ucec_tcga_pan_can_atlas_2018_clinical_data.csv",head=TRUE)
cli31<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/uvm_tcga_pan_can_atlas_2018_clinical_data.tsv",sep="\t",head=TRUE)

###Patient ID	Sample ID, Cancer Studies	Cancer Type	Cancer Type Detailed	Center of sequencing	Diagnosis Age Oncotree Code Neoplasm Disease Stage American Joint Committee on Cancer CodePerson GenderRace Category,Center of sequencing,American Joint Committee on Cancer Metastasis Stage Code
#li17a<-subset(cli17,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))

cli1a<-subset(cli1,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli2a<-subset(cli2,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli3a<-subset(cli3,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli4a<-subset(cli4,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli5a<-subset(cli5,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli6a<-subset(cli6,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli7a<-subset(cli7,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli8a<-subset(cli8,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli9a<-subset(cli9,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli10a<-subset(cli10,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli11a<-subset(cli11,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli12a<-subset(cli12,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli13a<-subset(cli13,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli14a<-subset(cli14,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli15a<-subset(cli15,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli16a<-subset(cli16,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli17a<-subset(cli17,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli18a<-subset(cli18,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli20a<-subset(cli20,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli21a<-subset(cli21,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli22a<-subset(cli22,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli23a<-subset(cli23,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli24a<-subset(cli24,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli25a<-subset(cli25,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli26a<-subset(cli26,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli27a<-subset(cli27,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli28a<-subset(cli28,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli30a<-subset(cli30,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))
cli31a<-subset(cli31,select=c("Patient.ID","Sample.ID","Cancer.Studies","Cancer.Type","Cancer.Type.Detailed","Center.of.sequencing","Diagnosis.Age","Oncotree.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Person.Gender","Race.Category","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code"))

cli.TCGA<-rbind(cli1a,cli2a,cli3a,cli4a,cli5a,cli6a,cli7a,cli8a,cli9a,cli10a,cli11a,cli12a,cli13a,cli14a,cli15a,cli16a,cli17a,cli18a,cli20a,cli21a,cli22a,cli23a,cli24a,cli25a,cli26a,cli27a,cli28a,cli30a,cli31a) 


##  grep -v UTR data_mut1.matrix | grep -v Intron | grep -v Silent | grep -v Targeted_Region | grep -v RNA  > data_mut2.matrix


##  grep -v UTR data_mut1.matrix | grep -v Intron | grep -v Silent | grep -v Targeted_Region | grep -v RNA  > data_mut2.matrix

##  grep -v UTR data_mut1.matrix | grep -v Intron | grep -v Silent | grep -v Targeted_Region | grep -v RNA  > data_mut2.matrix

mut.dataT<- read.table("/scratch/sccs/GENIE/Project-earlyonset/xTCGA/data_mutTCGA1.matrix",head=TRUE,sep="\t")  

mut.dataT$Tumor_Sample_Barcode <- substr(mut.dataT$Tumor_Sample_Barcode,1,12)

Ctype<- as.data.frame(table(cli.TCGA$Cancer.Type))
Ctype <- Ctype[Ctype$Freq>50,]
#Ctype <- Ctype[Ctype$Var1 %in% geneList$Cancertype,]
mut.dataT$Patient.ID <-mut.dataT$Tumor_Sample_Barcode
j <-0
#TCGAdata<-matrix("na",120,25) 
out.m<-matrix("na",2500,7)
for(add in c(2:10,12:14,17:18,20))
{
#	add<-18
	type <- Ctype$Var1[add]
        x1<-cli.TCGA[cli.TCGA$Cancer.Type %in% type,]
	x1<- x1[!is.na(x1$Diagnosis.Age),]
	x1<- x1[!duplicated(x1$Patient.ID),]
        sub.mut<- mut.dataT[mut.dataT$Patient.ID %in% as.character(x1$Patient.ID),]
        A.matrix <- merge(sub.mut,x1,by.x = "Patient.ID")
        A.matrix$AGE<-cut(A.matrix$Diagnosis.Age, breaks = c(0,49,100),labels=c(0,1))
        gene<- as.data.frame(table(A.matrix$Hugo_Symbol))
	gene <- gene[gene$Var1 %in% geneList1,]
        gene <- gene[gene$Freq  > 2,]
	print(type)
        for(g1 in gene$Var1)
        {
		#print(".")
                Asub <- A.matrix
                Asub$y<-ifelse((Asub$Hugo_Symbol %in% g1) &  (Asub$Variant_Classification %in% mutation),1,0)
                Asub.1 <- Asub[Asub$y == 1,]
                Asub.2 <- Asub[!(Asub$Patient.ID %in% Asub.1$Patient.ID),]
                Asub.1<-Asub.1[!duplicated(Asub.1$Patient.ID),]
                Asub.2<-Asub.2[!duplicated(Asub.2$Patient.ID),]

                Bsub<-rbind(Asub.1,Asub.2)
 #               xx<-sd(Bsub$Person.Gender)
                if(add  %in% c(2,3, 5, 6, 7, 9 ,10, 14 ,18))
                {	
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing + Person.Gender + Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code + Race.Category +American.Joint.Committee.on.Cancer.Metastasis.Stage.Code  ,data = Bsub,family="binomial"))
                }
                if(add  == 19 )
                {	
			# no gender
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing + Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code + Race.Category +American.Joint.Committee.on.Cancer.Metastasis.Stage.Code  ,data = Bsub,family="binomial"))
                }
                if(add  == 12  )
                {	
			# no race
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing + Person.Gender + Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code  +American.Joint.Committee.on.Cancer.Metastasis.Stage.Code  ,data = Bsub,family="binomial"))
                }
                if(add  %in% c(8, 15, 17))
                {	
			# no stage & no meta
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing + Person.Gender  + Race.Category  ,data = Bsub,family="binomial"))
                }
                if(add  == 13 )
                {	
			# no stage + no gender
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing  + Race.Category   ,data = Bsub,family="binomial"))
                }
		if(add  == 4 )
                {	
			# no stage + no gender
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing  + Race.Category +American.Joint.Committee.on.Cancer.Metastasis.Stage.Code  ,data = Bsub,family="binomial"))
                }
                if(add == 1   )
                {
			# no meta 
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing + Person.Gender + Race.Category +Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code ,data = Bsub,family="binomial"))
                        #Aout <- summary(glm(y ~AGE  + CENTER + PRIMARY_RACE + SAMPLE_TYPE_DETAILED,data = Bsub,family="binomial"))
                }
                if(add == 16)
                {
			# no geder + mo meta + no stage + race
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing  ,data = Bsub,family="binomial"))
                        #Aout <- summary(glm(y ~AGE  + CENTER + PRIMARY_RACE + SAMPLE_TYPE_DETAILED,data = Bsub,family="binomial"))
                }
                if(add == 20)
                {
			# no race+ mo meta + no stage 
                        Aout <- summary(glm(y ~relevel(AGE,ref='1')  + Center.of.sequencing + Race.Category ,data = Bsub,family="binomial"))
                        #Aout <- summary(glm(y ~AGE  + CENTER + PRIMARY_RACE + SAMPLE_TYPE_DETAILED,data = Bsub,family="binomial"))
                }
                j<- j+1
		sum.coef<- Aout$coef

		est<-exp(sum.coef[2,1])
		upper.ci<-exp(sum.coef[2,1]+1.96*sum.coef[2,2])
		lower.ci<-exp(sum.coef[2,1]-1.96*sum.coef[2,2])
		out.m[j,]<- c(paste(type),g1,Aout$coeff[2,c(1,4)],est,lower.ci,upper.ci)
                #out.m[j,]<- c(paste(type),g1,Aout$coeff[2,c(1,4)])
		#### for statistical number: T refers to typical ####
		xtitle<-paste(type,g1) #### need to define allout.sig first ####
		if(xtitle %in% geneList$ID)
		{
		Pmut.T<- Asub.1[Asub.1$AGE==1,]  
		Pmut.E<- Asub.1[Asub.1$AGE==0,]  
		
		PmutT.table<- as.data.frame(table(Pmut.T$Variant_Classification)[c(1:4,6:8,11:12)])
		PmutE.table<- as.data.frame(table(Pmut.E$Variant_Classification)[c(1:4,6:8,11:12)])
		
		Pnonmut.T<-dim(Asub.2[Asub.2$AGE==1,])[1]  
		Pnonmut.E<-dim(Asub.2[Asub.2$AGE==0,])[1]  
		jjj <- jjj +1
		GENIETdata[jjj,] <- c("TCGA",paste(type),g1,PmutE.table$Freq, Pnonmut.E,PmutT.table$Freq,Pnonmut.T)		
		write.table(Asub.1,paste("TCGA",paste(type,g1,sep="."),sep="-"),sep="\t",quote=F,row.names=FALSE)
		}
        }
}

GENIETdata1 <- GENIETdata[1:111,] # need to check the matrix 
out.m1<-out.m[1:574,]  # need to chck the matrix 
write.table(out.m1,"TCGA_DiffGENE_EarlyvsTypical_byCancertype2020.txt",quote=F,sep="\t")   
out.m1<- read.table("TCGA_DiffGENE_EarlyvsTypical_byCancertype2020.txt",head=TRUE,sep="\t")   

out.m1<-as.data.frame(out.m1)
names(out.m1) <- c("Cancertype","Gene","Beta","P","OR","95CI1","95CI2")
out.m1$P<- as.numeric(as.character(out.m1$P))
out.m1$ID<- paste(out.m1$Cancertype,out.m1$Gene)
#out.TCGA.sig<- out.m1[out.m1$P < 0.05,]

#geneList$pdj<-padj 
#names(geneList) <- c("Cancertype","Gene","Beta","P","padj")  
#geneList$ID<- paste(geneList$Cancertype,geneList$Gene)
allout.sig<-merge(geneList,out.m1,by.x = "ID",by.y = "ID",  all.x = TRUE )
write.table(allout.sig,"GENIE_TCGA_EarlyvsTypical_byCancertype2020.txt",quote=F,sep="\t")   
sig<-read.table("GENIE_TCGA_EarlyvsTypical_byCancertype2020.txt",sep="\t")   
sig<-sig[!is.na(sig$P.y) & (sig$P.y < 0.05),]
##### II.END  ######
##### III.plot figures  ######
# check Pancreatic Cancer	KRAS
GENIETdata1<-as.data.frame(GENIETdata1)
names(GENIETdata1)<-c("Study","Cancertype","Gene",paste(mutation),"Early_nonMut",paste(mutation),"Typical_nonMut")
write.table(GENIETdata1,"GENIE_TCGA_sigGene.sta",quote=F,sep="\t")   
GENIE_TCGA<-read.table("GENIE_TCGA_sigGene.sta",sep="\t",head=TRUE)  
GENIE_TCGA$ID <- paste(GENIE_TCGA$Cancertype,GENIE_TCGA$Gene)
### barplot #####
library(colourlovers)
palette3 <- clpalette('92095') 


#palette1 <- clpalette('113451')
#palette2 <- clpalette('92095')
#palette3 <- clpalette('629637')
#palette4 <- clpalette('694737')


G_p<- GENIE_TCGA[GENIE_TCGA$ID %in% sig$ID,]
G_p1<- subset(G_p,select=c(24,1:13))
G_p2<- subset(G_p,select=c(24,1:3,14:23))
G_p1$Mut1<- apply(G_p1[,5:13],1,sum) # total mutaionts
G_p1$ALL<- apply(G_p1[,c(14:15)],1,sum) # sample size
G_p1$P1<- G_p1$Missense_Mutation/G_p1$ALL
G_p1$P2<- apply(G_p1[,5:8],1,sum)/G_p1$ALL 
G_p1$P3<- apply(G_p1[,10:11],1,sum)/G_p1$ALL 
G_p1$P4<- G_p1[,12]/G_p1$ALL 
G_p1$P5<- G_p1[,13]/G_p1$ALL 

G_p2$Mut1<- apply(G_p2[,5:13],1,sum) # total mutaionts
G_p2$ALL<- apply(G_p2[,c(14:15)],1,sum) # sample size
G_p2$P1<- G_p2$Missense_Mutation.1/G_p2$ALL
G_p2$P2<- apply(G_p2[,5:8],1,sum)/G_p2$ALL 
G_p2$P3<- apply(G_p2[,10:11],1,sum)/G_p2$ALL 
G_p2$P4<- G_p2[,12]/G_p2$ALL 
G_p2$P5<- G_p2[,13]/G_p2$ALL 

G_p1a<- G_p1[1:26,]
G_p1b<- G_p1[27:52,]

G_p2a<- G_p2[1:26,]
G_p2b<- G_p2[27:52,]
for(i in 1:26)
{
	print("*")
	data <- t(rbind( G_p1a[i,c(17:21)], G_p2a[i,c(17:21)]))
	pdf(paste("GENIE",paste(G_p1a$ID[i],".pdf",sep="."),sep="."))
	barplot(data , beside=F ,ylim=c(0,1.2*max(sum(G_p1a[i,17:21]),sum(G_p2a[i,17:21]))), col = swatch(palette3)[[1]],  border="NA")
	dev.off()
	data <- t(rbind( G_p1b[i,c(17:21)], G_p2b[i,c(17:21)]))
##### main Function ######
##### section I. GENIE data analysis ######
#3'UTR
#5'UTR
#Frame_Shift_Del
#Frame_Shift_Ins
#In_Frame_Del
#In_Frame_Ins
#Intron
cliP<- read.table("/scratch/sccs/GENIE/Project-earlyonset/data_clinical_patient.txt",head=TRUE,sep="\t")
cliS<- read.csv("/scratch/sccs/GENIE/Project-earlyonset/data_clinical_sample.csv",head=TRUE)

cli<-merge(cliP,cliS,by.x = "PATIENT_ID")
cli$AGE_AT_SEQ_REPORT <- gsub(">","",cli$AGE_AT_SEQ_REPORT,perl =TRUE)
cli$AGE <- gsub("<","",cli$AGE_AT_SEQ_REPORT,perl =TRUE)
cli$AGE <- as.numeric(as.character(cli$AGE))

#patients<- read.table("/scratch/sccs/GENIE/Project-earlyonset/phs001337.v1.pht009621.v1.p1.AACR_GENIE_Sample.MULTI.txt",head=TRUE,sep="\t")
#cli <- cli[cli$SAMPLE_ID %in% patients$SAMPLE_ID,]
#Ctype<- as.data.frame(table(cli$ONCOTREE_CODE)) 
Ctype<- as.data.frame(table(cli$CANCER_TYPE)) 
Ctype <- Ctype[Ctype$Freq > 200,]
mut.dat$SAMPLE_ID <-mut.dat$Tumor_Sample_Barcode
j <-0
jjj <-0
GENIETdata<-matrix("na",5000,23)
out.m<-matrix("na",5000,7)
for(type in Ctype$Var1)
{
	x1<-cli[cli$CANCER_TYPE %in% type,]
	x1<-x1[!duplicated(x1$SAMPLE_ID),]
	x1<-x1[!duplicated(x1$PATIENT_ID),]
