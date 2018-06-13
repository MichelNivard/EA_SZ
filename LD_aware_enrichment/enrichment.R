library(data.table)
library(gplots)  


########################**Example (BIP) for calculating slope and intercept for different traits**########################

# python ~/ldsc/munge_sumstats.py --merge-alleles ~/Downloads/eur_w_ld_chr/w_hm3.snplist --sumstats pgc.bip.full.2012-04.txt --a1 a1 --a2 a2 --signed-sumstats or,1 --N 16731 --out BIP
# python ~/ldsc/ldsc.py --h2 BIP.sumstats.gz --ref-ld-chr ~/Downloads/eur_w_ld_chr/ --w-ld-chr ~/Downloads/eur_w_ld_chr/

#########################Pre-calculated values for all used traits are provided in "Slopes_etc_values.txt"#################

####Pre-calcultaed LD scorse for Hapmap3 snplist were downloaded along with the LDSC software https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
####LD score of 105 SNPs used in this study are provided in "LD_Score_of_Best_Bud_105.txt" file




###Calculate expected chi-square stats

SNP_with_score=read.delim("LD_Score_of_Best_Bud_105.txt", stringsAsFactors=F)
Slopes=read.delim("Slopes_etc_values.txt", stringsAsFactors=F)


SNP_with_score$Height=(Slopes[1,2]*Slopes[1,4]*SNP_with_score$L2/1173569)+Slopes[1,3]
SNP_with_score$SZ=(Slopes[2,2]*Slopes[2,4]*SNP_with_score$L2/1173569)+Slopes[2,3]
SNP_with_score$BIP=(Slopes[3,2]*Slopes[3,4]*SNP_with_score$L2/1173569)+Slopes[3,3]
SNP_with_score$IBD=(Slopes[4,2]*Slopes[4,4]*SNP_with_score$L2/1173569)+Slopes[4,3]
SNP_with_score$ADHD_C=(Slopes[5,2]*Slopes[5,4]*SNP_with_score$L2/1173569)+Slopes[5,3]
SNP_with_score$Alzheimer=(Slopes[6,2]*Slopes[6,4]*SNP_with_score$L2/1173569)+Slopes[6,3]
SNP_with_score$Autism_C=(Slopes[7,2]*Slopes[7,4]*SNP_with_score$L2/1173569)+Slopes[7,3]
SNP_with_score$CAD=(Slopes[8,2]*Slopes[8,4]*SNP_with_score$L2/1173569)+Slopes[8,3]
SNP_with_score$Menopause=(Slopes[9,2]*Slopes[9,4]*SNP_with_score$L2/1173569)+Slopes[9,3]
SNP_with_score$Smoke_CPD=(Slopes[10,2]*Slopes[10,4]*SNP_with_score$L2/1173569)+Slopes[10,3]
SNP_with_score$MDD=(Slopes[11,2]*Slopes[11,4]*SNP_with_score$L2/1173569)+Slopes[11,3]
SNP_with_score$Neuroticism=(Slopes[12,2]*Slopes[12,4]*SNP_with_score$L2/1173569)+Slopes[12,3]
SNP_with_score$Fasting_Ins=(Slopes[13,2]*Slopes[13,4]*SNP_with_score$L2/1173569)+Slopes[13,3]
SNP_with_score$BMI=(Slopes[14,2]*Slopes[14,4]*SNP_with_score$L2/1173569)+Slopes[14,3]
SNP_with_score$SWB=(Slopes[15,2]*Slopes[15,4]*SNP_with_score$L2/1173569)+Slopes[15,3]
SNP_with_score$DS=(Slopes[16,2]*Slopes[16,4]*SNP_with_score$L2/1173569)+Slopes[16,3]
SNP_with_score$Menarche=(Slopes[17,2]*Slopes[17,4]*SNP_with_score$L2/1173569)+Slopes[17,3]
SNP_with_score$ICV=(Slopes[18,2]*Slopes[18,4]*SNP_with_score$L2/1173569)+Slopes[18,3]
SNP_with_score$CHIC=(Slopes[19,2]*Slopes[19,4]*SNP_with_score$L2/1173569)+Slopes[19,3]
SNP_with_score$F_PIns=(Slopes[20,2]*Slopes[20,4]*SNP_with_score$L2/1173569)+Slopes[20,3]
SNP_with_score$BW=(Slopes[21,2]*Slopes[21,4]*SNP_with_score$L2/1173569)+Slopes[21,3]
SNP_with_score$BL=(Slopes[22,2]*Slopes[22,4]*SNP_with_score$L2/1173569)+Slopes[22,3]


#******************************************************************************************************************************************################
###Calculate observed chi-square stats
#Z value from Pvalue 

Height=read.delim("../Slope_Intercept/GWAS_Downloads/Height/HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt", stringsAsFactors=F)
Height$Z=qnorm(1 - (Height$p/2))

data=fread("../Slope_Intercept/GWAS_Downloads/SZ/SCZ2_ckqny.scz2snpres.txt")
data$Z=qnorm(1 - (data$p/2))

BP=read.delim("../Slope_Intercept/GWAS_Downloads/BIP/pgc.bip.full.2012-04.txt", stringsAsFactors=F, sep="")
BP$Z=qnorm(1 - (BP$pval/2))

IBD=fread("../Slope_Intercept/GWAS_Downloads/IBD/EUR.IBD.gwas_info03_filtered.assoc")
IBD$Z=qnorm(1 - (IBD$P/2))

ADHD_C=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/ADHD Cross Disorder/pgc.cross.ADD4.2013-05.txt", stringsAsFactors=F, sep="")
ADHD_C$Z=qnorm(1 - (ADHD_C$pval/2))

Alzheimer=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Alzheimer/IGAP_stage_1.txt", stringsAsFactors=F)
Alzheimer$Z=qnorm(1 - (Alzheimer$Pvalue/2))

Autism_C=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Autism Cross Disorder/pgc.cross.AUT8.2013-05.txt", stringsAsFactors=F, sep="")
Autism_C$Z=qnorm(1 - (Autism_C$pval/2))

CAD=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Coronary Artery Disease/CARDIoGRAM_GWAS_RESULTS.txt", stringsAsFactors=F)
CAD$Z=qnorm(1 - (CAD$pvalue/2))

Menopause=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Menopause/Menopause_HapMap2_DayNG2015_18112015.txt.gz", stringsAsFactors=F)
Menopause$Z=qnorm(1 - (Menopause$p/2))

Smoke_CPD=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Tabacco/tag.cpd.tbl.gz", stringsAsFactors=F)
Smoke_CPD$Z=qnorm(1 - (Smoke_CPD$P/2))

MDD=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/MDD/pgc.mdd.full.2012-04.txt",stringsAsFactors=F, sep="")
MDD$Z=qnorm(1 - (MDD$pval/2))

Neuroticism=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Neuroticism/Neuroticism_Full.txt", stringsAsFactors=F)
Neuroticism$Z=qnorm(1 - (Neuroticism$P/2))

Fasting_Ins=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Fasting_Insulin/Insulin.txt", stringsAsFactors=F)
Fasting_Ins$Z=qnorm(1 - (Fasting_Ins$MainP/2))

BMI=read.delim("../Slope_Intercept/GWAS_Downloads/BMI/BMI_SNP_gwas_mc_merge_nogc.tbl.uniq.txt", stringsAsFactors=F)
BMI$Z=qnorm(1 - (BMI$p/2))

SWB=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/SWB_SSGAC/SWB_Full.txt.gz", stringsAsFactors=F)
SWB$Z=qnorm(1 - (SWB$Pval/2))

DS=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/DS_SSGAC/DS_Full.txt.gz" , stringsAsFactors=F)
DS$Z=qnorm(1 - (DS$Pval/2))

Menarche=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Menarche/Menarche_Nature2014_GWASMetaResults_17122014.txt", , stringsAsFactors=F)
Menarche$Z=qnorm(1 - (Menarche$GWAS_P/2))

ICV=fread("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/ICV/ENIGMA2_ICV_Combined_GenomeControlled_Jan23_v2.txt")
ICV$Z=qnorm(1 - (ICV$Pvalue/2))

CHIC=fread("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/CHIC_ssgac/CHIC_Summary_Benyamin2014.txt")
CHIC$Z=qnorm(1 - (CHIC$P/2))

Fasting_ProIns=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Fasting_Insulin/MAGIC_ln_fastingProinsulin.txt.gz", stringsAsFactors=F)
Fasting_ProIns$Z=qnorm(1 - (Fasting_ProIns$pvalue/2))

BW=fread("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Brith_weight/BW3_EUR_summary_stats.txt")
BW$Z=qnorm(1 - (BW$p/2))

BL=read.delim("../Slope_Intercept/GWAS_Downloads/Latest GWAS 27092016/Brith_weight/EGG-GWAS-BL.txt.gz",stringsAsFactors=F)
BL$Z=qnorm(1 - (BL$P/2))

#******************************************************************************************************************************************################


##Put observed Z values in the dataframe that contains expected chi-square stats
SNP_with_score$Height_Z=Height[(match(SNP_with_score$SNP2.Name,Height$MarkerName)),"Z"]
SNP_with_score$SZ_Z=data[(match(SNP_with_score$SNP2.Name,data$snpid)),Z]
SNP_with_score$BIP_Z=BP[(match(SNP_with_score$SNP2.Name,BP$snpid)),"Z"]
SNP_with_score$IBD_Z=IBD[(match(SNP_with_score$SNP2.Name,IBD$SNP)),Z]
SNP_with_score$ADHD_C_Z=ADHD_C[(match(SNP_with_score$SNP2.Name,ADHD_C$snpid)),"Z"]
SNP_with_score$Alzheimer_Z=Alzheimer[(match(SNP_with_score$SNP2.Name,Alzheimer$MarkerName)),"Z"]
SNP_with_score$Autism_C_Z=Autism_C[(match(SNP_with_score$SNP2.Name,Autism_C$snpid)),"Z"]
SNP_with_score$CAD_Z=CAD[(match(SNP_with_score$SNP2.Name,CAD$SNP)),"Z"]
SNP_with_score$Menopause_Z=Menopause[(match(SNP_with_score$SNP2.Name,Menopause$MarkerName)),"Z"]
SNP_with_score$Smoke_CPD_Z=Smoke_CPD[(match(SNP_with_score$SNP2.Name,Smoke_CPD$SNP)),"Z"]
SNP_with_score$MDD_Z=MDD[(match(SNP_with_score$SNP2.Name,MDD$snpid)),"Z"]
SNP_with_score$Neuroticism_Z=Neuroticism[(match(SNP_with_score$SNP2.Name,Neuroticism$SNP)),"Z"]
SNP_with_score$Fasting_Ins_Z=Fasting_Ins[(match(SNP_with_score$SNP2.Name,Fasting_Ins$Snp)),"Z"]
SNP_with_score$BMI_Z=BMI[(match(SNP_with_score$SNP2.Name,BMI$SNP)),"Z"]
SNP_with_score$SWB_Z=SWB[(match(SNP_with_score$SNP2.Name,SWB$MarkerName)),"Z"]
SNP_with_score$DS_Z=DS[(match(SNP_with_score$SNP2.Name,DS$MarkerName)),"Z"]
SNP_with_score$Menarche_Z=Menarche[(match(SNP_with_score$SNP2.Name,Menarche$MarkerName)),"Z"]
SNP_with_score$ICV_Z=ICV[(match(SNP_with_score$SNP2.Name,ICV$RSID)),Z]
SNP_with_score$CHIC_Z=CHIC[(match(SNP_with_score$SNP2.Name,CHIC$SNP)),Z]
SNP_with_score$Fasting_ProIns_Z=Fasting_ProIns[(match(SNP_with_score$SNP2.Name,Fasting_ProIns$snp)),"Z"]
SNP_with_score$BW_Z=BW[(match(SNP_with_score$SNP2.Name,BW$rsid)),Z]
SNP_with_score$BL_Z=BL[(match(SNP_with_score$SNP2.Name,BL$RSID)),"Z"]



##making enrichemnt data frame
Enrichment_frame=data.frame(((SNP_with_score$Height_Z^2)/SNP_with_score$Height),((SNP_with_score$SZ_Z^2)/SNP_with_score$SZ),((SNP_with_score$BIP_Z^2)/SNP_with_score$BIP), ((SNP_with_score$IBD_Z^2)/SNP_with_score$IBD), ((SNP_with_score$ADHD_C_Z^2)/SNP_with_score$ADHD_C),((SNP_with_score$Alzheimer_Z^2)/SNP_with_score$Alzheimer), ((SNP_with_score$Autism_C_Z^2)/SNP_with_score$Autism_C), ((SNP_with_score$CAD_Z^2)/SNP_with_score$CAD), ((SNP_with_score$Menopause_Z^2)/SNP_with_score$Menopause), ((SNP_with_score$Smoke_CPD_Z^2)/SNP_with_score$Smoke_CPD), ((SNP_with_score$MDD_Z^2)/SNP_with_score$MDD), ((SNP_with_score$Neuroticism_Z^2)/SNP_with_score$Neuroticism), ((SNP_with_score$Fasting_Ins_Z^2)/SNP_with_score$Fasting_Ins), ((SNP_with_score$BMI_Z^2)/SNP_with_score$BMI), ((SNP_with_score$SWB_Z^2)/SNP_with_score$SWB), ((SNP_with_score$DS_Z^2)/SNP_with_score$DS), ((SNP_with_score$Menarche_Z^2)/SNP_with_score$Menarche),((SNP_with_score$ICV_Z^2)/SNP_with_score$ICV),((SNP_with_score$CHIC_Z^2)/SNP_with_score$CHIC),((SNP_with_score$Fasting_ProIns_Z^2)/SNP_with_score$F_PIns),((SNP_with_score$BW_Z^2)/SNP_with_score$BW),((SNP_with_score$BL_Z^2)/SNP_with_score$BL))

rownames(Enrichment_frame)=SNP_with_score$SNP2.Name
#Making "Inf" value to NA
SNP_with_score[4,"Height_Z"]=NA

dat <- as.matrix(Enrichment_frame)
colnames(dat) <- c("Height","SZ","BIP","IBD", "ADHD_C", "Alzheimer", "Autism_C", "CAD", "Menopause", "Smoke_CPD", "MDD", "Neuroticism", "Fasting_Ins", "BMI", "SWB", "DS", "Menarche", "ICV", "CHIC", "Fast_PIns", "BW", "BL")


dat2=log(dat)
dat2[(dat2==Inf)]=NA
dat3=dat2


##LD-aware enrichment test for each SNP
lo=list()
for (i in 1:22){lo[[i]]=pchisq(SNP_with_score[,i+24]^2, df=Slopes[i,3], ncp=(Slopes[i,2]*Slopes[i,4]*SNP_with_score$L2/1173569), lower.tail = FALSE, log.p = FALSE)}
	
##multiple correction bonferroni
Pval_mat=do.call(rbind,lapply(lo,matrix,ncol=105,byrow=TRUE))
Pval_mat2=(t(apply(Pval_mat,1,function(x) {p.adjust(x,method="bonferroni")})))
Pval_mat2[(Pval_mat2<=0.05)]="*"
Pval_mat2[(!Pval_mat2=="*")]=""

##These will be used later just to output how many SNPs are below threshold
sig_chi_0.05=(lapply(lo,function(x)length(which(x<0.05))))
sig_chi_0.01=(lapply(lo,function(x)length(which(x<0.01))))


colnames(dat3)=c("Height","Schizophrenia","Bipolar disorder","Inflammatory bowel disease", "ADHD", "Alzheimer", "Autism", "Coronary artery disease", "Menopause", "Smoking (cigarettes per day)", "Major depressive disorder", "Neuroticism", "Fasting insulin", "Body mass index", "Subjective well-being", "Depressive symptoms", "Age at menarche", "Intracraneal volume", "Childhood intelligence", "Fasting proinsulin", "Birth weight", "Birth length")

pal=colorRampPalette(c("blue","white","red"))
col_my=pal(18)
    
##Exclude UC CD 
pdf("Enrichment_Heatmap_log_Z_from_P_with_stars_adjusted_0.05.pdf", width=30, height=10)
heatmap.2(t(dat3),colCol=c(rep("black",90),"red",rep("black",2),"red"),scale="none", cellnote=Pval_mat2, notecex=3, notecol="black", col=col_my,key=TRUE, symkey=FALSE, density.info="none", trace="none", cexCol=1.5,cexRow=1.5, dendrogram="row", Rowv=T, Colv=F, keysize=1,margin=c(12,18))
   
dev.off()







##LD-aware enrichment test for trait

lo_all=list()
for (i in 1:22){lo_all[[i]]=pchisq(sum(SNP_with_score[,i+24]^2, na.rm=T), df=(Slopes[i,3]* length(which(!is.na(SNP_with_score[,i+24]^2)))), ncp=sum(Slopes[i,2]*Slopes[i,4]*SNP_with_score[(which(!is.na(SNP_with_score[,i+24]^2))),2]/1173569), lower.tail = FALSE)}

names(lo_all)=colnames(dat)
enrich_all=data.frame(Pvalue=unlist(lo_all))
enrich_all$Bonferroni=p.adjust(enrich_all$Pvalue, method="bonferroni")
enrich_all$Number_of_SNPs_below_0.05 <- unlist(sig_chi_0.05)
enrich_all$Number_of_SNPs_below_0.01 <- unlist(sig_chi_0.01)
enrich_all$Total_number_of_SNPs=unlist(lapply(lo,function(x)length(which(!is.na(x)))))
enrich_all$Bonferroni=p.adjust(enrich_all$Pvalue, method="bonferroni")

write.table(enrich_all, file="Enrichment_trait_overall_pvalues.txt", sep="\t", quote=F)




