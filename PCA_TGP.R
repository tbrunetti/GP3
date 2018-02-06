args <- commandArgs(trailingOnly = T)
load_file <- args[1] # GENESIS PCA R object scanAnnot dataframe
cohort <- args[2]
outDir_path <- args[4]
libraru("BiocGeneric", lib.loc = args[3])
library("Biobase", lib.loc = args[3])
library("GENESIS", lib.loc = args[3])
library("GWASTools", lib.loc = args[3])
TGP_pop_file = args[5] # full path including file name of sup_sub TGP population file
center_population = args[6] # must be a population (case sensitive) in the TGP population file or the cohort using


# loaded from genesis setup script in GWAS_analysis_pipeline.py
load(paste(load_file)) 
get_dataframe <- pData(scanAnnot)
population_name = cohort # must match one listed in tgp_combined table in population column
TGP_info <- read.table(paste(TGP_pop_file), header=T)
tgp_combined <-  merge(get_dataframe, TGP_info, all.x = T)
colnames(tgp_combined)[colnames(tgp_combined)=="Sup_pop"] <- "population"
tgp_combined$population <- as.character(tgp_combined$population)
tgp_combined$population[is.na(tgp_combined$population)] <- population_name


# generate PDF of PCA plots and outlier calculations
pdf(paste(outDir_path,population_name, "_TGP_PCA_plots.pdf", sep=''),onefile=T)
par(mfrow=c(1,1))
par(mar=c(5,4.5,0,2))
par(omi=c(.5,0,.5,0))
fig=c(3,9,1,4)/10

#Plot first 200 Eigenvalues from mypcair
data <- data.frame(as.matrix(mypcair$values))
data['eigenvalues'] <- data.frame(as.matrix(mypcair$values))
data['proportion_var_explained'] <-data.frame((as.matrix(mypcair$values/mypcair$sum.values)))
data['PC'] <- data.frame((as.numeric(rownames(data))))
subset_data <- subset(data,PC<201)
subset_data['new_eigensum'] <-data.frame(colSums(subset_data)['eigenvalues'])
subset_data['subset_Prop_variance'] <- data.frame(subset_data$eigenvalues/subset_data$new_eigensum)
plot(subset_data$PC, subset_data$subset_Prop_variance, col="black", pch=20, cex=1,ylab="Proportion of Variance Explained",xlab="PCs",xaxt="n")
axis(1,at=subset_data$PC,labels=subset_data$PC)

#Plot first 20 Eigenvalues from mypcair
data <- data.frame(as.matrix(mypcair$values))
data['eigenvalues'] <- data.frame(as.matrix(mypcair$values))
data['proportion_var_explained'] <-data.frame((as.matrix(mypcair$values/mypcair$sum.values)))
data['PC'] <- data.frame((as.numeric(rownames(data))))
subset_data <- subset(data,PC<21)
subset_data['new_eigensum'] <-data.frame(colSums(subset_data)['eigenvalues'])
subset_data['subset_Prop_variance'] <- data.frame(subset_data$eigenvalues/subset_data$new_eigensum)
plot(subset_data$PC, subset_data$subset_Prop_variance, col="black", pch=20, cex=1,ylab="Proportion of Variance Explained",xlab="PCs",xaxt="n")
axis(1,at=subset_data$PC,labels=subset_data$PC)

# Assign experimental vs TGP
tgp_combined$group[tgp_combined$population==population_name] <- population_name
tgp_combined$group[tgp_combined$population!=population_name] <- "TGP"
tgp_combined$group <- as.factor(tgp_combined$group)
tgp_combined$population <- as.factor(tgp_combined$population)

subset_dataframe = subset(tgp_combined, population==population_name)

cols = list()
potential_TGP_cols <- c("blue", "green", "orange", "yellow", "red", "brown")
color_order <- c(levels(as.factor(tgp_combined$population)))
index_of_cohort <- which(color_order %in% population_name)
total_pops <- length(color_order)
added_cohort = 0 # turns to 1 if cohort has been added to list so can utilize all potential TGP colors
# color order for scatter plots
for (i in 1:total_pops){
  if (i == index_of_cohort){ # population cohort will always be in black
    cols <- c(cols, "black")
    added_cohort <- added_cohort +1
  } else if (added_cohort ==1) {
    cols <- c(cols, potential_TGP_cols[i-1])
  } else{
    cols <- c(cols, potential_TGP_cols[i])
  }
}

cols <- unlist(cols)
# scatter PC1 vs PC2
plot(tgp_combined$pc1, tgp_combined$pc2, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 2")
legend("topright", levels(tgp_combined$population),pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc1, subset_dataframe$pc2, col='black', pch=17, cex=0.5)

# scatter PC1 vs PC3
plot(tgp_combined$pc1, tgp_combined$pc3, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 3")
legend("topright", levels(tgp_combined$population),pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc1, subset_dataframe$pc3, col='black', pch=17, cex=0.5)

# scatter PC1 vs PC4
plot(tgp_combined$pc1, tgp_combined$pc4, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 4")
legend("topright", levels(tgp_combined$population),pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc1, subset_dataframe$pc4, col='black', pch=17, cex=0.5)

# scatter PC1 vs PC5
plot(tgp_combined$pc1, tgp_combined$pc5, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 5")
legend("topright", levels(tgp_combined$population),pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc1, subset_dataframe$pc5, col='black', pch=17, cex=0.5)


# scatter PC2 vs PC3
plot(tgp_combined$pc2, tgp_combined$pc3, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 3")
legend("topright", levels(tgp_combined$population),pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc2, subset_dataframe$pc3, col='black', pch=17, cex=0.5)


# scatter PC2 vs PC4
plot(tgp_combined$pc2, tgp_combined$pc4, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 4")
legend("topright", levels(tgp_combined$population),pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc2, subset_dataframe$pc4, col='black', pch=17, cex=0.5)

# scatter PC2 vs PC5
plot(tgp_combined$pc2, tgp_combined$pc5, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 5")
legend("topright", levels(tgp_combined$population),pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc2, subset_dataframe$pc5, col='black', pch=17, cex=0.5)


# scatter PC3 vs PC4
plot(tgp_combined$pc3, tgp_combined$pc4, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 3",ylab="Principal Component 4")
legend("topright", levels(tgp_combined$population), pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc3, subset_dataframe$pc4, col='black', pch=17, cex=0.5)


# scatter PC3 vs PC5
plot(tgp_combined$pc3, tgp_combined$pc5, col=cols[as.numeric(tgp_combined$population)],
     pch=c(17,20)[as.numeric(tgp_combined$group)], cex=0.5,xlab="Principal Component 3",ylab="Principal Component 5")
legend("topright", levels(tgp_combined$population),pch=rep(c(21), total_pops),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", levels(tgp_combined$group),pch=c(17,20),pt.cex=0.5,cex=0.5)
points(subset_dataframe$pc3, subset_dataframe$pc5, col='black', pch=17, cex=0.5)



#Boxplots centered on mean African ancestry
PCA_pops_center = subset(tgp_combined,population==center_population)


#boxplot for PC1
m11=mean(PCA_pops_center$pc1)
s11=sd(PCA_pops_center$pc1)
boxplot(pc1 ~ population, data = tgp_combined,ylab = "PC1", xlab=paste("PC1, +/- 6 SD on",center_population, sep=" "),cex=0.5,outline=F)
stripchart(pc1 ~ population, data = tgp_combined,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m11, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m11-(6*s11), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m11+(6*s11), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC2
m12=mean(PCA_pops_center$pc2)
s12=sd(PCA_pops_center$pc2)
boxplot(pc2 ~ population, data = tgp_combined,ylab = "PC2", xlab=paste("PC2, +/- 6 SD on", center_population, sep=" "),cex=0.5,outline=F)
stripchart(pc2 ~ population, data = tgp_combined,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m12, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m12-(6*s12), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m12+(6*s12), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC3
m13=mean(PCA_pops_center$pc3)
s13=sd(PCA_pops_center$pc3)
boxplot(pc3 ~ population, data = tgp_combined,ylab = "PC3", xlab=paste("PC3, +/- 6 SD on", center_population, sep=" "),cex=0.5,outline=F)
stripchart(pc3 ~ population, data = tgp_combined,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m13, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m13-(6*s13), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m13+(6*s13), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC4
m14=mean(PCA_pops_center$pc4)
s14=sd(PCA_pops_center$pc4)
boxplot(pc4 ~ population, data = tgp_combined,ylab = "PC4", xlab=paste("PC4, +/- 6 SD on", center_population, sep=" "),cex=0.5,outline=F)
stripchart(pc4 ~ population, data = tgp_combined,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m14, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m14-(6*s14), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m14+(6*s14), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC5
m15=mean(PCA_pops_center$pc5)
s15=sd(PCA_pops_center$pc5)
boxplot(pc5 ~ population, data = tgp_combined,ylab = "PC5", xlab=paste("PC5, +/- 6 SD on", center_population, sep=" "),cex=0.5,outline=F)
stripchart(pc5 ~ population, data = tgp_combined,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m15, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m15-(6*s15), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m15+(6*s15), col = "black", lty=2, lwd=1,xpd=F)

dev.off()


#Determine outliers
outliers_only_cohort = subset(tgp_combined,population==population_name)
outliers_only_cohort$PC1_outlier[as.numeric(outliers_only_cohort$pc1)<(m11-(6*s11)) | as.numeric(outliers_only_cohort$pc1)>(m11+(6*s11))] <- 1
outliers_only_cohort$PC1_outlier[as.numeric(outliers_only_cohort$pc1)>=(m11-(6*s11)) & as.numeric(outliers_only_cohort$pc1)<=(m11+(6*s11))] <- 0
outliers_only_cohort$PC2_outlier[as.numeric(outliers_only_cohort$pc2)<(m12-(10*s12)) | as.numeric(outliers_only_cohort$pc2)>(m12+(10*s12))] <- 1
outliers_only_cohort$PC2_outlier[as.numeric(outliers_only_cohort$pc2)>=(m12-(10*s12)) & as.numeric(outliers_only_cohort$pc2)<=(m12+(10*s12))] <- 0
outliers_only_cohort$PC3_outlier[as.numeric(outliers_only_cohort$pc3)<(m13-(6*s13)) | as.numeric(outliers_only_cohort$pc3)>(m13+(6*s13))] <- 1
outliers_only_cohort$PC3_outlier[as.numeric(outliers_only_cohort$pc3)>=(m13-(6*s13)) & as.numeric(outliers_only_cohort$pc3)<=(m13+(6*s13))] <- 0
outliers_only_cohort$PC4_outlier[as.numeric(outliers_only_cohort$pc4)<(m14-(6*s14)) | as.numeric(outliers_only_cohort$pc4)>(m14+(6*s14))] <- 1
outliers_only_cohort$PC4_outlier[as.numeric(outliers_only_cohort$pc4)>=(m14-(6*s14)) & as.numeric(outliers_only_cohort$pc4)<=(m14+(6*s14))] <- 0
outliers_only_cohort$PC5_outlier[as.numeric(outliers_only_cohort$pc5)<(m15-(6*s15)) | as.numeric(outliers_only_cohort$pc5)>(m15+(6*s15))] <- 1
outliers_only_cohort$PC5_outlier[as.numeric(outliers_only_cohort$pc5)>=(m15-(6*s15)) & as.numeric(outliers_only_cohort$pc5)<=(m15+(6*s15))] <- 0

cohort_all_outliers = subset(outliers_only_cohort,PC1_outlier==1 | PC2_outlier==1 | PC3_outlier==1 | PC4_outlier==1 | PC5_outlier==1)
write.table(cohort_all_outliers,paste(outDir_path,population_name,"_outliers.txt", sep=''),row.names=F,col.names=T,quote=F,sep="\t")

cohort_all_outliers <- read.table(paste(outDir_path,population_name,"_outliers.txt", sep=''),header=T,all=T,sep="\t")
cohort_all_outliers2 = subset(cohort_all_outliers,select=c("scanID","pheno","PC1_outlier","PC2_outlier"))

cohort_all_outliers2$EXTRACOL <- 1
cohort2 <- merge(outliers_only_cohort,cohort_all_outliers2,all=T)
cohortfin <- cohort2[is.na(cohort2$EXTRACOL),]
cohortfin$EXTRACOL <- NULL
cohort_final = subset(cohortfin,select=c("scanID","pc1","pc2","pc3","pc4","pc5"))
write.table(cohort_final,paste(outDir_path,population_name,"_PCs_No_Outliers.txt", sep=''),row.names=F,col.names=T,quote=F,sep="\t")

