#source("https://bioconductor.org/biocLite.R")
#biocLite("GENESIS")
args <- commandArgs(trailingOnly = T)
load_file <- args[1] # GENESIS PCA R object scanAnnot dataframe
cohort <- args[2]
outDir_path <- args[4]
library("GENESIS", lib.loc = args[3])


# loaded from genesis setup script in GWAS_analysis_pipeline.py
load(paste(load_file))
get_dataframe <- pData(scanAnnot)

population_name = cohort


pdf(paste(outDir_path, population_name, "_individual_PCA_plots.pdf", sep=''),onefile=T)
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

# Assign group name
pcs_with_pops <- get_dataframe
pcs_with_pops$population <- population_name

# color order for scatter plots
cols=c("blue")

# scatter PC1 vs PC2
plot(pcs_with_pops$pc1, pcs_with_pops$pc2, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 1",ylab="Principal Component 2")
abline(v=mean(pcs_with_pops$pc1), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc1)+(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(6*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(6*sd(pcs_with_pops$pc1)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc2), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc2)+(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(6*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(6*sd(pcs_with_pops$pc2)), col='darkgrey')
legend("bottomright", population_name ,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)

# scatter PC1 vs PC3
plot(pcs_with_pops$pc1, pcs_with_pops$pc3, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 1",ylab="Principal Component 3")
abline(v=mean(pcs_with_pops$pc1), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc1)+(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(6*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(6*sd(pcs_with_pops$pc1)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc3), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc3)+(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(6*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(6*sd(pcs_with_pops$pc3)), col='darkgrey')
legend("topright", population_name ,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)

# scatter PC1 vs PC4
plot(pcs_with_pops$pc1, pcs_with_pops$pc4, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 1",ylab="Principal Component 4")
abline(v=mean(pcs_with_pops$pc1), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc1)+(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(6*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(6*sd(pcs_with_pops$pc1)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc4), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc4)+(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(6*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(6*sd(pcs_with_pops$pc4)), col='darkgrey')

legend("bottomright", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC1 vs PC5
plot(pcs_with_pops$pc1, pcs_with_pops$pc5, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 1",ylab="Principal Component 5")
abline(v=mean(pcs_with_pops$pc1), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc1)+(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(6*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(6*sd(pcs_with_pops$pc1)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc5), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc5)+(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(6*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(6*sd(pcs_with_pops$pc5)), col='darkgrey')
legend("bottomright", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)

# scatter PC2 vs PC3
plot(pcs_with_pops$pc2, pcs_with_pops$pc3, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 2",ylab="Principal Component 3")
abline(v=mean(pcs_with_pops$pc2), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc2)+(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(6*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(6*sd(pcs_with_pops$pc2)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc3), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc3)+(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(6*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(6*sd(pcs_with_pops$pc3)), col='darkgrey')

legend("topleft", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC2 vs PC4
plot(pcs_with_pops$pc2, pcs_with_pops$pc4, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 2",ylab="Principal Component 4")
abline(v=mean(pcs_with_pops$pc2), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc2)+(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(6*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(6*sd(pcs_with_pops$pc2)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc4), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc4)+(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(6*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(6*sd(pcs_with_pops$pc4)), col='darkgrey')
legend("bottomleft", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC2 vs PC5
plot(pcs_with_pops$pc2, pcs_with_pops$pc5, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 2",ylab="Principal Component 5")
abline(v=mean(pcs_with_pops$pc2), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc2)+(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(6*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(6*sd(pcs_with_pops$pc2)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc5), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc5)+(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(6*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(6*sd(pcs_with_pops$pc5)), col='darkgrey')
legend("bottomleft", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC3 vs PC4
plot(pcs_with_pops$pc3, pcs_with_pops$pc4, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 3",ylab="Principal Component 4")
abline(v=mean(pcs_with_pops$pc3), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc3)+(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(6*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(6*sd(pcs_with_pops$pc3)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc4), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc4)+(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(6*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(6*sd(pcs_with_pops$pc4)), col='darkgrey')
legend("bottomleft", population_name, pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC3 vs PC5
plot(pcs_with_pops$pc3, pcs_with_pops$pc5, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 3",ylab="Principal Component 5")
abline(v=mean(pcs_with_pops$pc3), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc3)+(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(6*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(6*sd(pcs_with_pops$pc3)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc5), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc5)+(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(6*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(6*sd(pcs_with_pops$pc5)), col='darkgrey')
legend("bottomleft", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)




#Boxplots centered on mean African ancestry
PCA_pop_for_centering = subset(pcs_with_pops,population==population_name)


#boxplot for PC1
m11=mean(PCA_pop_for_centering$pc1)
s11=sd(PCA_pop_for_centering$pc1)
boxplot(pc1 ~ population, data = pcs_with_pops,ylab = "PC1", xlab=paste("PC1, +/- 6 SD on",population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc1 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m11, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m11-(6*s11), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m11+(6*s11), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC2
m12=mean(PCA_pop_for_centering$pc2)
s12=sd(PCA_pop_for_centering$pc2)
boxplot(pc2 ~ population, data = pcs_with_pops,ylab = "PC2", xlab=paste("PC2, +/- 6 SD on", population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc2 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m12, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m12-(6*s12), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m12+(6*s12), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC3
m13=mean(PCA_pop_for_centering$pc3)
s13=sd(PCA_pop_for_centering$pc3)
boxplot(pc3 ~ population, data = pcs_with_pops,ylab = "PC3", xlab=paste("PC3, +/- 6 SD on", population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc3 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m13, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m13-(6*s13), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m13+(6*s13), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC4
m14=mean(PCA_pop_for_centering$pc4)
s14=sd(PCA_pop_for_centering$pc4)
boxplot(pc4 ~ population, data = pcs_with_pops,ylab = "PC4", xlab=paste("PC4, +/- 6 SD on", population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc4 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m14, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m14-(6*s14), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m14+(6*s14), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC5
m15=mean(PCA_pop_for_centering$pc5)
s15=sd(PCA_pop_for_centering$pc5)
boxplot(pc5 ~ population, data = pcs_with_pops,ylab = "PC5", xlab=paste("PC5, +/- 6 SD on", population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc5 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m15, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m15-(6*s15), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m15+(6*s15), col = "black", lty=2, lwd=1,xpd=F)


#Determine outliers
indi_cohort = subset(pcs_with_pops,population==population_name)
indi_cohort$PC1_outlier[as.numeric(indi_cohort$pc1)<(m11-(6*s11)) | as.numeric(indi_cohort$pc1)>(m11+(6*s11))] <- 1
indi_cohort$PC1_outlier[as.numeric(indi_cohort$pc1)>=(m11-(6*s11)) & as.numeric(indi_cohort$pc1)<=(m11+(6*s11))] <- 0
indi_cohort$PC2_outlier[as.numeric(indi_cohort$pc2)<(m12-(10*s12)) | as.numeric(indi_cohort$pc2)>(m12+(10*s12))] <- 1
indi_cohort$PC2_outlier[as.numeric(indi_cohort$pc2)>=(m12-(10*s12)) & as.numeric(indi_cohort$pc2)<=(m12+(10*s12))] <- 0
indi_cohort$PC3_outlier[as.numeric(indi_cohort$pc3)<(m13-(6*s13)) | as.numeric(indi_cohort$pc3)>(m13+(6*s13))] <- 1
indi_cohort$PC3_outlier[as.numeric(indi_cohort$pc3)>=(m13-(6*s13)) & as.numeric(indi_cohort$pc3)<=(m13+(6*s13))] <- 0
indi_cohort$PC4_outlier[as.numeric(indi_cohort$pc4)<(m14-(6*s14)) | as.numeric(indi_cohort$pc4)>(m14+(6*s14))] <- 1
indi_cohort$PC4_outlier[as.numeric(indi_cohort$pc4)>=(m14-(6*s14)) & as.numeric(indi_cohort$pc4)<=(m14+(6*s14))] <- 0
indi_cohort$PC5_outlier[as.numeric(indi_cohort$pc5)<(m15-(6*s15)) | as.numeric(indi_cohort$pc5)>(m15+(6*s15))] <- 1
indi_cohort$PC5_outlier[as.numeric(indi_cohort$pc5)>=(m15-(6*s15)) & as.numeric(indi_cohort$pc5)<=(m15+(6*s15))] <- 0

indi_cohort_outliers = subset(indi_cohort,PC1_outlier==1 | PC2_outlier==1 | PC3_outlier==1 | PC4_outlier==1 | PC5_outlier==1)
write.table(indi_cohort_outliers,paste(outDir_path,population_name,"_outliers_based_on_individual_pcs.txt", sep=''),row.names=F,col.names=T,quote=F,sep="\t")

indi_cohort_outliers <- read.table(paste(outDir_path,population_name,"_outliers_based_on_individual_pcs.txt", sep=''),header=T,all=T,sep="\t")
indi_cohort_outliers2 = subset(indi_cohort_outliers,select=c("scanID","pheno","PC1_outlier","PC2_outlier"))

indi_cohort_outliers2$EXTRACOL <- 1
cohort_2 <- merge(indi_cohort,indi_cohort_outliers2,all=T)
cohort_fin <- cohort_2[is.na(cohort_2$EXTRACOL),]
cohort_fin$EXTRACOL <- NULL
cohort_final = subset(cohort_fin,select=c("scanID","pc1","pc2","pc3","pc4","pc5"))
write.table(cohort_final,paste(outDir_path,population_name,"_PCs_No_Outliers_based_on_individual_pcs.txt", sep=''),row.names=F,col.names=T,quote=F,sep="\t")


dev.off()

