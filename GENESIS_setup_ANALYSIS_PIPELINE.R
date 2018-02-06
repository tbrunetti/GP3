#source("https://bioconductor.org/biocLite.R")
#biocLite("GWASTools")

# read in command line arugments -- POSITIONAL!
args <- commandArgs(trailingOnly = T)
filePrefix <- args[1] # for KING and PLINK files
phenoFile <- args[2] # specially designed for GENESIS (essentially the 2nd and 6th column of fam file; need full name and path)
pcmat_num <- as.integer(args[4]) # used to calculate pcrelate
print(pcmat_num)
fullKin <- args[5]
library(GWASTools, lib.loc=args[3])

###Read in Barbados genotype data for PCAs
#biocLite("SNPRelate")
library("SNPRelate", lib.loc=args[3])
snpgdsBED2GDS(bed.fn = paste(filePrefix, ".bed", sep=""), bim.fn = paste(filePrefix, ".bim", sep=""), fam.fn = paste(filePrefix, ".fam", sep=""), 
              out.gdsfn = paste(filePrefix, ".gds", sep=""))

#biocLite("GENESIS")
library("GENESIS", lib.loc=args[3])
file.kin0 <- paste(filePrefix, ".kin0", sep="")
file.kin <- paste(filePrefix, ".kin", sep="")
geno <- GdsGenotypeReader(filename = paste(filePrefix, ".gds", sep=""))
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)

###Run PC analysis
#use .kin and .kin0 file
if(fullKin == "TRUE"){
  print("Both .kin0 and .kin will be used to generate pcair and pcrelate metrics")
  Kingmat <- king2mat(file.kin0=file.kin0,file.kin=file.kin,type="kinship",iids = iids)
  mypcair <- pcair(genoData = genoData, kinMat = Kingmat,divMat = Kingmat, v=200)
  mypcrel <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:pcmat_num],training.set = mypcair$unrels)
  #pcMat is not the number of PCs you have but instead the number of different admixture populations
}else{
  print("Only .kin0 will be used to generate pcair and pcrelate metrics")
  Kingmat <- king2mat(file.kin0=file.kin0,file.kin=NULL,type="kinship",iids = iids)
  mypcair <- pcair(genoData = genoData, kinMat = Kingmat,divMat = Kingmat, v=200)
  mypcrel <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:pcmat_num],training.set = mypcair$unrels)
  #pcMat is not the number of PCs you have but instead the number of different admixture populations 
}

pheno <- as.vector(as.matrix(read.table(phenoFile,header=F,na.string="NA")['V2']))
pheno <- pheno - 1

scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = mypcrel$sample.id,pc1 = mypcair$vectors[,1],pc2 = mypcair$vectors[,2],
	pc3 = mypcair$vectors[,3],pc4 = mypcair$vectors[,4],pc5 = mypcair$vectors[,5],pc6 = mypcair$vectors[,6],
	pc7 = mypcair$vectors[,7],pc8 = mypcair$vectors[,8], pc9 = mypcair$vectors[,9],
	pc10 = mypcair$vectors[,10],pc11 = mypcair$vectors[,11],pc12 = mypcair$vectors[,12],
	pc13 = mypcair$vectors[,13],pc14 = mypcair$vectors[,14],pc15 = mypcair$vectors[,15],
	pc16 = mypcair$vectors[,16],pc17 = mypcair$vectors[,17],pc18 = mypcair$vectors[,18],
	pc19 = mypcair$vectors[,19],pc20 = mypcair$vectors[,20],pc21 = mypcair$vectors[,21],
	pc22 = mypcair$vectors[,22],pc23 = mypcair$vectors[,23],pc24 = mypcair$vectors[,24],
	pc25 = mypcair$vectors[,25],pc26 = mypcair$vectors[,26],pc27 = mypcair$vectors[,27],
	pc28 = mypcair$vectors[,28],pc29 = mypcair$vectors[,29],pc30 = mypcair$vectors[,30],
	pc31 = mypcair$vectors[,31],pc32 = mypcair$vectors[,32],pc33 = mypcair$vectors[,33],
	pc34 = mypcair$vectors[,34],pc35 = mypcair$vectors[,35],pc36 = mypcair$vectors[,36],
	pc37 = mypcair$vectors[,37],pc38 = mypcair$vectors[,38],pc39 = mypcair$vectors[,39],
	pc40 = mypcair$vectors[,40],pc41 = mypcair$vectors[,41],pc42 = mypcair$vectors[,42],
	pc43 = mypcair$vectors[,43],pc44 = mypcair$vectors[,44],pc45 = mypcair$vectors[,45],
	pc46 = mypcair$vectors[,46],pc47 = mypcair$vectors[,47],pc48 = mypcair$vectors[,48],
	pc49 = mypcair$vectors[,49],pc50 = mypcair$vectors[,50],pheno = pheno))

covMatList <- list("Kin" = pcrelateMakeGRM(mypcrel))

# creates a binary file -- can open in R with load()
save(scanAnnot, covMatList, mypcair, file = paste(filePrefix, '_GENESIS', sep=""))
