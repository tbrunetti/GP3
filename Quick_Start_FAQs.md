## Running Pipeline with Defaults

If you run the minimum number of requirements like so:
```
chunky run run_GWAS_analysis_pipeline.py -inputPLINK <myPlink.bed/.ped> -phenoFile <sample_sheet_template.xlsx>
```  

This will run through all the default parameters of the pipeline from HWE through graphing PCAs (noTGP).  By running the above command the following flow will occur:  
1.  Output directory is current working directory and the project name is set to date and time stamp **[--outDir, --projectName]**
2.  Split out all ethnicies separately (the subsequent steps are recalculated for each individual ethnic group listed in the sample_sheet_template.xlsx)  
**FOR EACH ETHNIC GROUP CALCULATE AND FILTER:**  
3.  remove snps below 97% call rate **[--snpMiss]**  
4.  remove snps with HWE less than 1e-6 **[--hweThresh]**  
5.  remove samples below 97% call rate **[--sampleMiss]**  
6.  LD prune snps using method *indep* with VIF of 2 and window size of 50, step size of 5 **[--LDmethod, --VIF, --rsq, --windowSize, --stepSize]**  
7.  remove snps with MAF <= 0.05 **[--maf]**  
8.  remove samples with method *meanStd* with  excess heterozygosity of more than 3 standard deviations from the mean for each population **[--hetMethod, --hetTresh, --hetThreshMin, --het\_std]**  
9.  perform IBD and remove true duplicate samples  
10. Remove outliers specified by user (if available) **[--outliers]**   
11. PCA **without** TGP for each group using GENESIS (includes running KING) and using admixture population number estimation of 5 **[--pcmat]**  
12. graph each group with PCA screen plots and boxplots with mean and standard deviations centered around each input cohort **[--centerPop]**  
13. generate a list of samples to remove  


## Running Pipeline with TGP with Defaults

You can run the default steps above with TGP by specifying the **--TGP** flag at the command line prompt:
```
chunky run run_GWAS_analysis_pipeline.py -inputPLINK <myPlink.bed/.ped> -phenoFile <sample_sheet_template.xlsx> --TGP
```  

This will run through all the default parameters of the pipeline from HWE through graphing PCAs (only with TGP).  By running the above command the following flow will occur:  
1.  Output directory is current working directory and the project name is set to data and time stamp **[--outDir, --projectName]**  
2.  Split out all ethnicies separately (the subsequent steps are recalculated for each individual ethnic group listed in the sample_sheet_template.xlsx)  
**FOR EACH ETHNIC GROUP CALCULATE AND FILTER:**  
3.  remove snps below 97% call rate **[--snpMiss]**  
4.  remove snps with HWE less than 1e-6 **[--hweThresh]**  
5.  remove samples below 97% call rate **[--sampleMiss]**  
6.  LD prune snps using method *indep* with VIF of 2 and window size of 50, step size of 5 **[--LDmethod, --VIF, --rsq, --windowSize, --stepSize]**  
7.  remove snps with MAF <= 0.05 **[--maf]**  
8.  remove samples with method *meanStd* with  excess heterozygosity of more than 3 standard deviations from the mean for each population **[--hetMethod, --hetTresh, --hetThreshMin, --het\_std]**  
9.  perform IBD and remove true duplicate samples  
10. Remove outliers specified by user (if available) **[--outliers]**  
11. Merge TGP samples with input samples and only take common snps between groups **[--TGP]**  
12. PCA **including** TGP for each group using GENESIS (includes running KING) and using admixture population number estimation of 5 **[--pcmat]**  
13. graph each group with PCA screen plots and boxplots with mean and standard deviations centered around each input cohort **[--centerPop]**   
14. generate a list of samples to remove  



## Pipeline Arguments

| Argument | Usage | Type | Default | Explanation |
| :---: | :---: | :---: | :---: | :---: |
| -inputPLINK | required | str ending in .bed/.ped |NA | Full path to PLINK file ends in .bed or .ped, __whitespace characters are not allowed__ |
| -phenoFile | required | NA | NA | Full path to populated phenotype file: sample_sheet_template.xlsx | 
| --config | only required if .json not in .home of chunky | str | search in chunky .home directory | configuration file produced after running `chunky config run_GWAS_analysis_pipeline.py` | 
| --outDir | optional | str | current working directory | Full path to an already existing directory or location where you would like GP3 to build the project | 
| --projectName | optional | str | year-month-date-hour-min-sec | Name of project to be created in the outDir location __whitespace characters are not allowed__ | 
| --startStep | optional | str | hwe | str | options: if --TGP not specified -> [hwe, LD, maf, het, ibd, PCA_indi, or PCA_indi_graph] if --TGP is set ->[hwe, LD, maf, het, ibd, outlier_removal, PCA_TGP, or PCA_TGP_graph]] The part of the pipeline you would like to start with | 
| --endStep | optional | str | PCA_indi_graph (noTGP) or PCA_TGP_graph (if --TGP used) | options: if --TGP not set -> [hwe, LD, maf, het, ibd, PCA_indi, or PCA_indi_graph] if --TGP is set -> [hwe, LD, maf, het, ibd, outlier_removal, PCA_TGP, or PCA_TGP_graph]] Point of the pipeline where you would like to stop analysis.  __This step is inclusive!__ | 
| --hweThresh | optional | float | 1e-6 | Filters out SNPs that are smaller than this threshold due to liklihood of genotyping error | 
| --LDmethod | optional | str | indep | options:[indep, indep-pairwise or indep-pairphase] Method to calculate linkage disequilibrium.  See PLINK documentation for more information. | 
| --VIF | optional | int | 2 | variant inflation factor for indep method LD pruning method only; indep-pairwise or indep-pairphase method will not use VIF | 
| --rsq | optional | float | 0.50 | any floating point number between 0.0-1.0; r-squared threshold for indep-pairwise or indep-pairphase LD pruning method; indep method will not use rsq | 
| --windowSize | optional | int | 50 | any integer; the window size in kb for LD analysis | 
| --stepSize | optional | int | 5 | any integer; variant count to shift window after each iteration | 
| --maf | optional | float | 0.05 | any floating point number between 0.0-1.0; filter remaining LD pruned variants by MAF, any MAF below set threshold is filtered out | 
| --hetMethod | optional | str | meanStd | options:  minMax or meanStd; method to use to determine heterozygosity.  minMax filter based on the parameters --hetThresh as the max F-inbreeding coefficient and --hetThreshMin for the minimum F-inbreeding coeffient, which by default are 0.10 and -0.10, respectively.  The meanStd filter method calculates a het_score: 1-[observed[HOM]/total] and then filters out any samples that are more than 3 std deviations from the mean het_score. The number of standard deviations from the mean can be changed using the --het_std parameter | 
| --het_std | optional | int or float | 3 | any floating point number or integer; if using hetMethod=meanStd you can determine how many standard deviations aways from the mean is allowable for heterozygosity. Setting to 3 is interpreted as +/-3 standard deviations away from the mean of the het_score, calculated as 1-[observed(HOM)/total] | 
| --hetThresh | optional | float | 0.10 | any floating point number; filter out samples where inbreeding coefficient is greater than threshold (heterozygosity filtering); only used when method minMax for --hetMethod is selected |
| --hetThreshMin | optional | float | -0.10 | any floating point number; filter out samples where inbreeding coefficient is samller than min threshold set (heterozygosity filtering); only used when method minMax for hetThresh is selected | 
| --sampleMiss | optional | float | 0.03 | any floating point number between 0.0-1.0; Maximum missingness of genotype call in sample before it should be filtered out. Where 0 is no missing, and 1 is all missing (0.03 is interpreted as 3 percent of snp calls are missing in a sample) | 
| --snpMiss | optional | float | 0.03 | any floating point number between 0.0-1.0; Maximum missingness of genotype call in a SNP cluster before the SNP should be filtered out.  Where 0 is no missing, and 1 is all missing (0.03 is interpreted as 3 percent of sample calls are missing in a snp) | 
| --TGP | optional | flag | NA | specifying this flag means to generate PCA plots with TGP data merged into the given cohort data set for the 5 superpopulations in TGP (AFR, AMR, EAS, EUR, SAS)
| --centerPop | optional | str | myGroup | options: literally the string myGroup or available TGP group merged into input dataset; when using the TGP flag, you have the option to specify which population cohort that PCs should be centered around for boxplots.  By default this is set to your group(s) listed in the sample sheet.  You can pick a TGP super population listed in the TGP_Sub_and_SuperPopulation_info.txt file. __CASE SENSITIVE!__ | 
| --outliers | optional | str | None | A txt file of FID and IID, tab-delimited and one sample per line, that are outliers that should be removed from the sample set (PCA outlier removal); Use original names (original FID and IID), not renamed 1-n for GENESIS formatting | 
| --pcmat | optional | int | 5 | any integer; Number of predicted admixture populations in dataset to be used in GENESIS calculation for PCA | 
| --reanalyze | optional | flag | NA | by adding this flag, it means you are going to pass a dataset through the pipeline that has already been partially/fully analyzed by this pipeline. WARNING! May over write exisiting data!! __required if using --startStep argument OR if using --endStep arguments on an already existing project__ |  


## Pipeline Help Flag

For more help you can use the **-h** or **--help** flag at the command line prompt and it will give you a detailed list of arguments that can be used:
```
chunky run run_GWAS_analysis_pipeline.py -h
```
The brackets in bold following the step numbers listed above, indicate the arguments you can call at the command prompt to change the default values and default methods.  The **-h** or **--help** will list all the input options for using and changing these parameters.

