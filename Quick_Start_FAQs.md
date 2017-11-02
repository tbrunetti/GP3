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
6.  LD prune snps using method *indep* with rsq of 0.50 and window size of 50, step size of 5 **[--LDmethod, --VIF, --rsq, --windowSize, --stepSize]**  
7.  remove snps with MAF <= 0.05 **[--maf]**  
8.  remove samples with method *minMax* with  excess heterozygosity of F-inbreeding coefficient > 0.10 or F-inbreeing coefficient < -0.10 **[--hetMethod, --hetTresh, --hetThreshMin, --het\_std]**  
9.  perform IBD and remove true duplicate samples  
10. PCA **without** TGP for each group using GENESIS (includes running KING)
11. graph each group with PCA screen plots and boxplots  
12. generate a list of samples to remove  


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
6.  LD prune snps using method *indep* with rsq of 0.50 and window size of 50, step size of 5 **[--LDmethod, --VIF, --rsq, --windowSize, --stepSize]**  
7.  remove snps with MAF <= 0.05 **[--maf]**  
8.  remove samples with method *minMax* with excess heterozygosity of F-inbreeding coefficient > 0.10 or F-inbreeding coefficient < -0.10 **[--hetMethod, --hetThresh, --hetThreshMin, --het\_std]**  
9.  perform IBD and remove true duplicate samples  
10. Merge TGP samples with input samples and only take common snps between groups **[--TGP]**  
11. PCA **including** TGP for each group using GENESIS (includes running KING)  
12. graph each group with PCA screen plots and boxplots with mean and standard deviations centered around each input cohort **[--centerPop]**   
13. generate a list of samples to remove


## Pipeline Help Flag

For more help you can use the **-h** or **--help** flag at the command line prompt and it will give you a detailed list of arguments that can be used:
```
chunky run run_GWAS_analysis_pipeline.py -h
```
The brackets in bold following the step numbers listed above, indicate the arguments you can call at the command prompt to change the default values and default methods.  The **-h** or **--help** will list all the input options for using and changing these parameters.

