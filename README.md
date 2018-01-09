# GP3
**G**WAS **P**re-**P**rocessing **P**ipeline

## Overview and Purpose
------------------------
An automated pipeline to pre-process GWAS data to determine samples to remove prior to input into imputation pipelines and association testing pipelines.  This should be used after initial round of QC/filtering has been performed (i.e. removing SNPs and samples that fail due to poor snp/sample quality from idats).

<p align="center">
<img src="https://github.com/tbrunetti/GP3/blob/master/images/GP3_pipeline_workflow.png" />
</p>

## Software Requirements
------------------------
The following are the minimum software requirements:
* Python version 2.7 (https://www.python.org/)
* R version 3.2 or better (https://cran.r-project.org/)
* PLINK version 1.9 or better (https://www.cog-genomics.org/plink2)
* KING software package (http://people.virginia.edu/~wc9c/KING/manual.html)
* chunkypipes (http://chunky-pipes.readthedocs.io/en/stable/getting_started.html)
* virtualenv -- ONLY required if admin rights are not granted (https://virtualenv.pypa.io/en/stable/) 


__*--R libraries that need to be installed manually--*__  
The following list of R libraries, including their dependencies must be installed and functional:  
  * GENESIS (http://bioconductor.org/packages/release/bioc/html/GENESIS.html)
  * GWASTools (http://bioconductor.org/packages/release/bioc/html/GWASTools.html)


__*--Software Requirements that can be installed automatically--*__  
The following list of Python libraries are required but the pipeline can automatically install them if pip is available:
  * SciPy stack, in particular the following packages: (https://scipy.org/)
    * pandas
    * numpy
    * matplotlib
  * Statistics (https://pypi.python.org/pypi/statistics)
  * Pillow (https://pypi.python.org/pypi/Pillow/3.4.2)
  * pyFPDF (http://pyfpdf.readthedocs.io/en/latest/index.html)
  * pyPDF2 (https://pypi.python.org/pypi/PyPDF2/1.26.0)


## User Generated File Requirements
-----------------------------------
There are two files that are minimally required in order to run the pipeline:
* Input PLINK file either in .bed or .ped format
* Populated sample_sheet_template.xlsx  



## Installation of virtual environment, chunkypipes, and pipeline
------------------------------------------------------------------
Please [click here](https://github.com/tbrunetti/GP3/blob/master/Installation_Instructions.md) for detailed instructions on setting up a virtual environment for shared systems or for installation on systems with root privledges.


__ALREADY INSTALLED CHUNKYPIPES AND PIPELINE?__  Click [here](https://github.com/tbrunetti/GP3/blob/master/Quick_Start_FAQs.md) for quick start.



## Output Files
---------------
If you navigate to your output directory you should notice a new directory matching the project name you specified at the time of the run.  Navigate into this directory and you should see new directories based on the ethnic group names you specified in your sample_sheet.xlsx as well as a set of PDFs.  These PDFs are the ones promised above in the diagram. If navigate into one of the directories of your ethnic group you will notice several PLINK files that were generated at each step of the pipeline.  
  
Addtionally, here are the notable final files of interest if the __--TGP__ flag is specified:  
  1. \<ethnic group name\>_all_samples_to_remove_from_original.txt
  2. \<ethnic group name\>_all_steps_completed_TGP_final followed by the following suffixes:
    * .bed
    * .bim
    * .fam
    * .kin
    * .kin0
    * .gds
  3. \<ethnic group name\>_all_steps_completed_TGP_final_GENESIS.Rdata
  4. \<ethnic group name\>_all_steps_completed_TGP_final_phenoGENESIS.txt  
  5. \<ethnic group name\>_TGP_PCA_plots.pdf
  

If __no --TGP flag is specified__ here are the final output file names:  
  1. \<ethnic group name\>_all_samples_to_remove_from_original.txt
  2. \<ethnic group name\>_all_steps_completed_final followed by the following suffixes:
    * .bed
    * .bim
    * .fam
    * .kin
    * .kin0
    * .gds
  3. \<ethnic group name\>_all_steps_completed_TGP_final_GENESIS.Rdata
  4. \<ethnic group name\>_all_steps_completed_TGP_final_phenoGENESIS.txt
  5. \<ethnic group name\>_all_steps_completed_final_GENESIS_sample_key_file.txt  
  6. \<ethnic group name\>_individual_PCA_plots.pdf


## Questions?
-------------
For more information, please visit the [Wiki](https://github.com/tbrunetti/GP3/wiki) on or contact me (tbrunetti) and I would be happy to address any issues. 





