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


## Questions?
-------------
For more information, please visit the [Wiki](https://github.com/tbrunetti/GP3/wiki) on or contact me (tbrunetti) and I would be happy to address any issues. 





