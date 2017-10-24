# GP3
GWAS Pre-Processing Pipeline

## Overview and Purpose
------------------------
An automated pipeline to pre-process GWAS data for input into imputation pipelines and association testing pipelines.  This should be used after initial round of QC/filtering has been performed (i.e. removing SNPs and samples that fail due to poor snp/sample quality from idats).


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
The following list of R libraries must be installed and functional:


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


## Installing and Running Pipeline on HPC or without sudo privileges
---------------------------------------------------------------------
1. Install, Create, Activate Virutal Environment
2. Clone Repository (only needs to be done once)
3. Activate chunkypipes, configure and install pipeline (only needs to be performed once)
4. Modifying the Slurm script
5. Run pipeline

### Installation, Creation, and Activation of Virtual Environment (installAtion and creation only needs to ever be done once)
------------------------------------------------------------------------------------------------------------------------------
For users running the pipeline on a HPC system or without administrative privileges, a virtual environment must be created in order to utilize the pipeline to its maximum capabilities and efficiency.  The reason being, certain packages and dependencies can be installed using the pipeline to provide the user with an easy-to-use experience.  Many of these packages will install on your system globally, which cannot be done on shared distributed file systems such as HPC or if administrative rights are not granted on the computer you are using.  To bypass this, you will need to create a Python Virtual Environment.  Below are the steps for installing and creating this environment on your HPC.  NOTE:  THIS ONLY EVER NEEDS TO BE DONE ONE TIME!

**1.  Download and Install virtualenv**
----------------------------------------
For the full documentation of of virtualenv please refer to the following website:  https://virtualenv.pypa.io/en/stable/installation/  Typing in the following commands below will download virtualenv from the Web and create a virtual environment call myVE under the virtualenv-15.0.0 directory.  The name myVE can be changed to any full path file name of your choosing so that the virtual environment does not need to be listed in the virtualenv-15.0.0 directory.

```
$ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-15.0.0.tar.gz
$ tar xvfz virtualenv-15.0.0.tar.gz
$ cd virtualenv-15.0.0
$ python virtualenv.py myVE
```
Now you have created a virtual environment call myVE in the directory virtualenv-15.0.0.  Next you need to activate the environment.  In order to activate your virutal environment type in the commmand below:
```
$ source bin/activate
```
Upon successful activation, your commmand prompt should now read something like this:
```
(myVE)user@myaddress$  
```
where the parantheses is the name of the virtual environment you just created and appears to the far left to let you know you have activated and are working in the virtual environment.  

Now anytime you would like to work in the virutal environment make you use
```
$ source bin/activate
```
and see the parenthesis and it is all set up!  Note, however, although you every only need to install the virtual environment once, each time you log back into your system you will need to activate the environment if you want to use it.  In order to logout or deactivate the virtual environment just type the following into your command prompt:
```
$ deactivate
```
Now you are working in your regular environment.  All the data you generated in the virtual environment is saved and you can access it like any other directory or folder in your regular working environment.

It is also **critical** on HPC systems to load the python module **PRIOR TO ACTIVATION** of your virtual environment.  

**2.  Clone Repository (only needs to be done once)**
------------------------------------------------------

First make sure you have loaded the Python module on your HPC and then activated your virutal environment.  Next, within your virutual environment clone the this git repository.
```
$ module load <python module on HPC>
$ cd /path/to/myVE
$ source bin/activate
(myVE)$ git clone https://github.com/tbrunetti/GWAS_Analysis_Pipeline.git
(myVE)$ cd GWAS_Analysis_Pipeline
```
The pipeline has been successfully cloned into your virtual environment.  Next, we will need to install and configure the pipeline using chunkypipes.  

**3. Activate chunkypipes, configure and install pipeline (only needs to be performed once)**
----------------------------------------------------------------------------------------------
With your virutal environment still active, install chunkypipes into your virtual environment.
```
(myVE)$ pip install chunkypipes
(myVE)$ chunky init
```
If chunkypipes was sucessfully installed and initiated the follow message will appear:  
```
> ChunkyPipes successfully initialized at /home/user
```
Next, we will used chunkypipes to install the GWAS pipeline.
```
(myVE)$ chunky install run_GWAS_analysis_pipeline.py
```
If the pipeline has been successfull installed the following message will appear:
```
> Pipeline run_GWAS_analysis_pipeline.py successfully installed
```
chunkypipes will then ask the user if they would like to install all the pipeline dependencies.  If the user is confident that these programs have already been downloaded and are functional OR they are using module load on the HPC for all these dependencies, the user can decline by typing in 'n' followed by enter. Pressing any other letter followed by enter will tell the pipeline to automatically install all the dependencies which is strongly recommended.
```
Attempting to install the following dependencies:
pandas
numpy
matplotlib
fpdf
Pillow
pypdf2
statistics
xlrd

Proceed with dependency installation? [y/n] 
```
Lastly, we must configure the pipeline.  Note, the configuration only needs to be performed once unless the paths of the required software and files has changed.
```
(myVE)$ chunky configure run_GWAS_analysis_pipeline.py
```
This will prompt the user for the full path locations of the PLINK, KING, 1000 Genomes LD pruned data, 
```
Full path to KING executable []:
Full path to directory where R libraries are stored []:
Full path to PLINK executable []:
Full path to PLINK BED file format LD pruned phase 3 1000 genomes files []:
Configuration file successfully written
```  


**4.   Modifying the Slurm script**
------------------------------------
Now that the pipeline had been successfully installed and configured the user can modify the provided slurm script.  The slurm script is the only script that the user will have to submit in order to run the pipeline.  This script is called *run_analysis_HPC_slurm.sh*.  All other python scripts, shell scripts, and R scripts are taken care of on the backend.  

Any of the #SBATCH lines can be modified based on the specs and allocations of the HPC.  
```
#!/bin/bash

#SBATCH --time=1440
#SBATCH --ntasks=5
#SBATCH --mem=100000
#SBATCH --job-name=prototyping_pipeline
#SBATCH --output=prototyping_pipeline.log
#SBATCH -p defq
```
Although the script  has a time set of 24 hours (1440 min), when running the pipeline from the very beginning throught PCA calculations (the most time intensive part) runs ~1700 samples on the MEGA chip (~1.5 million SNPS per sample post initial QC cleanup) in about ~10 hours.  Although this time will shorten once the code is parallelized.  

It is recommended not to set the memory parameter too low since calculating the PCs and running GENESIS may require large amounts of memory.  To give the user an idea of how much to allocate, this pipeline was prototyped on the MEGA chip (~1.5 million snps per sample post initial QC cleanup) across ~1700 samples and it runs from the very beginning of the pipeline though PCA calculation step using less than 100GB of memory.  The same amount of memory allocation was used to run the GENanalysis step of the pipeline, although this is can probably be significantly scaled down.  
  
The eval statment can be commented out if the --account flag is not going to be used when submitting the script via sbatch.
```
eval "#SBATCH --account=TICR=${USER}"
```
The slurm script should also be modified to load the appropriate modules on the HPC:
```
module load <python location>
module load <R location>
module load <torque location>
module load <libpng location>
module load <freetype location>
module load <gcc location>
```
Replace <module location> with the location of the module on the users' HPC.  Some of the modules such as libpng and freetype may be automatically unloaded on some HPC systems since they are standard Linux libraries.  

Next the user must specify to use the Python virutal environment.  To do this modify the following lines in the slurm script:
```
cd <into your virtual environment>
source bin/activate
cd <into cloned repository>
```
Replace <into your virutal environment> with the full path into your virutal environment (it should be the location where you see the bin directory of the virutal environment) and replace <into cloned respository> with the full path into the clone git repository located within the virutal environment directory (i.e. change directories into GWAS_Analysis_Pipeline).  

Finally, the last part of the slurm file is the command to run the pipeline. To run the pipeline, the user must call 'chunky run' followed by the name of the pipeline, which is run_GWAS_analysis_pipeline.py
```
chunky run run_GWAS_analysis_pipeline.py -inputPLINK <plink file.bed or .ped> -phenoFile <sample_sheet_template.xlsx> --outDir <directory to output all results> --projectName <give run a name>
```
To minimally run this pipeline two parameters are required:  -inputPLINK and -phenoFile.  Replace <plink file.bed or .ped> with the full path to the PLINK files you want to analyze.  Replace <sample_sheet_template.xlsx> with the full path or renamed and populated sample_sheet_template.xlsx file.  Although --outDir and --projectName are not required it is **HIGHLY RECOMMENDED** that the user provides a pre-existing directory to redirect results to as it produces several files.  Additionally the projectName can be any string to identify the current run.  This will create a subdirectory within the output directory with the name of the run so all the results do not mix with other files or runs in the output directory.  By default, if this command is called as is, it will make sure the project name does not already exist and create a new subdirectory within the output directory and run the pipeline starting at Hardy-Weinberg equilibrium all the way through the PCA calculations.  **THE PIPELINE WILL NOT RUN THE ASSOCIATION ANALYSIS!!!**  The reason being, GENESIS requires that for association testing the user knows what PCs wil be used for correcting the data.  Since this requires the user to look at the data after the PCs are generated, the slurm script must be submitted again with the **EXACT SAME RUN COMMAND AS BEFORE** except 3 new parameters should be added: --usePCs, --reanalyze and --startStep GENanalysis need to be added as follows:
```
chunky run run_GWAS_analysis_pipeline.py -inputPLINK <plink file.bed or .ped> -phenoFile <sample_sheet_template.xlsx> --outDir <directory to output all results> --projectName <give run a name> --startStep GENanalysis --reanalyze  --usePCs pc1,pc2,pc3,pc4,pc5,pc6,pc7
```
This tells the pipeline to use the same run as before except to start at the GENESIS analysis step and go through the association testing analysis.  The reanalyze flag tells the pipeline that the project already exists and to not overwrite any data and not to exit the system.  Although the --usePCs parameter is optional, it is recommended that the user list the PCs that they would like corrected.  If the parameter is not used, it will default to correcting for the first 3 PCs in the set.  A maximum of the first 20 PCs can be used and corrected for.


This call can be modified with any of the optional parameters that the pipeline offers.  For more information and details on all the parameters please refer to the Pipeline Parameters section.



## Installing and Running Pipeline with on personal system with sudo privileges
-------------------------------------------------------------------------------
1. Clone Repository (only needs to be done once)
2. Activate chunkypipes, configure and install pipeline (only need to be performed once)
3. Run pipeline  

Note:  If you choose **NOT** to run this pipeline on a HPC with slurm, the pipeline cannot utilize the GENanalysis step.  This is because due to potentially large space and memory requirements for this part of the analysis, the power of the HPC is utilized to streamline this process.  However, the individual R scripts have been provided so that the user can download these and modify and run them on their own personal computer without the HPC dependency.  


## Pipeline Set-Up and Configuration Logic
-------------------------------------------
<p align="center">
<img src="https://github.com/tbrunetti/GWAS_Analysis_Pipeline/blob/master/how_to_run_flowchat.png" />
</p>

