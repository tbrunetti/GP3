from chunkypipes.components import *
import collections
import multiprocessing
import math
import datetime
import sys
import os
import shutil
import time

class Pipeline(BasePipeline):
	def dependencies(self):
		return ['pandas', 'numpy', 'matplotlib', 'fpdf', 'Pillow', 'pypdf2', 'statistics', 'xlrd']

	def description(self):
		return 'Pipeline made for analyzing GWAS data after QC cleanup'

	def configure(self):
		return {
			'plink':{
				'path':'Full path to PLINK executable'
			},
			'king':{
				'path':'Full path to KING executable'
			},
			'thousand_genomes':{
				'path': 'Full path to PLINK BED file format LD pruned phase 3 1000 genomes file'
			},
			'R_libraries':{
				'path': 'Full path to directory where R libraries are stored'
			},
			'R':{
				'path': 'Full path to R/Rscript executable'
			},
			'TGP_populations':{
				'path': 'Full path including file name of the TGP_Sub_and_SuperPopulation_info.txt file'
			}
		}

	def add_pipeline_args(self, parser):
		# should other input options be available??
		parser.add_argument('-inputPLINK', required=True, type=str, help='Full path to PLINK file ending in .BED or .PED, no whitespaces in file name!')
		parser.add_argument('-phenoFile', required=True, type=str, help='Full path to phenotype file, see argument readme for more details on format')
		parser.add_argument('--outDir', default=os.getcwd(), type=str, help='[default=current working directory] Full path of existing directory to output results, no spaces/whitespaces!')
		parser.add_argument('--projectName', default=str(datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")), type=str, help='[default=date/time stamp] Name of project, no spaces/whitespaces')
		parser.add_argument('--startStep', default='hwe', type=str, help='[default=hwe, options: if --TGP not set -> [hwe, LD, maf, het, ibd, PCA_indi, or PCA_indi_graph] if --TGP is set ->\
																		[hwe, LD, maf, het, ibd, PCA_TGP, or PCA_TGP_graph]] The part of the pipeline you would like to start with')
		parser.add_argument('--endStep', default=None, type=str, help='[default=PCA_indi_graph (noTGP) or PCA_TGP_graph (if --TGP used), options: if --TGP not set -> [hwe, LD, maf, het, ibd, PCA_indi, or PCA_indi_graph] if --TGP is set -> \
																		[hwe, LD, maf, het, ibd, PCA_TGP, or PCA_TGP_graph]] Point of the pipeline where you would like to \
																		stop analysis')
		parser.add_argument('--hweThresh', default=1e-6, help='[default=1e-6, options: any scientic notation value] Filters out SNPs that are smaller than this threshold due to liklihood of genotyping error')
		parser.add_argument('--LDmethod', default='indep', type=str, help='[default=indep, options:indep, indep-pairwise or indep-pairphase] Method to calculate linkage disequilibrium')
		parser.add_argument('--VIF', default=2, type=int, help='[default=2, options: any INT] variant inflation factor for indep method LD pruning')
		parser.add_argument('--rsq', default=0.50, type=float, help='[default=0.50, options: any FLOAT between 0.0-1.0] r squared threshold for indep-pairwise or indep-pairphase LD pruning method')
		parser.add_argument('--windowSize', default=50, type=int, help='[default=50, options: any INT] the window size in kb for LD analysis')
		parser.add_argument('--stepSize', default=5, type=int, help='[default=5, options: any INT] variant count to shift window after each interation')
		parser.add_argument('--maf', default=0.05, type=float, help='[default=0.05, option: any FLOAT between 0.0-1.0], filter remaining LD pruned variants by MAF, any MAF below set threshold is filtered out')
		parser.add_argument('--hetMethod', default='meanStd', type=str, help='[default=meanStd, options:  minMax or meanStd], method to use to determine heterozygosity.  minMax filter based on \
																			the parameters --hetThresh as the max F-inbreeding coefficient and --hetThreshMin for the minimum F-inbreeding coeffient, which \
																			by default are 0.10 and -0.10, respectively.  The meanStd filter method calculates a het_score: \
																			1-[observed[HOM]/total] and then filters out any samples that are more than 3 std deviations from the mean het_score. \
																			The number of standard deviations from the mean can be changed using the --het_std parameter')
		parser.add_argument('--het_std', default=3, type=float, help='[default=3, options: any FlOAT or INT] if using hetMethod=meanStd you can determine how many standard deviations aways from the mean is allowable for heterozygosity.  \
																	Setting to 3 is interpreted as +/-3 standard deviations away from the mean of the het_score, calculated as 1-[observed(HOM)/total]')
		parser.add_argument('--hetThresh', default=0.10, type=float, help='[default=0.10, options: any FlOAT], filter out samples where inbreeding coefficient is greater than threshold (heterozygosity filtering)')
		parser.add_argument('--hetThreshMin', default=-0.10, type=float, help='[default= -0.10, options: any FlOAT] filter out samples where inbreeding coefficient is samller than min threshold set (heterozygosity filtering)')
		parser.add_argument('--reanalyze', action='store_true', help='by adding this flag, it means you are going to pass a dataset through the pipeline that has already been partially/fully analyzed by this pipeline. WARNING! May over write exisiting data!!')
		parser.add_argument('--sampleMiss', default=0.03, type=float, help='[default: 0.03, options: any FLOAT between 0.0-1.0] Maximum missingness of genotype call in sample before it should be filtered out.  Float between 0-1, where 0 is no missing, and 1 is all missing (0.03 is interpreted as 3 percent of calls are missing)')
		parser.add_argument('--snpMiss', default=0.03, type=float, help='[default: 0.03], options: any FLOAT between 0.0-1.0] Maximum missingness of genotype call in a SNP cluster before it should be filtered out.  Float between 0-1, where 0 is no missing, and 1 is all missing (0.03 is interpreted as 3 percent of calls are missing)')
		parser.add_argument('--TGP', action='store_true', help='specifying this flag means to generate PCA plots with TGP data merged into the given cohort data set for the 5 superpopulations in TGP (AFR, AMR, EAS, EUR, SAS)')
		parser.add_argument('--centerPop', default='myGroup', type=str, help="[default: myGroup, options: myGroup or available TGP group merged into input dataset] when using the TGP flag, you have the option to specify which population cohort that PCs should be centered around for boxplots.  By default this is set to your group(s).  You can pick a TGP super population listed in the TGP_Sub_and_SuperPopulation_info.txt file. CASE SENSITIVE!")
		parser.add_argument('--outliers', default=None, type=str, help="A txt file of FID and IID, tab-delimited and one sample per line, that are outliers that should be removed from the sample set (PCA outlier removal); Use original names (original FID and IID), not renamed 1-n for GENESIS formatting")
		parser.add_argument('--pcmat', default=5, type=int, help='[default=5, INT] Number of predicted admizture populations in dataset to be used in GENESIS calculation for PCA')
	
	
	@staticmethod
	def check_steps(order, start, stop):
		if start == 'hwe' and stop == None:
			return order
		else: # reanalyze flag should be used here
			index_of_start =order.index(start)
			if stop == None:
				return order[index_of_start:]
			else:
				index_of_end = order.index(stop)
				return order[index_of_start:index_of_end+1]


	# checks the type of PLINK file input
	@staticmethod
	def check_plink_format(plinkFile, plink):
		if plinkFile[-4:].lower() == '.ped':
			print "Input .ped, converting to binary"
			convert = subprocess.call([str(plink), '--file', str(plinkFile[:-4]), '--make-bed', '--out', str(plinkFile[:-4])])
		
		elif plinkFile[-4:].lower() == '.bed':
			print "Input seems to follow bed format"

		else:
			sys.exit("Error!! Input not recognized, please input .ped or .bed PLINK file" )


	# create files to input into plink for samples to keep
	@staticmethod
	def ethnic_plinks_lists(phenotype, plinkFileName, famFile, outDir):
		import pandas as pd

		phenotype_table = pd.read_excel(phenotype, sheetname="Sheet2", header=0, converters={'FID':str, 'IID':str}) 
		phenotype_table['IID'] = phenotype_table['IID'].str.strip() # remove whitespace for proper merging
		phenotype_table['FID'] = phenotype_table['FID'].str.strip() # remove whitespace for proper merging
		
		# update fam file with phenotype information
		original_fam = pd.read_table(famFile, delim_whitespace=True, names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'], converters={'FID':str, 'IID':str, 'Gender':str}) 
		original_fam['IID'] = original_fam['IID'].str.strip() # remove whitespace for proper merging
		original_fam['FID'] = original_fam['FID'].str.strip() # remove whitespace for proper merging

		merged_dataframe = original_fam.merge(phenotype_table, how='left', on=['FID', 'IID']) # only merge information where FID and IID are available in the original fam file

		race_subsets = list(set(list(merged_dataframe['Race']))) # gets all racial groups listed in phenotype file
		race_subsets_cleaned = [i for i in race_subsets if str(i)!='nan'] # removes samples that do not have a race listed 
		make_plinks = {} # stores file name of samples to keep for each group
		

		for ethnic_group in race_subsets_cleaned:
			os.mkdir(outDir + '/' + '_'.join(ethnic_group.split()))
			keep_file = open(outDir + '/' + '_'.join(ethnic_group.split()) + '/' + plinkFileName +'_'+str('_'.join(ethnic_group.split())) + '_keepIDs.txt', 'w')
			
			# need to check for last part of conditional in case FID is missing, remove the row...causes problems with PLINK
			subset_by_race = merged_dataframe.loc[(merged_dataframe['Race'] == ethnic_group) & (pd.isnull(merged_dataframe['FID'].str.strip()) == False)]

			subset_by_race[['FID', 'IID']].to_csv(keep_file.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>
			merged_dataframe[['FID', 'IID', 'PAT', 'MAT', 'Gender', 'Phenotype']].to_csv(famFile, sep=' ', index=False, header=False)
			make_plinks['_'.join(ethnic_group.split())]=keep_file.name
			keep_file.flush()
			keep_file.close()


		return make_plinks


	def run_pipeline(self, pipeline_args, pipeline_config):
		import pandas as pd
		import numpy as np
		sys.path.append(".")
		import subprocess
		import summary_stats
		import statistics as stats
		from fpdf import FPDF
		import PyPDF2

		# keeps track of all files that can be deleted after the pipeline finishes
		stage_for_deletion = []

		reduced_plink_name = pipeline_args['inputPLINK'].split('/')[-1][:-4] #only get plink file name not full absolute path and removes suffix
		
		# specifying output location and conflicting project names of files generated	
		try:
			if pipeline_args['reanalyze'] == True:
				outdir = pipeline_args['outDir']+'/'+pipeline_args['projectName']
				print "Reanalyzing data from an existing project"
			else:
				os.stat(pipeline_args['outDir']+'/'+pipeline_args['projectName'])
				sys.exit("project already exists!!")
		except:
			print "Making new directory called "+str(pipeline_args['projectName']) + ' located in ' + str(pipeline_args['outDir'])
			outdir = pipeline_args['outDir']+'/'+pipeline_args['projectName'] # new output directory
			os.mkdir(outdir)
			#self.settings.logger.set(
			#	destination=pipeline_args['outDir'] + '/' + pipeline_args['projectName'] +'/stdout.log',
			#	destination_stderr=pipeline_args['outDir'] + '/' + pipeline_args['projectName'] + '/stderr.log')
 
		
		# determine what steps need to be performed
		# NOTE: outlier_removal is an optional step where if the --outlier argument is not specified, the step is skipped
		if pipeline_args['TGP'] == True:
			step_order = self.check_steps(
				order = ['hwe', 'LD', 'maf', 'het', 'ibd', 'outlier_removal', 'PCA_TGP', 'PCA_TGP_graph'], 
				start = pipeline_args['startStep'],
				stop = pipeline_args['endStep']
				)
		else:
			step_order = self.check_steps(
				order = ['hwe', 'LD', 'maf', 'het', 'ibd', 'outlier_removal', 'PCA_indi', 'PCA_indi_graph'], 
				start = pipeline_args['startStep'],
				stop = pipeline_args['endStep']
				)			


		# initialize PLINK and KING software
		general_plink = Software('plink', pipeline_config['plink']['path'])
		general_king = Software('king', pipeline_config['king']['path'])
		# make sure plink input format is correct
		self.check_plink_format(
			plinkFile = pipeline_args['inputPLINK'],
			plink = pipeline_config['plink']['path']
			)

		
		if pipeline_args['reanalyze'] == False:
			keep_files = self.ethnic_plinks_lists(
				phenotype = pipeline_args['phenoFile'],
				plinkFileName = reduced_plink_name,
				famFile = pipeline_args['inputPLINK'][:-4] + '.fam',
				outDir = outdir
				)


			# will make separate plink files for each ethnic group and use autosomes only
			for key, value in keep_files.iteritems():
				general_plink.run(
					Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
					Parameter('--keep', value),
					Parameter('--autosome'),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/' + str(key) + '/' + reduced_plink_name + '_' + str(key))
					)


		while len(step_order) != 0:
			
			# hardy-weinberg equilibrium filtering
			if step_order[0] == 'hwe':
				print "running HWE step"
				print multiprocessing.cpu_count()
				hwe_passing = {}
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories),
							Parameter('--geno', pipeline_args['snpMiss']),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name + '_snpMissFiltered_' + directories)
							)

						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_snpMissFiltered_' + directories),
							Parameter('--hardy'),
							Parameter('--hwe', pipeline_args['hweThresh']),
							Parameter('midp'),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_hweFiltered')
							)
						
						total_snps_analyzed_hwe = subprocess.check_output(['wc', '-l',  outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.bim'])
						total_snps_passing_hwe = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_hweFiltered.bim'])
						hwe_passing[directories] = [total_snps_analyzed_hwe.split()[0]] + [total_snps_passing_hwe.split()[0]] # store total analyzed and passing for hwe step
						hwe_passing[directories] = hwe_passing[directories] + [outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_hweFiltered.hwe']
				
				hwe_stats = summary_stats.hwe(dictHWE=hwe_passing, thresh=pipeline_args['hweThresh'], outDir = outdir)
				hwe_stats.output(outdir + '/' + pipeline_args['projectName'] + '_hweStats.pdf', 'F') # output results to PDF
				
				for directories in os.listdir(outdir): # delete individual hwe images 
					if (os.path.isdir(os.path.join(outdir, directories))):
						subprocess.call(['rm', '-rf', 'hwe_'+str(directories)+'.png'])
				
				step_order.pop(0)

			
			# LD pruning
			elif step_order[0] == 'LD':
				print "running LD pruning step"
				ld_passing = {}

				if pipeline_args['LDmethod'] == 'indep': # gets lists of variants to keep and to prune using VIP
					for directories in os.listdir(outdir):
						if (os.path.isdir(os.path.join(outdir, directories))):
							samples_missing = open(outdir + '/'+ directories + '/' + reduced_plink_name + '_' + directories + '_failedSampleMiss.txt', 'w') # file to record samples failing gentype call threshold
							# below plink call removes samples that fail overall genotyping call threshold
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
								Parameter('--mind', pipeline_args['sampleMiss']),
								Parameter('--make-bed'),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered_sampleMissFiltered')
								)
							# below plink call returns file of the missingness calls PRIOR to removing samples failing genotyping call threshold
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
								Parameter('--missing'),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered')
								)
							# below plink call for LD pruning
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered_sampleMissFiltered'),
								Parameter('--'+pipeline_args['LDmethod'], str(pipeline_args['windowSize'])+'kb', str(pipeline_args['stepSize']), str(pipeline_args['VIF'])),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories)
								)
							
							missingness_table = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered.imiss', delim_whitespace=True)
							get_missing = missingness_table.loc[missingness_table['F_MISS'].astype(float) > pipeline_args['sampleMiss']]
							get_missing[['FID', 'IID']].to_csv(samples_missing.name, sep='\t', index=False, header=False)
							samples_missing.close()

							total_snps_analyzed_ld = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered_sampleMissFiltered.bim'])
							total_snps_passing_ld = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.in'])
							ld_passing[directories] = [total_snps_analyzed_ld.split()[0]] + [total_snps_passing_ld.split()[0]] # store total analyzed and passing for LD pruning step
							
							# creates new PLINK files with excluded variants removed
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered_sampleMissFiltered'),
								Parameter('--exclude', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.out'),
								Parameter('--make-bed'),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned')
								)
							
				
				else:
					for directories in os.listdir(outdir): # get lists of variants to keep and to prune using rsq
						if (os.path.isdir(os.path.join(outdir, directories))):
							samples_missing = open(outdir + '/'+ directories + '/' + reduced_plink_name + '_' + directories + '_failedSampleMiss.txt', 'w') # file to record samples failing gentype call threshold
							# below plink call removes samples that fail overall genotyping call threshold
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'), 
								Parameter('--mind', pipeline_args['sampleMiss']),
								Parameter('--make-bed'),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered_sampleMissFiltered')
								)
							# below plink call returns file of the missingness calls PRIOR to removing samples failing genotyping call threshold
							general_plink.run(
                            	Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
                            	Parameter('--missing'),
                            	Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered')
                            	)
							# below plink call for LD pruning
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
								Parameter(pipeline_args['LDmethod'], str(pipeline_args['windowSize'])+'kb', str(pipeline_args['stepSize']), str(pipeline_args['rsq'])),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories)
								)

							missingness_table = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered.imiss', delim_whitespace=True)
							get_missing = missingness_table.loc[missingness_table['F_MISS'].astype(float) > pipeline_args['sampleMiss']]
							get_missing[['FID', 'IID']].to_csv(samples_missing.name, sep='\t', index=False, header=False)
							samples_missing.close()

							total_snps_analyzed_ld = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered_sampleMissFiltered.bim'])
							total_snps_passing_ld = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.in'])
							ld_passing[directories] = [total_snps_analyzed_ld] + [total_snps_passing_ld] # store total analyzed and passing for LD pruning step
							
							# creates new PLINK files with excluded variants removed
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered_sampleMissFiltered'),
								Parameter('--exclude', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.out'),
								Parameter('--make-bed'),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned')
								)
					
				ld_stats = summary_stats.pruning(dictLD=ld_passing)
				ld_stats.output(outdir + '/' + pipeline_args['projectName'] + '_ldStats.pdf', 'F') # output results to PDF
				step_order.pop(0)

			# filters pruned variants by MAF
			elif step_order[0] == 'maf':
				print "running maf step"
				list_of_merge_maf_greater_thresh = open(outdir + '/' + reduced_plink_name + '_greater_mafs.txt', 'w')
				list_of_merge_maf_less_thresh = open(outdir + '/' + reduced_plink_name + '_small_mafs.txt', 'w')
				
				maf_passing = {}
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						
						# filters out variants below set maf threshold
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned'),
							Parameter('--maf', str(pipeline_args['maf'])),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf')
							)

						total_snps_analyzed_maf = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned.bim'])
						total_snps_greater_maf = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.bim'])
						maf_passing[directories] = [total_snps_analyzed_maf.split()[0]] + [total_snps_greater_maf.split()[0]]

						# stores name of file for future merging (must be space delimited ordered bed, bim, fam files)
						list_of_merge_maf_greater_thresh.write(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.bed ' + 
							outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.bim ' +
							outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.fam' + '\n')

						# filters out variants set maf threshold
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned'),
							Parameter('--max-maf', str(pipeline_args['maf'])),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf')
							)
						
						total_snps_less_maf = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf.bim'])
						maf_passing[directories] = maf_passing[directories] + [total_snps_less_maf.split()[0]]
	

						# stores name of file for future merging
						list_of_merge_maf_less_thresh.write(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf.bed ' +
							outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf.bim ' +
							outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf.fam' +'\n')
					
				
						# recalculate freq stats for each set (het and ibd are sensitive to maf, so make sure freq is recalculated)
						general_plink.run(
							Parameter('--bfile', (outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf')),
							Parameter('--freq'),
							Parameter('--out', (outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + 'freq_greater_than_'+str(pipeline_args['maf'])+'_maf'))
							)	
				
				list_of_merge_maf_greater_thresh.flush() # push out buffer
				list_of_merge_maf_less_thresh.flush() # push out buffer
				
				maf_stats = summary_stats.minor_allele_freq(dictMAF=maf_passing, thresh=pipeline_args['maf'])
				maf_stats.output(outdir + '/' + pipeline_args['projectName'] + '_mafStats.pdf', 'F') # output results to PDF

				step_order.pop(0)
				
		

			elif step_order[0] == 'het':
				print "checking heterozygosity"
				het_pdf_merge = PyPDF2.PdfFileMerger() # PDF merge object
				clean_up = []
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):

						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf'),
							Parameter('--read-freq', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + 'freq_greater_than_'+str(pipeline_args['maf'])+'_maf.frq'),
							Parameter('--het'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf')
							)


						het_dataframe = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.het', delim_whitespace=True)
						# het_score calculates the heterozygosity score by the following calculation: 1-[oberserved(HOM)/total]
						het_dataframe['het_score'] = 1-(het_dataframe['O(HOM)'].astype(float)/het_dataframe['N(NM)'].astype(float))
						samples_failing_het, het_pdf = summary_stats.heterozygosity(het_method= pipeline_args['hetMethod'], std=pipeline_args['het_std'], het_dataframe = het_dataframe, thresh = pipeline_args['hetThresh'], minThresh=pipeline_args['hetThreshMin'], population=directories, outDir = outdir)
						
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf'),
							Parameter('--remove', samples_failing_het),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered')
							)
					
						
						het_pdf.output(outdir + '/' + pipeline_args['projectName'] + '_' + directories + '_hetStats.pdf', 'F') # output results to PDF
						clean_up.append(outdir + '/' + pipeline_args['projectName'] + '_' + directories + '_hetStats.pdf')
						het_pdf_merge.append(outdir + '/' + pipeline_args['projectName'] + '_' + directories + '_hetStats.pdf')
				
				het_pdf_merge.write(outdir + '/' + pipeline_args['projectName'] + '_hetStats.pdf')
				het_pdf_merge.close()

				for delete_pdfs in clean_up: # remove individual pdfs after merging with final pdf
					subprocess.call(['rm', '-rf', delete_pdfs])
				
				step_order.pop(0)
				


			# determine the pairwise relationship of all merged samples
			# DUPLICATES MUST BE REMOVED BEFORE USING GENESIS PIPELINE!
			elif step_order[0] == 'ibd':
				# only gets run on the bed file with MAF > specified tresh
				print "running PLINK ibd step"
				ibd_pdf_merge = PyPDF2.PdfFileMerger() # PDF merge object
				clean_up = []
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):

						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered'),
							Parameter('--genome'),
							Parameter('--read-freq', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + 'freq_greater_than_'+str(pipeline_args['maf'])+'_maf.frq'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered')
							)
						
						ibd_results = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered.genome', delim_whitespace=True)
						relatedness_stats, remove_samples = summary_stats.relatedness(ibd_dataframe = ibd_results, population=directories, outDir=outdir+'/'+directories)

						# remove samples that are duplicates
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered'),
							Parameter('--remove', remove_samples),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed')
							)
						
						relatedness_stats.output(outdir + '/' + pipeline_args['projectName'] + '_' + directories + '_relatedness.pdf', 'F') # output results to PDF
						clean_up.append(outdir + '/' + pipeline_args['projectName'] + '_' + directories + '_relatedness.pdf')
						ibd_pdf_merge.append(outdir + '/' + pipeline_args['projectName'] + '_' + directories + '_relatedness.pdf')
				
				ibd_pdf_merge.write(outdir + '/' + pipeline_args['projectName'] + '_relatedness.pdf') # merge all pds together
				ibd_pdf_merge.close()

				for delete_pdfs in clean_up: # deletes individual pdfs already concatenated together in final ibd pdf
					subprocess.call(['rm', '-rf', delete_pdfs])

				step_order.pop(0)
	
			# remove outliers based input file list
			# if no list is provided (i.e outlier option is not specified) then just creates copy of file renamed
			# user will notice that the file outliers_removed_directory.txt is empty if no list is provided since no outliers present
			elif step_order[0] == 'outlier_removal':
				if pipeline_args['outliers'] == None:
					for directories in os.listdir(outdir):
						if (os.path.isdir(os.path.join(outdir, directories))):
							print "Skipping outlier removal step"
							outlier_list = open(outdir + '/' + directories + '/' + 'outliers_removed_' + directories + '.txt', 'w')
							# copy file with new name, suffix outliers_removed
							shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed.bed',
											outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.bed')
							shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed.bim',
											outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.bim')
							shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed.fam',
											outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.fam')
							outlier_list.close()
							
					step_order.pop(0)
				else:
					for directories in os.listdir(outdir):
						if (os.path.isdir(os.path.join(outdir, directories))):
							print "removing outliers post IBD step as specified by user input"
							outlier_list = open(outdir + '/' + directories + '/' + 'outliers_removed_' + directories + '.txt', 'w')
							input_list = pd.read_table(pipeline_args['outliers'], names = ['FID', 'IID'], dtype=str)
							fam_file = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed.fam', 
										names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'], dtype=str, delim_whitespace=True)

							ids_found = set(list(input_list[['FID', 'IID']].itertuples(index=False, name=None))).intersection(set(list(fam_file[['FID', 'IID']].itertuples(index=False, name=None))))

							for samples in ids_found:
								outlier_list.write(samples[0] + '\t' + samples[1] + '\n')

							outlier_list.flush()
							outlier_list.close()

							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed'),
								Parameter('--remove', outlier_list.name),
								Parameter('--make-bed'),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed')
								)

					step_order.pop(0)

			# PCA with 1000 genomes
			elif step_order[0] =='PCA_TGP':

				print "merging with 1000_genomes"
				phenoFiles = {}
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						keep_these_snps_in_1000 = open(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_all_passing_snps_from_project.txt', 'w')
						bim_file = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.bim', names=['chr', 'SNP_name', 'pos1', 'pos2', 'allele1', 'allele2'])
						bim_file_1000 = pd.read_table(pipeline_config['thousand_genomes']['path'][:-4]+'.bim',  names=['chr', 'SNP_name', 'pos1', 'pos2', 'allele1', 'allele2'])
						extract_for_PCA = list(set(list(bim_file_1000['SNP_name'])) & set(list(bim_file['SNP_name'])))

						keep_these_snps_in_1000.write('\n'.join(list(extract_for_PCA))) # input file for PLINK extraction of snps in 1000 genomes
						keep_these_snps_in_1000.close() # flushes and closes file

						no_suffix = pipeline_config['thousand_genomes']['path'][:-4]

						general_plink.run(
							Parameter('--bfile', no_suffix),
							Parameter('--extract', keep_these_snps_in_1000.name),
							Parameter('--make-bed'),
							Parameter('--out', no_suffix + '_extracted_from_passing_' + directories +'_SNPs')
							)
						

						general_plink.run(
            				Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed'),
                    		Parameter('--extract', keep_these_snps_in_1000.name),
                			Parameter('--make-bed'),
                			Parameter('--out', outdir + '/' + directories + '/extracted_from_passing_' + directories +'_SNPs_for_use_with_1000_merge')
                            )


						general_plink.run(
							Parameter('--bfile',  outdir + '/' + directories + '/extracted_from_passing_' + directories +'_SNPs_for_use_with_1000_merge'),
							Parameter('--bmerge',no_suffix + '_extracted_from_passing_' + directories +'_SNPs.bed', no_suffix + '_extracted_from_passing_' + directories +'_SNPs.bim', no_suffix + '_extracted_from_passing_' + directories +'_SNPs.fam'),
							Parameter('--allow-no-sex'),
							Parameter('--make-bed'),
							Parameter('--out',  outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +'_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen')
							)
				

				
				print "running KING step"
				phenoFiles = {}
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						all_samples_to_remove = open(outdir + '/' + directories + '/' + directories + '_all_samples_to_remove_from_original_TGP.txt', 'w')
						phenoFile_Genesis = open(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen_phenoGENESIS.txt', 'w')
						phenoFiles[directories] = phenoFile_Genesis.name
						# run KING and output file as -b prefix name ending in .kin, .kin0 for each group
						general_king.run(
							Parameter('-b', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen.bed'),
							Parameter('--prefix', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen')
							)
						
						# generate phenotype table for input into GENESIS setup analysis pipeline WITHOUT 1000 genomes
						pheno_Genesis = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen.fam', delim_whitespace=True, names = ['FID,', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'])
						pheno_Genesis[['IID', 'AFF']].to_csv(phenoFile_Genesis.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>
						phenoFile_Genesis.close()

						original_fam_file = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '.fam', names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'], dtype=str, delim_whitespace=True)
						original_fam_file['tuple_name'] = list(zip(original_fam_file.FID, original_fam_file.IID))
						pheno_Genesis['tuple_name'] = list(zip(pheno_Genesis.FID, pheno_Genesis.IID))

						pre_filter = set(list(original_fam_file['tuple_name']))
						post_filter = set(list(pheno_Genesis['tuple_name']))
						remove_samples_list = list(pre_filter.difference(post_filter))
						for fails in remove_samples_list:
							all_samples_to_remove.write(str(fails[0]) + '\t' + str(fails[1]) + '\n')

						all_samples_to_remove.flush()
						all_samples_to_remove.close()


				print "running PCA step"
				
				processes = []
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						#Popen should launch jobs in parallel
						processes.append(subprocess.Popen([os.path.join(pipeline_config['R']['path'], 'Rscript'), 'GENESIS_setup_ANALYSIS_PIPELINE.R', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen', phenoFiles[directories], pipeline_config['R_libraries']['path'], str(pipeline_args['pcmat'])]))

				for job in processes:
					job.wait() # wait for all parallel jobs to finish before proceeding to next step


				# copies final output *_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed* in another file renamed to ease of searching final file name
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen.bed', outdir + '/' + directories + '/' + directories + '_all_steps_completed_TGP_final.bed')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen.bim', outdir + '/' + directories + '/' + directories + '_all_steps_completed_TGP_final.bim')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen.fam', outdir + '/' + directories + '/' + directories + '_all_steps_completed_TGP_final.fam')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen.kin', outdir + '/' + directories + '/' + directories + '_all_steps_completed_TGP_final.kin')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen.kin0', outdir + '/' + directories + '/' + directories + '_all_steps_completed_TGP_final.kin0')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen.gds', outdir + '/' + directories + '/' + directories + '_all_steps_completed_TGP_final.gds')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen_GENESIS', outdir + '/' + directories + '/' + directories + '_all_steps_completed_TGP_final_GENESIS.Rdata')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen_phenoGENESIS.txt', outdir + '/' + directories + '/' + directories + '_all_steps_completed_TGP_final_phenoGENESIS.txt')
				

				step_order.pop(0)

			elif step_order[0] == "PCA_TGP_graph":
				
				print "graphing PCA plots"

				graphing_processes = []
				for directories in os.listdir(outdir):					
					if (os.path.isdir(os.path.join(outdir, directories))):
						if pipeline_args['centerPop'] == 'myGroup':
							center_pop = str(directories)
						else:
							center_pop = pipeline_args['centerPop']
						#Popen should launch jobs in parallel for producing PDFs of graphs for each cohort
						graphing_processes.append(subprocess.Popen([os.path.join(pipeline_config['R']['path'], 'Rscript'), 'PCA_TGP.R', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_thousGen_GENESIS', str(directories), pipeline_config['R_libraries']['path'], outdir + '/' + directories + '/',
																		pipeline_config['TGP_populations']['path'], center_pop]))

				for pdfs in graphing_processes:
					pdfs.wait() # wait for all parallel jobs to finish

				step_order.pop(0)

	
			# PCA without TGP, only on individual cohort level
			# prepares GENESIS R object output for direct input into GENESIS association pipeline and renames sample IDs to be compatible with GENESIS association pipeline
			# performs the KING, GENESIS setup and PC calculation, and generates PCA graphs and screen plots for each cohort
			elif step_order[0] =='PCA_indi':

				print "running KING step"
				phenoFiles = {}
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						all_samples_to_remove = open(outdir + '/' + directories + '/' + directories + '_all_samples_to_remove_from_original.txt', 'w')

						fam_key = open(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_phenoGENESIS_number_ids_included.txt', 'w')
						phenoFile_Genesis = open(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_phenoGENESIS.txt', 'w')
						phenoFiles[directories] = phenoFile_Genesis.name
						# generate phenotype table for input into GENESIS setup analysis pipeline WITHOUT 1000 genomes
                                                pheno_Genesis = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.fam', delim_whitespace=True, names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'], dtype=str)

						pheno_Genesis.insert(6, 'numberID', range(1, len(pheno_Genesis)+1))
						pheno_Genesis.loc[(pheno_Genesis['AFF']!='1') & (pheno_Genesis['AFF']!='2'), 'AFF']='NA'
						pheno_Genesis[['numberID', 'AFF']].to_csv(phenoFile_Genesis.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>
						phenoFile_Genesis.close()
						pheno_Genesis[['FID', 'numberID', 'PAT', 'MAT', 'SEX', 'AFF']].to_csv(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.fam', sep='\t', index=False, header=False)

						pheno_Genesis[['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF', 'numberID']].to_csv(fam_key.name, sep='\t', index=False)
						# run KING and output file as -b prefix name ending in .kin, .kin0 for each group
						general_king.run(
							Parameter('-b', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.bed'),
							Parameter('--prefix', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed')
							)


						original_fam_file = pd.read_table(outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '.fam', names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'], dtype=str, delim_whitespace=True)
						original_fam_file['tuple_name'] = list(zip(original_fam_file.FID, original_fam_file.IID))
						pheno_Genesis['tuple_name'] = list(zip(pheno_Genesis.FID, pheno_Genesis.IID))

						pre_filter = set(list(original_fam_file['tuple_name']))
						post_filter = set(list(pheno_Genesis['tuple_name']))
						remove_samples_list = list(pre_filter.difference(post_filter))
						for fails in remove_samples_list:
							all_samples_to_remove.write(str(fails[0]) + '\t' + str(fails[1]) + '\n')

						all_samples_to_remove.flush()
						all_samples_to_remove.close()

				print "running PCA step"
				
				processes = []
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						#Popen should launch jobs in parallel
						processes.append(subprocess.Popen([os.path.join(pipeline_config['R']['path'], 'Rscript'), 'GENESIS_setup_ANALYSIS_PIPELINE.R', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed', phenoFiles[directories], pipeline_config['R_libraries']['path'], str(pipeline_args['pcmat'])]))

				for job in processes:
					job.wait() # wait for all parallel jobs to finish before proceeding to next step

				
				# copies final output *_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed* in another file renamed to ease of searching final file name
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.bed', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final.bed')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.bim', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final.bim')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.fam', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final.fam')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.kin0', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final.kin0')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.gds', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final.gds')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_GENESIS', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final_GENESIS.Rdata')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_phenoGENESIS_number_ids_included.txt', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final_GENESIS_sample_key_file.txt')
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_phenoGENESIS.txt', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final_phenoGENESIS.txt')
					try:
						shutil.copyfile(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed.kin', outdir + '/' + directories + '/' + directories + '_all_steps_completed_final.kin')
					except IOError:
						print("WARNING: " + outdir + '/' + directories + '/' + directories + '_all_steps_completed_final.kin')

				
				step_order.pop(0)
		
			elif step_order[0] == "PCA_indi_graph":

				print "graphing PCA plots"

				graphing_processes = []
				for directories in os.listdir(outdir):
					if (os.path.isdir(os.path.join(outdir, directories))):
						#Popen should launch jobs in parallel for producing PDFs of graphs for each cohort
						graphing_processes.append(subprocess.Popen([os.path.join(pipeline_config['R']['path'], 'Rscript'), 'PCA_indi.R', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories +  '_maf_greater_thresh_hetFiltered_dups_removed_outliers_removed_GENESIS', str(directories), pipeline_config['R_libraries']['path'], outdir + '/' + directories + '/']))

				for pdfs in graphing_processes:
					pdfs.wait() # wait for all parallel jobs to finish

				step_order.pop(0)

		print "writing results to PDF"
		paramsThresh = summary_stats.parameters_and_thresholds(params=pipeline_args)


		# output PDFs -- need to make this compatible with --reanalyze
		# put these under each of the steps
		paramsThresh.output(outdir + '/' + pipeline_args['projectName'] + '_parameters_and_thresholds.pdf', 'F') # output results to PDF

