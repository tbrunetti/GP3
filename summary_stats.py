import pandas as pd
import statistics as stats
import matplotlib as mpl
mpl.use('Agg') # bypass X11 if using server where X11 is not supported
import matplotlib.pyplot as plt
import numpy as np

from fpdf import FPDF

def parameters_and_thresholds(params):
	
	pdf = FPDF()
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Parameters and Thresholds", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', '', 16)
	for key in params:
		if key not in ['inputPLINK', 'phenoFile', 'outDir', 'projectName', 'config']:
			pdf.multi_cell(0, 8, str(key)+':     '+str(params[key]), 0, 1, 'L')

	return pdf



def hwe(dictHWE, thresh, outDir):
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Hardy-Weinberg Equilibrium", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_x(20)
	pdf.set_font('Arial', '', 12)
	pdf.multi_cell(0, 5, 'Hardy-Weinberg equilibrium is only used to remove SNPs with extreme p-values are that are likely \
		to occur due to sequencing, genotyping, or study-design errors.  This calculation is sensitive to different ethinic groups \
		and races.  Therefore, it is independently calculated for each ethnic group.  The current p-value threshold that was used to determine \
		whether a SNP was removed was  ' + str(thresh) + '.  This calculation will only consider founders; nonfounders are ignored.' , 0, 1, 'J')
	pdf.multi_cell(0, 5, '\n', 0, 1, 'J')
	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)
	
	# iterate through all ethnic groups for HWE stats
	for key, value in dictHWE.iteritems():
		pdf.multi_cell(0, 8, str(key), 1, 'L', True)
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs analyzed:  ' +  str(value[0]), 1, 1, 'L')
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs Passing:  ' +  str(value[1]) + ' (' + str("%.2f" % round((float(value[1])/float(value[0]))*100, 2)) + '%)', 1, 1, 'L')
		pdf.multi_cell(0, 8, '\n\n', 0, 1, 'J')

	
	# NOTE hweFile is before filtering by HWE threshold and prior to removal of SNPs failing threshold
	# these plot the observed vs expected from pre-filter HWE and the associated p-values
	# red fail threhold
	for key, value in dictHWE.iteritems():
		hweFile_dataframe = pd.read_table(value[2], delim_whitespace=True)	
		figure = plt.figure(1)
		num_phenos = len(list(set(list(hweFile_dataframe['TEST'])))) # for use in automating number of subplots to include in figure
		for phenotypes in list(set(list(hweFile_dataframe['TEST']))):
			pheno_subset = hweFile_dataframe.loc[hweFile_dataframe['TEST'] == phenotypes]
			colors = np.where(pheno_subset.P < thresh, 'r', 'k')
			plt.subplot(220 + num_phenos)
			plt.scatter(pheno_subset['E(HET)'], pheno_subset['O(HET)'], c=colors, s=8)
			plt.xlabel('expected(het)', fontsize=8)
			plt.ylabel('observed(het)', fontsize=8)
			plt.title(phenotypes + ':  observed(het) vs expected(het) of HWE', fontsize=8)
			num_phenos = num_phenos - 1
		plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
		plt.savefig(outDir+'/'+'hwe_'+str(key)+'.png')
		plt.close()
		pdf.add_page()
		pdf.set_margins(20, 10, 20)
		pdf.set_font('Arial', 'B', 14)
		pdf.set_x(20)
		pdf.multi_cell(0, 10, "HWE Plots for "+str(key) +" population", 0, 1, 'L')
		pdf.line(20, 32, 190, 32)
		pdf.set_x(20)
		pdf.image(outDir+'/'+'hwe_'+str(key)+'.png', x=10, y=50, w=190, h=150)
			
	
	return pdf

def pruning(dictLD):
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "LD Pruning", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', '', 12)
	pdf.multi_cell(0, 5, 'Pruning SNPs based upon linkage equilibrium is senstive to race/ethnicity.  Therefore, LD-pruning is performed \
		independently on each ethnic group in the data set.', 0, 1, 'J')
	pdf.multi_cell(0, 5, '\n', 0, 1, 'J')
	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)
	
	# iterate through all ethnic groups for LD pruning stats
	for key, value in dictLD.iteritems():
		pdf.multi_cell(0, 8, str(key), 1, 'L', True)
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs analyzed:  ' +  str(value[0]), 1, 1, 'L')
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs Passing:  ' +  str(value[1]) + ' (' + str("%.2f" % round((float(value[1])/float(value[0]))*100, 2)) + '%)', 1, 1, 'L')
		pdf.multi_cell(0, 8, '\n\n', 0, 1, 'J')


	return pdf

def minor_allele_freq(dictMAF, thresh):
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Minor Allele Frequency", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)

	# iterate through all ethnic groups for LD pruning stats
	for key, value in dictMAF.iteritems():
		pdf.multi_cell(0, 8, str(key), 1, 'L', True)
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs analyzed:  ' +  str(value[0]), 1, 1, 'L')
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs MAF >= ' + str(thresh) + ': ' +  str(value[1]) + ' (' + str("%.2f" % round((float(value[1])/float(value[0]))*100, 2)) + '%)', 1, 1, 'L')
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs MAF <= ' + str(thresh) + ': ' +  str(value[2]) + ' (' + str("%.2f" % round((float(value[2])/float(value[0]))*100, 2)) + '%)', 1, 1, 'L')
		pdf.multi_cell(0, 8, '\n\n', 0, 1, 'J')


	return pdf

def heterozygosity(het_method, std, het_dataframe, thresh, minThresh, population, outDir):
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Heterozygosity", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)

	sample_fails = open(outDir + '/' + population + '/samples_failing_heterozygosity.txt', 'w')
	#het_method of setting min and max thresholds based on F score inbreeding coefficient
	if het_method == 'minMax':
		fail_het = het_dataframe.loc[((het_dataframe['F'] > thresh) | (het_dataframe['F'] < minThresh))]
		fail_het[['FID', 'IID']].to_csv(sample_fails.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>


		# sorts inbreeing coeffiecient score to rank for ploting
		het_dataframe.sort_values(by='F', ascending=True, axis=0, inplace=True)
		het_dataframe = het_dataframe.reset_index()
		het_dataframe['rank'] = het_dataframe.index
		
		pdf.set_font('Arial', '', 12)
		pdf.multi_cell(0, 5, 'Heterozygosity is calculated on the merged sample space.  Samples with extreme levels of heterozygosity as measured by \
			the F inbreeding coefficient, may be indicative of a sample exhibiting poor quality.  The current inbreeding coefficient threshold is set \
			to a maximum of ' + str(thresh) + ' and a minimum of ' + str(minThresh) +'.  Any F value above the maximum threshold or below the minimum threshold is filtered out.', 0, 1, 'J')
		pdf.multi_cell(0, 5, '\n', 0, 1, 'J')
		pdf.set_font('Arial', 'B', 16)
		pdf.set_fill_color(200)
		pdf.multi_cell(0, 8, 'Total Number Samples Analyed: '+ str(len(het_dataframe.index)), 1, 'L', True)
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number Samples Failing: '+ str(len(fail_het.index)), 1, 1, 'L')

		het_dataframe.plot(kind='scatter', x='rank', y='F', title='Ranked Inbreeding Coefficient scores', s=7)
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+ str(population) + '_heterozygosity_plot.png')
		plt.close()
		pdf.image(outDir+'/'+ str(population) + '_heterozygosity_plot.png', x=10, y=130, w=190, h=150)

		sample_fails.flush()
		sample_fails.close()
	
	elif het_method == 'meanStd':
		het_mean = het_dataframe['het_score'].mean() # calculates the mean of the calculated het_score across all samples
		het_std = het_dataframe['het_score'].std() # calculated the standard deviation of the calculated het_score across all samples
		max_value = het_mean + (int(std)*het_std) # gets the max acceptable het_score give the number of standard deviations allowed from the mean set by the user
		min_value = het_mean - (int(std)*het_std) # gets the min acceptable het_score give the number of standard deviations allowed from the mean set by the user
		het_dataframe['Group'] = 'Before'
		fail_het = het_dataframe.loc[((het_dataframe['het_score'] > max_value) | (het_dataframe['het_score'] < min_value))]
		fail_het[['FID', 'IID']].to_csv(sample_fails.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>
		pass_het = het_dataframe.loc[((het_dataframe['het_score'] <= max_value) | (het_dataframe['het_score'] >= min_value))]
		pass_het['Group'] = 'After'
		
		combine = pandas.concat(['het_dataframe', 'pass_het']) # concantenates dataframes before and after filterings to make boxplot comparisons by Group
		het_dataframe.plot(kind='box', x='Group', y='het_score', title='Heterozygosity Scores')
		plt.tight_layout(pad=2, w_pad=2, h_pad=2)
		plt.savefig(outDir+'/'+ str(population) +'_heterozygosity_plot.png')
		plt.close()

		# generate PDF of heterozygosity results
		pdf.set_font('Arial', '', 12)
		pdf.multi_cell(0, 5, 'Heterozygosity is calculated on the merged sample space.  Samples with extreme levels of heterozygosity as measured by \
			the heterozygosity score calculated using the following formula: het_score = 1 - [O(HOM)/N(NM)], may be indicative of a sample exhibiting poor quality.  The current het_score thresholds \
			are set to filter out any samples exhibiting excess heterozygosity of +/- ' + str(std) + ' from the mean het_score.' , 0, 1, 'J')
		pdf.multi_cell(0, 5, '\n', 0, 1, 'J')
		pdf.set_font('Arial', 'B', 16)
		pdf.set_fill_color(200)
		pdf.multi_cell(0, 8, 'Total Number Samples Analyed: '+ str(len(het_dataframe.index)), 1, 'L', True)
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number Samples Failing: '+ str(len(fail_het.index)), 1, 1, 'L')
		pdf.image(outDir+'/'+ str(population) +'_heterozygosity_plot.png', x=10, y=130, w=190, h=150)
	
		sample_fails.flush()
		sample_fails.close()

	return sample_fails.name, pdf


def relatedness(ibd_dataframe, outDir):
	
	dups_text = open(outDir + '/' + 'duplicate_pairs.txt', 'w') # outputs pairs with Z0, Z1, Z2	score
	remove_dups = open(outDir + '/remove_all_duplicate_pairs.txt', 'w') # outputs duplicate samples for PLINK format removal
	
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Relatedness", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)


	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)
	pdf.multi_cell(0, 10, 'Total Number of Sample Pairs Analyzed:  ' +  str(len(ibd_dataframe.index)), 1, 'L', True)

	
	duplicates = ibd_dataframe.loc[ibd_dataframe['Z2'] > 0.97]
	parent_offspring = ibd_dataframe.loc[(ibd_dataframe['Z1'] > 0.97) & (ibd_dataframe['Z0'] < 0.05) & (ibd_dataframe['Z2'] < 0.05)]
	full_sibs = ibd_dataframe.loc[(ibd_dataframe['Z0'] < 0.40) & (ibd_dataframe['Z2'] > 0.16)]
	second_degree = ibd_dataframe.loc[(ibd_dataframe['Z0'] < 0.60) & (ibd_dataframe['Z1'] < 0.58) & (ibd_dataframe['Z2'] < 0.05)]
	unrelated = ibd_dataframe.loc[ibd_dataframe['Z0'] > 0.78]
	
	# format data so it can be put into usable format by PLINK --remove
	first_in_pair = duplicates[['FID1', 'IID1']]
	second_in_pair = duplicates[['FID2', 'IID2']]
	first_in_pair.columns = ['FID', 'IID']
	second_in_pair.columns = ['FID', 'IID']
	merged_removals = first_in_pair.merge(second_in_pair, how='outer', on=['FID', 'IID'])
	merged_removals[['FID', 'IID']].to_csv(remove_dups.name, sep='\t', index=False, header=False) # output file created to PLINK --remove


	pdf.set_font('Arial', '', 16)
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of duplicate pairs:  '+str(len(duplicates.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of parent-offspring pairs:  '+str(len(parent_offspring.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of full siblings pairs:  '+str(len(full_sibs.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of 2nd degree pairs:  '+str(len(second_degree.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of unrelated pairs:  '+str(len(unrelated.index)), 1, 1, 'L')


	dups_text.flush()
	dups_text.close()

	remove_dups.flush()
	remove_dups.close()

	return pdf, remove_dups.name

