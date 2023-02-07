#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

def region_methylation_covOfEachCpG(region_G_bed,outputF):
	'''
		The methylation level of each annotated genomic region in each samples was measured as
		the sum of the methylation level of every CpG site divided by the total number of the CpG sites.
		Genome Research 2013 Guo.

		Note: Require each CpG with coverage at least 3 !!!
	'''
	region_dict = {}
	inf = open(region_G_bed)
	for line in inf:
		ele = line.strip().split('\t')
		index = '@'.join(ele[0:3])
		if index not in region_dict:
			region_dict[index] = ['NA','NA']
			if ele[7]!='.' and ele[8]!='.' and int(ele[7])>=3: #minimum coverage is 3
				region_dict[index][0] = int(ele[7]) #total C
				region_dict[index][1] = int(ele[8]) #methyl C
		else:
			if ele[7]!='.' and ele[8]!='.' and int(ele[7])>=3:
				if region_dict[index][0] == "NA":
					region_dict[index][0] = int(ele[7]) #total C
					region_dict[index][1] = int(ele[8]) #methyl C
				else:
					region_dict[index][0] += int(ele[7]) #total C
					region_dict[index][1] += int(ele[8]) #methyl C
		# print('%s\t%s\t%s\n'%(index,region_dict[index][0],region_dict[index][1]))
	inf.close()

	outputF = open(outputF,'w')
	for key in region_dict:
		location = key.split('@')
		if region_dict[key][1]!='NA' and region_dict[key][0]!='NA':
			outputF.write('%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],float(region_dict[key][1])/float(region_dict[key][0])))
		else:
			outputF.write('%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],'NA'))

	outputF.close()

def region_methylation(region_G_bed, outputF):
	'''
		The methylation level of each annotated genomic region in each sample was measured as
		the number of methylated CpG by the total number of CpG sites.
		Genome Research 2013 Guo.

		Note: Coverage requirment is restricted to each annotated genomic region!!!
	'''
	region_dict = dict()
	inFH = open(region_G_bed, 'r')
	for line in inFH:
		ele = line.strip().split('\t')
		index = '@'.join(ele[0:3])
		if index not in region_dict:
			region_dict[index] = ['NA', 'NA'] # Initiation
			if ele[7]!='.' and ele[8]!='.':
				region_dict[index][0] = int(ele[7]) #total C
				region_dict[index][1] = int(ele[8]) #methyl C
		else:
			region_dict[index][0] += int(ele[7]) #total C
			region_dict[index][1] += int(ele[8]) #methyl C
	inFH.close()

	outFH = open(outputF, 'w')
	for key in region_dict:
		location = key.split('@')
		if region_dict[key][0]!='NA' and region_dict[key][1]!='NA':
			if region_dict[key][0]>=3:
				outFH.write('%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],float(region_dict[key][1])/float(region_dict[key][0])))
			else:
				outFH.write('%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],'NA'))
		else:
			outFH.write('%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],'NA'))
	outFH.close()


def region_methylation_coverage(region_G_bed, outputF):
	'''
		D: Aug-02-2021 16:12 Mon

		The methylation level of each annotated genomic region in each sample was measured as
		the number of methylated CpG by the total number of CpG sites.
		Genome Research 2013 Guo.

		Note: Coverage requirment is restricted to each annotated genomic region!!!
		Provide the total coverage information for each annotated genomic region!!!
	'''
	region_dict = dict()
	inFH = open(region_G_bed, 'r')
	for line in inFH:
		ele = line.strip().split('\t')
		index = '@'.join(ele[0:3])
		if index not in region_dict:
			region_dict[index] = ['NA', 'NA'] # Initiation
			if ele[7]!='.' and ele[8]!='.':
				region_dict[index][0] = int(ele[7]) #total C
				region_dict[index][1] = int(ele[8]) #methyl C
		else:
			region_dict[index][0] += int(ele[7]) #total C
			region_dict[index][1] += int(ele[8]) #methyl C
	inFH.close()

	outFH = open(outputF, 'w')
	for key in region_dict:
		location = key.split('@')
		if region_dict[key][0]!='NA' and region_dict[key][1]!='NA':
			if region_dict[key][0]==0:
				outFH.write('%s\t%s\t%s\t%s\t%d\n'%(location[0],location[1],location[2],'NA', 0))
			else:
				outFH.write('%s\t%s\t%s\t%s\t%d\n'%(location[0],location[1],location[2],float(region_dict[key][1])/float(region_dict[key][0]),region_dict[key][0]))
		else:
			outFH.write('%s\t%s\t%s\t%s\t%d\n'%(location[0],location[1],location[2],'NA', 0))
	outFH.close()


def region_methylationCR(region_G_bed,outputF,index_coverage,index_ratio):
	region_dict = {}
	inf = open(region_G_bed)
	for line in inf:
		ele = line.strip().split('\t')
		index = '@'.join(ele[0:3])
		if index not in region_dict:
			region_dict[index] = ['NA',0.0]
			if ele[3]!='.' and ele[index_coverage] !="NA" and ele[index_ratio]!="NA":
				if float(ele[index_coverage])>=3.0: #minimum coverage is 3
					region_dict[index][0] = float(ele[index_ratio])
		else:
			if ele[3]!='.' and ele[index_coverage] !="NA" and ele[index_ratio]!="NA":
				if float(ele[index_coverage])>=3.0:
					if region_dict[index][0] == "NA":
						region_dict[index][0] = float(ele[index_ratio])
					else:
						region_dict[index][0] += float(ele[index_ratio])
						region_dict[index][1] += 1.0 #methyl C
	inf.close()

	outputF = open(outputF,'w')
	for key in region_dict:
		location = key.split('@')
		if region_dict[key][0]!='NA' and region_dict[key][1] != 0.0:
			outputF.write('%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],float(region_dict[key][0])/float(region_dict[key][1])))
		else:
			outputF.write('%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],'NA'))

	outputF.close()

def region_methylation_pandas(region_G_bed,outputF):
	outFH = open(outputF,'a')
	methylF = pd.read_table(region_G_bed,sep='\t',header=None,usecols=list(range(0,3))+[6])
	methylF.index(methylF[[0,1,2]].apply(lambda x: '_'.join(str(value) for value in x),axis=1))
	regions = methylF.index.unique()
	for ele in range(len(regions)):
		each = methylF.loc[regions[ele]]
		if(each[6]=='.'):
			np.savetxt(outFH,np.array([ele,'NA']))
		else:
			np.savetxt(outFH,np.array([ele,each[6].mean(axis=0,skipna=True)]))
	outFH.close()


def region_methylation_label(region_G_bed,outputF):
	'''
		The methylation level of each annotated genomic region in each samples was measured as
		the sum of the methylation level of every CpG site divided by the total number of the CpG sites.
		Genome Research 2013 Guo.
	'''
	region_dict = {}
	inf = open(region_G_bed)
	for line in inf:
		ele = line.strip().split('\t')
		index = '@'.join(ele[0:4])
		if index not in region_dict:
			region_dict[index] = ['NA','NA']
			if ele[8]!='.' and ele[9]!='.' and int(ele[8])>=3: #minimum coverage is 3
				region_dict[index][0] = int(ele[8]) #total C
				region_dict[index][1] = int(ele[9]) #methyl C
		else:
			if ele[8]!='.' and ele[9]!='.' and int(ele[8])>=3:
				if region_dict[index][0] == "NA":
					region_dict[index][0] = int(ele[8]) #total C
					region_dict[index][1] = int(ele[9]) #methyl C
				else:
					region_dict[index][0] += int(ele[8]) #total C
					region_dict[index][1] += int(ele[9]) #methyl C
		# print('%s\t%s\t%s\n'%(index,region_dict[index][0],region_dict[index][1]))
	inf.close()

	
	outputF = open(outputF,'w')
	for key in region_dict:
		location = key.split('@')
		if region_dict[key][1]!='NA' and region_dict[key][0]!='NA':
			outputF.write('%s\t%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],location[3],float(region_dict[key][1])/float(region_dict[key][0])))
		else:
			outputF.write('%s\t%s\t%s\t%s\t%s\n'%(location[0],location[1],location[2],location[3],'NA'))

	outputF.close()

def unit(file_list):
	all_methyl = {}
	title = []
	for i in range(len(file_list)):
		methyl = {}
		title.append(file_list[i].split('/')[-1].split('.bed')[0])
		for line in open(file_list[i],'r'):
			line = line.strip().split('\t')
			if line[0] == '#chrom':
				continue
			else:
				if int(line[4])>=3:
					methyl[line[0]+'_'+line[1]+'_'+line[2]] = line[3]
				else:
					methyl[line[0]+'_'+line[1]+'_'+line[2]] = "NA"

		for k,v in methyl.items():
			if k in all_methyl:
				all_methyl[k][i] = v
			else:
				all_methyl[k] = ['NA']*len(file_list)
				all_methyl[k][i] = v

	outf = open("all.merged.values.txt",'w')
	print("#Chrom\tPos\tStrand\t"+'\t'.join(title), file=outf)
	for k,v in all_methyl.items():
		print(k.replace('_','\t')+'\t'+'\t'.join([str(t) for t in v]), file=outf)
	outf.close()

def main():
	functName = sys.argv[1]
	if functName == 'region_methylation_covOfEachCpG':
		region_methylation_covOfEachCpG(sys.argv[2], sys.argv[3])
	if functName == 'region_methylation':
		region_methylation(sys.argv[2], sys.argv[3])
	if functName == "region_methylation_coverage":
		region_methylation_coverage(sys.argv[2], sys.argv[3])
	if functName == 'region_methylationCR':
		region_methylationCR(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
	if functName == 'region_methylation_pandas':
		region_methylation_pandas(sys.argv[2], sys.argv[3])
	if functName == 'region_methylation_label':
		region_methylation_label(sys.argv[2], sys.argv[4])
	if functName == 'unit':
		unit(sys.argv[2:])

main()
