import sys
import os
import math
# from bx.bbi.bigwig_file import BigWigFile

def bisulfite_conversion_rate(input):
	"""This function is used to calculate the bisulfite conversion rate based on 
	lambda DNA."""
	logfile = open("bilsulfite_conversion_rate.txt",'a')
	for file in input:
		total_c = 0
        	methyl_c = 0
		for line in open(file,'r').xreadlines():
			line = line.strip().split('\t')
			if line[0] == 'chr':
				continue
			else:
				total_c += float(line[5])
				methyl_c += float(line[6])
	
		print >>logfile, '\t'.join([file,str(1-methyl_c/total_c)])
	logfile.close()

def bisulfite_conversion_rate_mcall(input):
	"""This function is used to calculate the bisulfite conversion rate based on 
	lambda DNA."""
	logfile = open("bilsulfite_conversion_rate_mcall.txt",'a')
	for file in input:
		total_c = 0
        	methyl_c = 0
		for line in open(file,'r').xreadlines():
			line = line.strip().split('\t')
			if line[0] == '#chrom':
				continue
			else:
				total_c += float(line[4])
				methyl_c += float(line[5])
	
		print >>logfile, '\t'.join([file,str(1-methyl_c/total_c)])
	logfile.close()

def average_mratio_coverage(input):
	"""Calculate the average CpG coverage of mratio output.
	"""
	trans = {'-':'R','+':'F'}
	
	logfile = open("average_mratio_coverage.txt",'w')
	for file in input:
		count_1, sum_1 = 0,0
		count_5, sum_5 = 0,0
		count_10, sum_10 = 0,0
		outf = open(file.split('.')[0]+'.bsmap.CpG.txt', 'w')
		print >>outf, 'chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT'
		for line in open(file,'r').xreadlines():
			line = line.strip().split('\t')
			if line[3] == 'CG':
				if float(line[5]) >= 1:
					count_1 += 1
					sum_1 += float(line[5])
				if float(line[5]) >= 5:
					count_5 += 1
					sum_5 += float(line[5])
				if float(line[5]) >= 10:
					count_10 += 1
					sum_10 += float(line[5])
				print >>outf, '\t'.join([line[0]+'.'+line[1], line[0], line[1], trans[line[2]], line[7], str(100*float(line[4])), str(100-100*float(line[4]))])
		outf.close()
		print >>logfile, '\t'.join([str(t) for t in [file,count_1,sum_1/count_1,count_5,sum_5/count_5,count_10,sum_10/count_10]])
	logfile.close()

def average_mcall_coverage(input):
	"""Calculate the average CpG coverage of mcall output.
	"""
	trans = {'-':'R','+':'F'}
	number = {'-':1,'+':1,'B':2}
	
	logfile = open("average_mcall_coverage.txt",'a')
	for file in input:
		count_1, sum_1 = 0,0
		count_5, sum_5 = 0,0
		count_10, sum_10 = 0,0
		outf = open(file[:-10]+'.bsmap.CpG.txt', 'w')
		print >>outf, 'chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT'
		for line in open(file,'r').xreadlines():
			line = line.strip().split('\t')
			if line[0][0] == '#':
				continue
			if float(line[4]) >= 1:
				count_1 += number[line[6]]
				sum_1 += float(line[4])
			if float(line[4]) >= 5:
				count_5 += number[line[6]]
				sum_5 += float(line[4])
			if float(line[4]) >= 10:
				count_10 += number[line[6]]
				sum_10 += float(line[4])
			print >>outf, '\t'.join([line[0]+'.'+str(int(line[1])+1), line[0], str(int(line[1])+1), line[6], line[4], line[3], str(1-float(line[3]))])
		outf.close()
		print >>logfile, '\t'.join([str(t) for t in [file,count_1,sum_1/count_1,count_5,sum_5/count_5,count_10,sum_10/count_10]])
	logfile.close()
	
def filter_unite_reads(file_list, loc = 3, hip = 0.999):
	"""Filter low coverage (< loc) and high coverage (> hip%) reads.
	Unite all samples based on the effective loci (> eff_loci).
	"""
	all_methyl = {}
	title = []
	for i in range(len(file_list)):
		methyl = {}
		title.append(file_list[i].split('/')[-1].split('.')[0]+'_coverage')
		title.append(file_list[i].split('/')[-1].split('.')[0]+'_ratio')	
		for line in open(file_list[i],'r').xreadlines():
			line = line.strip().split('\t')
			if line[0] == 'chrBase':
				continue
			else:
				# methyl[line[1]+'_'+line[2]] = [int(line[3]), float(line[4]), float(line[5])]
				methyl[line[1]+'_'+line[2]+'_'+line[3]] = [int(line[4]), float(line[5]), float(line[6])]
	
		coverage = zip(*methyl.values())[0]
		if int(len(coverage)*hip)+1 >= int(len(coverage)):
			hip_count = sorted(coverage)[-1]
		else:
			hip_count = sorted(coverage)[int(len(coverage)*hip)+1]
		for k,v in methyl.iteritems():
			if v[0] <= loc:
				continue
			if v[0] >= hip_count:
				continue
			if all_methyl.has_key(k):
				all_methyl[k][2*i+1] = v[1]
				all_methyl[k][2*i] = v[0]
			else:
				all_methyl[k] = ['NA']*len(file_list)*2
				all_methyl[k][2*i+1] = v[1]
				all_methyl[k][2*i] = v[0]
	
	outf = open("all.merged.values.txt",'w')
	# print >>outf, "chr\tpos\t"+'\t'.join(title)
	print >>outf, "chr\tpos\tstrand\t"+'\t'.join(title)
	for k,v in all_methyl.iteritems():
		print >>outf, k.replace('_','\t')+'\t'+'\t'.join([str(t) for t in v])
	outf.close()

def filter_unite_reads_samGbed(output, file_list, loc = 3, hip = 0.999):
	"""Filter low coverage (< loc) and high coverage (> hip%) reads.
	Unite all samples based on the effective loci (> eff_loci).
	"""
	all_methyl = {}
	title = []
	for i in range(len(file_list)):
		methyl = {}
		title.append(file_list[i].split('/')[-1].split('.')[0]+'_coverage')
		title.append(file_list[i].split('/')[-1].split('.')[0]+'_ratio')	
		for line in open(file_list[i],'r').xreadlines():
			line = line.strip().split('\t')
			if line[0] == '#chrom':
				continue
			else:
				methyl[line[0]+'_'+line[1]+'_'+line[6]] = [int(line[4]), float(line[3])]#methyl[#Chrom_Start_Strand] = [TotalC(Coverage), freqC(MethylationLevel)]
	
		coverage = zip(*methyl.values())[0]
		if int(len(coverage)*hip)+1 >= int(len(coverage)):
			hip_count = sorted(coverage)[-1]
		else:
			hip_count = sorted(coverage)[int(len(coverage)*hip)+1]
		for k,v in methyl.iteritems():
			if v[0] <= loc:
				continue
			if v[0] >= hip_count:
				continue
			if all_methyl.has_key(k):
				all_methyl[k][2*i+1] = v[1]
				all_methyl[k][2*i] = v[0]
			else:
				all_methyl[k] = ['NA']*len(file_list)*2
				all_methyl[k][2*i+1] = v[1]
				all_methyl[k][2*i] = v[0]
	
	outf = open(output,'w')
	# print >>outf, "chr\tpos\t"+'\t'.join(title)
	print >>outf, "#Chrom\tStart\tStrand\t"+'\t'.join(title)
	for k,v in all_methyl.iteritems():
		print >>outf, k.replace('_','\t')+'\t'+'\t'.join([str(t) for t in v])
	outf.close()

def add_CpG_density():
	"""Add CpG density info to all the detected CpG sites."""
	CpG = {}
	for line in open('/mnt/Storage/home/wangcf/annotations/mm9.CpG.density.txt','r').xreadlines():
		line = line.strip().split('\t')
		CpG[line[0]+'_'+line[1]] = line[2]

	outf = open("all.merged.values.density.txt",'w')
	for line in open("all.merged.values.txt",'r').xreadlines():
		line = line.strip().split('\t')
		if line[0] == 'chr':
			print >>outf, '\t'.join(line[:3]+['density']+line[3:])
		else:
			if CpG.has_key(line[0]+'_'+line[1]):
				print >>outf, '\t'.join(line[:3]+[CpG[line[0]+'_'+line[1]]]+line[3:])
			else:
				print "CpG not found:",line[:3]
	outf.close()

def remove_strand_info():
	"""Merge +/-/B position into single position for downstream analysis."""
	CpG = {}
	outf = open("all.merged.values.strandness.txt",'w')
	for line in open("all.merged.values.density.txt",'r').xreadlines():
		line = line.strip().split('\t')
		if line[0] == 'chr':
			print >>outf, '\t'.join(line[:2]+line[3:])
		else:
			if CpG.has_key(line[0]+'_'+line[1]):
				CpG[line[0]+'_'+line[1]].append(line[3:])
			else:
				CpG[line[0]+'_'+line[1]] = [line[3:]]

	for k,v in CpG.iteritems():
		if len(v) == 1:
			print >>outf, k.replace('_','\t')+'\t'+'\t'.join(v[0])
		else:
			v = zip(*v)
			output = [v[0][0]]
			for i in range(1,len(v),2):
				output.extend(average_single_site(v[i+1],v[i]))
			print >>outf, k.replace('_','\t')+'\t'+'\t'.join(output)
	outf.close()

def average_single_site(rl, cl):
	"""Average the mratio of the input list, import by average_all_site."""

	tc, mc = 0,0
	if isinstance(rl,str):
		return [cl, rl]
	for i in range(len(rl)):
		if not rl[i] == 'NA':
			tc += float(cl[i])
			mc += float(cl[i])*float(rl[i])
	if not tc == 0:
		return [str(tc),str(mc/tc)]
	else:
		return ["NA","NA"]

def main():
	parameter = sys.argv[1]
	if parameter == 'ConversionRate':
		bisulfite_conversion_rate(sys.argv[2:])
	if parameter == 'ConversionRate_mcall':
		bisulfite_conversion_rate_mcall(sys.argv[2:])
	elif parameter == 'MratioCoverage':
		average_mratio_coverage(sys.argv[2:])
	elif parameter == 'McallCoverage':
		average_mcall_coverage(sys.argv[2:])
	elif parameter == 'FilterUniteReads':
		filter_unite_reads(sys.argv[2:])
	elif parameter == 'UniteReadsSamGbed':
		filter_unite_reads_samGbed(sys.argv[2], sys.argv[3:])
	elif parameter == 'AddCpGDensity':
		add_CpG_density()
	elif parameter == 'RemoveStrandInfo':
		remove_strand_info()
	
if __name__ == "__main__":
	main()