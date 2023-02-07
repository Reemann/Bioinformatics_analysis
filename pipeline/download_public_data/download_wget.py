import os
import urllib.request, urllib.error, urllib.parse
from optparse import OptionParser
import re
import linecache

'''
	download sra from SRA/EBI/DDBJ
	python download.py -d dataSource -i dataId -t dataType -o outputInfo
	python download.py -d GEO -i GSMID -t PAIRED -o test
'''

def download(dataSource, dataId, outputInfo, dataType):
	if dataSource == 'GEO':
		# SRR_pre, SRR_mid, SRR = dataId.split(';')
		# dataId: GSM ID
		GSM_link = 'https://www.ncbi.nlm.nih.gov/sra/?term=%s'%dataId
		response = urllib.request.urlopen(GSM_link)
		html = response.read().decode('utf-8')
		node_all = [m.start() for m in re.finditer('run=SRR', html)]
		if len(node_all) == 0:
			print("NoData")
		elif len(node_all) == 1: # only one SRR file
			SRR = html[node_all[0]:node_all[0]+60].split('">')[0].split('run=')[1] #SRRnum
			SRRsix = str(SRR)[0:6]
			SRRthree = str(SRR)[9:].zfill(3)
			if dataType == 'PAIRED':
				if len(SRR) == 9:
					cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_1.fastq.gz -O %s_R1.fq.gz'%(SRRsix, SRR, SRR, outputInfo)
					cmd_R2 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_2.fastq.gz -O %s_R2.fq.gz'%(SRRsix, SRR, SRR, outputInfo)
				else:
					cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s_1.fastq.gz -O %s_R1.fq.gz'%(SRRsix, SRRthree, SRR, SRR, outputInfo)
					cmd_R2 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s_2.fastq.gz -O %s_R2.fq.gz'%(SRRsix, SRRthree, SRR, SRR, outputInfo)
				os.system(cmd_R1)
				os.system(cmd_R2)
			elif dataType == 'SINGLE':
				if len(SRR) == 9:
					cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s.fastq.gz -O %s_R1.fq.gz'%(SRRsix, SRR, SRR, outputInfo)
				else:
					cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s.fastq.gz -O %s_R1.fq.gz'%(SRRsix, SRRthree, SRR, SRR, outputInfo)
				os.system(cmd_R1)
			else:
				print('Please choose the type of data: PAIRED or SINGLE!!!')
			# print(SRR)
		else:
			for i in range(len(node_all)):
				SRR = html[node_all[i]:node_all[i]+60].split('">')[0].split('run=')[1] #SRRnum
				SRRsix = str(SRR)[0:6]
				SRRthree = str(SRR)[9:].zfill(3)
				if dataType == 'PAIRED':
					if len(SRR) == 9:
						cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_1.fastq.gz -O %s_part%s_R1.fq.gz'%(SRRsix, SRR, SRR, outputInfo, i+1)
						cmd_R2 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_2.fastq.gz -O %s_part%s_R2.fq.gz'%(SRRsix, SRR, SRR, outputInfo, i+1)
					else:
						cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s_1.fastq.gz -O %s_part%s_R1.fq.gz'%(SRRsix, SRRthree, SRR, SRR, outputInfo, i+1)
						cmd_R2 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s_2.fastq.gz -O %s_part%s_R2.fq.gz'%(SRRsix, SRRthree, SRR, SRR, outputInfo, i+1)
					os.system(cmd_R1)
					os.system(cmd_R2)
				elif dataType == 'SINGLE':
					if len(SRR) == 9:
						cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s.fastq.gz -O %s_part%s_R1.fq.gz'%(SRRsix, SRR, SRR, outputInfo, i+1)
					else:
						cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s.fastq.gz -O %s_part%s_R1.fq.gz'%(SRRsix, SRRthree, SRR, SRR, outputInfo, i+1)
					os.system(cmd_R1)
				else:
					print('Please choose the type of data: PAIRED or SINGLE!!!')
				# print(SRR)
			# cmd_cat_R1='cat %s_part*_R1.fq.gz > %s_R1.fq.gz && rm %s_part*_R1.fq.gz'%(outputInfo, outputInfo, outputInfo)
			cmd_cat_R1='cat %s_part*_R1.fq.gz > %s_R1.fq.gz'%(outputInfo, outputInfo)
			os.system(cmd_cat_R1)
			if dataType == 'PAIRED':
				# cmd_cat_R2='cat %s_part*_R2.fq.gz > %s_R2.fq.gz && rm %s_part*_R2.fq.gz'%(outputInfo, outputInfo, outputInfo)
				cmd_cat_R2='cat %s_part*_R2.fq.gz > %s_R2.fq.gz'%(outputInfo, outputInfo)
				os.system(cmd_cat_R2)
	elif dataSource == 'EBI': # dataId: ERR522956
		ERR_pre = dataId[0:6]
		ERR = dataId
		if dataType == 'PAIRED':
			cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_1.fastq.gz -O %s_R1.fq.gz'%(ERR_pre, ERR, ERR, outputInfo)
			cmd_R2 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_2.fastq.gz -O %s_R2.fq.gz'%(ERR_pre, ERR, ERR, outputInfo)
			os.system(cmd_R1)
		elif dataType == 'SINGLE':
			cmd_R1 = 'wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s.fastq.gz -O %s_R1.fq.gz'%(ERR_pre, ERR, ERR, outputInfo)
			os.system(cmd_R1)
		else:
			print('Please choose the type of data: PAIRED or SINGLE!!!')
	elif dataSource == 'DDBJ':
		ERA_pre, ERA, ERX, ERR = dataId.split(';')
		if dataType == 'PAIRED':
			# https://ddbj.nig.ac.jp/public/ddbj_database/biosample/
			cmd_R1 = 'wget -c ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/%s/%s/%s/%s_1.fastq.bz2 %s_R1.fq.bz2'%(ERA_pre, ERA, ERX, ERR, outputInfo)
			cmd_R2 = 'wget -c ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/%s/%s/%s/%s_2.fastq.bz2 %s_R2.fq.bz2'%(ERA_pre, ERA, ERX, ERR, outputInfo)
		elif dataType == 'SINGLE':
			cmd_R1 = 'wget -c ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/%s/%s/%s/%s.fastq.bz2 %s_R1.fq.bz2'%(ERA_pre, ERA, ERX, ERR, outputInfo)
		else:
			print('Please choose the type of data: PAIRED or SINGLE!!!')
	else:
		print('Please choose the source of data: GEO, EBI or DDBJ!!!')

def main():
	parser = OptionParser()
	parser.add_option('-d','--dataSource',dest='dataSource',help='Source of data: GEO, EBI or DDBJ')
	parser.add_option('-i','--dataId',dest='dataId',help='Download information')
	parser.add_option('-t', '--dataType', dest = 'dataType', help = 'Data type: PAIRED or SINGLE')
	parser.add_option('-o','--outputInfo',dest='outputInfo',help='Output file name')
	(options,args) = parser.parse_args()
	download(options.dataSource, options.dataId, options.outputInfo, options.dataType)

main()