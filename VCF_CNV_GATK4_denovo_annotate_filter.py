#!/bin/python

import sys
import os
import gzip 


## Get arguments from snakemake

#PED_file=snakemake.input[0]
database=snakemake.input[0]
vcf=snakemake.input[1]

### Create dictionary from the pedigree file provided

dic_ped={}
#IN=open(PED_file,"r")
#a=0
#for line in IN:
#    a+=1
#    if a>1:
#        field=line.split('\t')
#        if field[5].rstrip() == "2":
#            dic_ped[field[1]]=[field[2],field[3]]
#IN.close()

### Iterate over vcf files
 

if vcf.endswith('vcf.gz'):
	l=gzip.open(vcf,'rt')
elif vcf.endswith('vcf'):
	l=open(vcf,'r')
else:
	print(vcf,"hasn't a compatible extension: vcf or vcf.gz")
	exit()

OUT=gzip.open(vcf.replace(".vcf.gz","_filtered.vcf.gz"),"wt")

for line in l:
	if line[0] == "#":
		OUT.write(line)
		if line[0:6] == "#CHROM":
			name=line.split('\t')[9].rstrip()

	if line[0] != "#":
		field=line.split('\t')
		chr=field[0]
		cnvStart=int(field[1])
		cnvStop=int(field[7].split('=')[1])
		size_cnv=int(cnvStop)-int(cnvStart)
		info=field[8].split(':')
		index_GT=info.index("GT")
		data=field[9]
		data_GT=data[index_GT]
		
		if data_GT == "0":
			continue
		elif data_GT == "1":
			ALT="DEL"
		elif data_GT == "2":
			ALT="DUP"
		
		Numb_samples_match=0
		
		if name in dic_ped.keys():
			INDEX_PED=True
			Father_name=dic_ped[name][0]
			Mother_name=dic_ped[name][1]
			Father="FALSE"
			Father_key="NA"
			Mother="FALSE"
			Mother_key="NA"
		else:
			INDEX_PED=False
			Father="NA"
			Father_key="NA"
			Mother="NA"
			Mother_key="NA"
			
		IN=open(database,'r')
		chr_match=False
		
		for line_db in IN:
			if line_db[0] != "#":
				match=False
				
				field_db=line_db.split('\t')
				targetChr=field_db[0]
				targetStart=int(field_db[1])
				targetStop=int(field_db[2])
				cnv_type=field_db[3]
				number_samples=int(field_db[4])
				name_samples=field_db[5].rstrip().split(',')
				
				if targetChr == chr:
					chr_match=True
					
					if cnv_type == ALT and ((cnvStart <= targetStart and targetStop <= cnvStop) or (targetStart <= cnvStart and cnvStart <= targetStop and targetStop <= cnvStop) or (targetStart >= cnvStart and targetStop >= cnvStop and targetStart <= cnvStop) or (cnvStart >= targetStart and cnvStop <= targetStop)):
					
						if (targetStart <= cnvStart and targetStop >= cnvStop):
							match=True

						elif (targetStart <= cnvStart and targetStop > cnvStart and targetStop < cnvStop):
							if ((targetStop - cnvStart) > (0.7 * size_cnv)):
								match=True

						elif (targetStart >= cnvStart and targetStart < cnvStop and targetStop > cnvStart and targetStop <= cnvStop):
							if (size_cnv - ((targetStart - cnvStart) + (cnvStop - targetStop))) > (0.7 * size_cnv):
								match=True
						
						elif (targetStart > cnvStart and targetStart < cnvStop and targetStop > cnvStop):
							if (cnvStop - targetStart) > (0.7 * size_cnv):
								match=True

					if match:
						if name in name_samples:
							number_samples=number_samples-1
						Numb_samples_match+=number_samples
						if INDEX_PED:
							if Father_name in name_samples: 
								Father="TRUE"
								Father_key=str(targetChr)+":"+str(targetStart)+"-"+str(targetStop)
							if Mother_name in name_samples:
								Mother="TRUE"
								Mother_key=str(targetChr)+":"+str(targetStart)+"-"+str(targetStop)
							
				elif targetChr != chr and chr_match:
					break
					
		IN.close()
		
		## If no overlapping cnv detected in Mother or Father, indicate it as a de novo cnv
		if INDEX_PED and (Father == "TRUE" or Mother == "TRUE"):
			DE_NOVO="FALSE"
		else:
			DE_NOVO="TRUE"
		
		##create new info field of vcf
		NEW_INFO=field[7]+";NB_INTERNALDB="+str(Numb_samples_match)+";DE_NOVO="+DE_NOVO+";Father="+Father+";Father_key="+Father_key+";Mother="+Mother+";Mother_key="+Mother_key
		
		## Write new line in vcf
		OUT.write('\t'.join(field[0:4])+'\t<'+ALT+'>\t'+'\t'.join(field[5:7])+'\t'+NEW_INFO+"\t"+"\t".join(field[8:]))
		
OUT.close()
l.close()
