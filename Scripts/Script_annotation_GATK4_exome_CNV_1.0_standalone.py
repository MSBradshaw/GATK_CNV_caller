#!usr/bin/python # 

import subprocess
import sys
import re
import urllib
import gzip

### Get paths of files from snakemake###
VCF_file=snakemake.input[0]
OMIM=snakemake.input[1]
refseq=snakemake.input[2]
gnomAD_file=snakemake.input[3]
Max_freq=snakemake.input[4]
database=snakemake.input[5]
OUT_file=snakemake.output[0]

IN=gzip.open(VCF_file,'rt')
OUT=open(OUT_file,'w')

# Get Number of samples
DATABASE=open(database,'r')
for line in DATABASE:
    if line[0]=="#":
        field=line.split('\t')
        number_of_samples=int(field[5].rstrip().split('_')[0])
    else:
        break
DATABASE.close()

### Dictionary from OMIM file [gene:disease] ###
OMIM=gzip.open(OMIM,'rt')
OMIM_lib={}
for line in OMIM:
    field=line.split('\t')
    if field[0] in OMIM_lib.keys():
        OMIM_lib[field[0]].append(field[1][:-1])
    else:
        OMIM_lib[field[0]]=[field[1][:-1]]
OMIM.close()

### Dictionary from gnomAD constraint file ###
gnomAD=gzip.open(gnomAD_file,'rt')
gnomAD_lib={}
a=0
for line in gnomAD:
    a+=1
    if a == 1:
        continue
    field=line.split('\t')
    if field[0] in gnomAD_lib.keys():
        if float(field[1].rstrip()) > float(gnomAD_lib[field[0]]):
            gnomAD_lib[field[0]]=field[1].rstrip()
    else:
        gnomAD_lib[field[0]]=field[1].rstrip()
gnomAD.close()

gene_disease=''
a=0
cnv_del_count=0
cnv_dup_count=0

    
### Iterate over VCF file
for line in IN:
    a+=1
    if line[0] == "#":
        if line[0:6] == "#CHROM":
            field=line.split('\t')
            name=field[9].rstrip()
            OUT.write("SAMPLE\tSVTYPE\tKEY\tSIZE\tGENE\tpLI\tPATHO\tINTERNAL_DB\tCN\tNP\tQA\tQS\tQSE\tQSS"+"\n")
    else:
        field=line.split('\t')

        chr_cnv=field[0]
        cnv_start=int(field[1])
        SVTYPE=field[4]

        INFO=field[7].split(";")
        for i in INFO:
            if i.split("=")[0] == "END":
                cnv_end=int(i.split("=")[1])
            if i.split("=")[0] == "NB_INTERNALDB":
                internal_db=i.split("=")[1]
        cnv_size=cnv_end-cnv_start
        cnv_key=chr_cnv+":"+str(cnv_start)+"-"+str(cnv_end)
        
        ## Keep only cnvs found in less than the defined max frequency
        if float(internal_db)/float(number_of_samples) > float(Max_freq):
            continue
        
        format=field[8].split(":")
        for i in format:
            if i == "CN":
                CN_index=format.index(i)
            if i == "NP":
                NP_index=format.index(i)
            if i == "QA":
                QA_index=format.index(i)
            if i == "QS":
                QS_index=format.index(i)
            if i == "QSE":
                QSE_index=format.index(i)
            if i == "QSS":
                QSS_index=format.index(i)

        sample=field[9].split(":")
        CN=sample[CN_index]
        NP=sample[NP_index]
        QA=sample[QA_index]
        QS=sample[QS_index]
        QSE=sample[QSE_index]
        QSS=sample[QSS_index].rstrip()
              
        ### Lists type to int type 
        genes=[]
        OMIM_disease=[]
        OMIM_match=False
        gnomAD_pLI=[]
        gnomAD_match=False
        exon=[]
        refseq_file=gzip.open(refseq,'rt')
        for line in refseq_file:
            field_refseq=line.split('\t')
            chr_refseq=str(field_refseq[2])
            txstart=int(field_refseq[4])
            txend=int(field_refseq[5])
            name_gene=field_refseq[12]
            exon_start=field_refseq[9].split(',')
            exon_end=field_refseq[10].split(',')
            strand=str(field_refseq[3])
            if chr_refseq == chr_cnv:
                if (txstart > cnv_start and txend < cnv_end) or (txstart < cnv_start and cnv_start < txend and txend < cnv_end) or (txstart > cnv_start and txend > cnv_end and txstart < cnv_end) or (cnv_start > txstart and cnv_end < txend):                      
                    if not name_gene in genes:
                        if name_gene in OMIM_lib.keys():
                            OMIM_disease.append(','.join(OMIM_lib[name_gene]))
                            OMIM_match=True
                        if name_gene not in OMIM_lib.keys():
                            OMIM_disease.append('.')
                        
                        if name_gene in gnomAD_lib.keys():
                            gnomAD_pLI.append(float(gnomAD_lib[name_gene]))
                            gnomAD_match=True                     

                        genes.append(name_gene)
        if not OMIM_match:
            OMIM_disease="."
        if not gnomAD_match:
            gnomAD_pLI=-1
        elif len(gnomAD_pLI) > 1:
            gnomAD_pLI=round(max(gnomAD_pLI),2)
        else:
            gnomAD_pLI=round(gnomAD_pLI[0],2)
        refseq_file.close()    
   
        ### Writing to OUT file ###
        OUT.write(name+"\t"+SVTYPE+"\t"+cnv_key+"\t"+str(cnv_size)+"\t"+';'.join(genes)+'\t'+str(gnomAD_pLI)+'\t'+';'.join(OMIM_disease)+'\t'+internal_db+'\t'+CN+'\t'+NP+'\t'+QA+'\t'+QS+'\t'+QSE+'\t'+QSS+"\n")
         
    
OUT.close() 
IN.close() 

