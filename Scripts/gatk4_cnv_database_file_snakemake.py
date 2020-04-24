#!/bin/python

import sys
import os
import gzip


files_vcf=snakemake.input
database_file=snakemake.params[0]
new_database=snakemake.output[0]

## Create disctionary from database file if exist
print("Create a dictionary from the database")

dic={}
names_dic={}
if not os.path.exists(database_file):
    print("No database file detected")
else:
    database=open(database_file,'r')
    for line in database:
        if line[0] != "#":
            field=line.split('\t')
            chr_d=field[0]
            start_d=field[1]
            end_d=field[2]
            cnv_d=field[3]
            number_d=int(field[4].rstrip())
            names=field[5].rstrip().split(",")            
            dic[chr_d+'-'+start_d+'-'+end_d+'-'+cnv_d]=number_d
            names_dic[chr_d+'-'+start_d+'-'+end_d+'-'+cnv_d]=names
    database.close()

## Iterate over vcf files

for file in files_vcf:
    print("Analysis of file ",file)
    
    if file.endswith('vcf.gz'):
        l=gzip.open(file,'rt')
    elif file.endswith('vcf'):
        l=open(file,'r')
    else:
        print(file,"hasn't a compatible extension: vcf or vcf.gz")
        exit()
                
    for line in l:
        if line[0:6] == "#CHROM":
            name=line.split('\t')[9].rstrip()
        if line[0] != "#":
            field=line.split('\t')
            chr=field[0]
            start=field[1]
            end=field[7].split('=')[1]
            info=field[8].split(':')
            index_GT=info.index("GT")
            data=field[9]
            data_GT=data[index_GT]
            
            if data_GT == "0":
                continue
            elif data_GT == "1":
                cnv="DEL"
            elif data_GT == "2":
                cnv="DUP"
            
            key=chr+'-'+start+'-'+end+'-'+cnv
            
            if key in dic.keys():
                if name not in names_dic.keys():
                    dic[key]+=1
                    names_dic[key].append(name)
                else:
                    continue
            else:
                dic[key]=1
                names_dic[key]=[name]
    l.close()

#Write new database
print("Write new database file")

database=open(new_database,'w')
#database.write("#CHR\tSTART\tEND\tCNV\tNUMBER\tNAMES\n")
for i in dic.keys():
    field=i.split('-')
    chr=field[0]
    start=field[1]
    end=field[2]
    cnv=field[3]
    number=str(dic[i])
    names=",".join(names_dic[i])
    database.write(chr+"\t"+start+"\t"+end+"\t"+cnv+"\t"+number+"\t"+names+"\n")
database.close()
