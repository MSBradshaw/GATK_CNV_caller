#!usr/bin/python # 

import subprocess
import sys
import re
# import MySQLdb // for python2
# from MySQLdb import cursors
import mysql.connector
import urllib
import gzip

OMIM='/d2/Exome_analysis_BC/Exome_pipeline_datas/OMIM/OMIM_gene_patho.txt'
refseq='/d2/Exome_analysis_BC/Exome_pipeline_datas/humandb/hg19_refGene.txt'

### Path Canoes_file ###
VCF_file=snakemake.input[0]
OUT_file=snakemake.output[0]
#VCF_file=sys.argv[1]
#OUT_file=sys.argv[2]
IN=gzip.open(VCF_file,'rt')
OUT=open(OUT_file,'w')

### Dictionary from OMIM file [gene:disease] ###
OMIM=open(OMIM,'r')
OMIM_lib={}
for line in OMIM:
    field=line.split('\t')
    if field[0] in OMIM_lib.keys():
        OMIM_lib[field[0]].append(field[1][:-1])
    else:
        OMIM_lib[field[0]]=[field[1][:-1]]
OMIM.close()

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
            OUT.write("SAMPLE\tSVTYPE\tKEY\tSIZE\tGENE\tpLI\tPATHO\tINTERNAL_DB\tVARIANTS_(hom|het)"+'\t'+'AD_DUP'+'\t'+'RARE_VARIANTS'+'\t'+'FREQ RARE VARIANTS(<0.005)'+'\t'+"CN\tNP\tQA\tQS\tQSE\tQSS"+"\n")
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
        ## Keep only cnvs found in less than 3 individuals
        if int(internal_db) > 3:
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
        exon=[]
        refseq_file=open(refseq,'r')
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
                        genes.append(name_gene)
        if not OMIM_match:
            OMIM_disease="."
        refseq_file.close()    
   
    ### Variants from database annotation ###             
        ### Get Internal database datas ###
        # Database connection #
        try:
            db = mysql.connector.connect(host="127.0.0.1",user="bcogne",passwd="bencog15",database='Exome_db')
            db2 = mysql.connector.connect(host="127.0.0.1",user="bcogne",passwd="bencog15",database='External_db')
            
            cur_gnomad = db2.cursor(buffered=True)
            pLI=[]
            for i in genes:
                command = "SELECT pLI FROM gnomad_constraint WHERE gene='"+i+"';"
                cur_gnomad.execute(command)
                results = cur_gnomad.fetchall()
                for row in results:
                    pLI.append(str(round(row[0],2)))
            cur_gnomad.close()
            
            cur = db.cursor(buffered=True)
            command = "SELECT chr,pos,ref,alt,gene,loc,exon,zyg,AD,exonic_function FROM "+name+" WHERE loc='exonic' AND gene in (%s);" %('"'+'","'.join(genes)+'"')
            cur.execute(command)
            results = cur.fetchall() #... WHERE gene in ('gene1','gene2','..)
            nb_variants_hom = 0
            nb_variants_het = 0

            variants_rares=[] #chr/pos/ref/alt/loc/freq
            freq_list = []
            gene_list = []
            #OMIM_list = []
            nb_variants_list = []
            AD_ratio = []
           
            if 'DEL' in SVTYPE:
                       
                ### If query results exists ###
                if results:
                    for row in results:
                    ### Get names to each column selected : chr,pos,ref,... ###
                        chr = row[0]
                        pos = row[1]
                        ref = row[2]
                        alt = row[3]
                        gene = row[4]
                        loc = row[5]
                        zyg = row[7]
                        AD = row[8]
                        function = row[9]
                        if zyg == 'het':
                            nb_variants_het+=1
                        else:
                            nb_variants_hom+=1

                        
                        cur_exac = db2.cursor(buffered=True)
                        ### Get internal database frequency and choose frequency < 0.005 ###
                        command_exac = "SELECT freq FROM exac3 WHERE chr='"+str(chr)+"' AND pos="+str(pos)+" AND ref='"+str(ref)+"' AND alt='"+str(alt)+"';"
                        cur_exac.execute(command_exac)
                        if not cur_exac.rowcount:
                            freq_list.append('.')
                            variants_rares.append(str(gene)+'/'+str(chr)+'/'+str(pos)+'/'+str(ref)+'/'+str(alt)+'/'+str(zyg)+'/'+str(function))
                        else:
                            results_exac = cur_exac.fetchall()
                            for row in results_exac:
                                freq_exac = float(row[0])
                                if freq_exac < 0.005:
                                    freq_list.append(str(freq_exac))
                                    variants_rares.append(str(gene)+'/'+str(chr)+'/'+str(pos)+'/'+str(ref)+'/'+str(alt)+'/'+str(zyg)+'/'+str(function))
         
                        #cur_exac.close()
                    nb_variants_list.append(str(nb_variants_hom)+'|'+str(nb_variants_het))    
                  
            if 'DUP' in SVTYPE:
               
                ### If query results exists ###
                if results:
                    for row in results:
                    ### Get names to each column selected : chr,pos,ref,... ###
                        chr = row[0]
                        pos = row[1]
                        ref = row[2]
                        alt = row[3]
                        gene = row[4]
                        loc = row[5]
                        zyg = row[7]
                        AD = row[8]
                        function = row[9]
                        cur_exac = db2.cursor(buffered=True)
                        ### Get internal database frequency and choose frequency < 0.005 ###
                        command_exac = "SELECT freq FROM exac3 WHERE chr='"+str(chr)+"' AND pos="+str(pos)+" AND ref='"+str(ref)+"' AND alt='"+str(alt)+"';"
                        cur_exac.execute(command_exac)
                        if not cur_exac.rowcount:
                            freq_list.append('.')
                            variants_rares.append(str(gene)+'/'+str(chr)+'/'+str(pos)+'/'+str(ref)+'/'+str(alt)+'/'+str(zyg)+'/'+str(function))
                        
                        else:
                            results_exac = cur_exac.fetchall()
                            for row in results_exac:
                                freq_exac = float(row[0])
                                if freq_exac < 0.005:
                                    freq_list.append(str(freq_exac))
                                    variants_rares.append(str(gene)+'/'+str(chr)+'/'+str(pos)+'/'+str(ref)+'/'+str(alt)+'/'+str(zyg)+'/'+str(function))   
                        cur_exac.close()
                        AD_ratio.append(str(AD))  
            ### Writing to OUT file ###
            OUT.write(name+"\t"+SVTYPE+"\t"+cnv_key+"\t"+str(cnv_size)+"\t"+';'.join(genes)+'\t'+";".join(pLI)+"\t"+';'.join(OMIM_disease)+'\t'+internal_db+'\t'+''.join(nb_variants_list)+'\t'+'|'.join(AD_ratio)+'\t'+';'.join(variants_rares)+'\t'+';'.join(freq_list)+'\t'+CN+'\t'+NP+'\t'+QA+'\t'+QS+'\t'+QSE+'\t'+QSS+"\n")
         
            db2.commit()
            db.commit()    
            # Disconnection from server #
            cur.close() 
            
        except:
            exit()
    
OUT.close() 
IN.close() 

# ### Dictionary from OUT file [gene:count gene] ### ### To replace/change if a database is created ###
# ### count_gene count times a gene is observed in population ###
# OUT_final=open(Canoes_file+'anno-final.tab','w')

# OUT=open(Canoes_file+'anno.tab','r')
# count_gene_lib={}
# for line in OUT:
    # count_gene=[]
    # field=line.split('\t')
    # field_gene=field[4].split(';')
    # for gene in field_gene:
        # if gene in count_gene_lib.keys():
            # count_gene_lib[gene]+=1
        # else:
            # count_gene_lib[gene]=1 
    # if field_gene==[]:
        # count_gene[field_gene]=''


# OUT.close()
# OUT=open(Canoes_file+'anno.tab','r')
# b=0

# for line in OUT:
    # b+=1
    # field=line.split('\t')
    # if b==1:
            # OUT_final.write('\t'.join(field[0:6])+'\t'+'FREQ GENE'+'\t'+'\t'.join(field[7:]))
            # pass
    # else:  
        # count_gene=[]
        # gene_list=field[4].split(';')    
        # for gene in gene_list:
            # if gene in count_gene_lib.keys():
                # count_gene.append(str(count_gene_lib[gene]))
        # OUT_final.write('\t'.join(field[0:6])+'\t'+';'.join(count_gene)+'\t'+'\t'.join(field[7:]))

# OUT.close()
# OUT_final.close()

