# SV-exome

This is a pipeline to analyze CNVs from a cohort of bam or cram files with GATK4.

It runs the COHORT mode. It means that you need a batch of exomes, at least 30 sequenced in the same conditions sounds reasonable for beginning, but you can try with less.

The pipeline contains several steps:
1) Count of the number of reads per intervals
2) The annotation and the filtration of intervals with gc content and extreme read counts
3) The calculation of the ploidy for all autosomal and sex chromosomes, it allows to call CNVs on chrX while reducing sex biases
4) The creation of a model and the calling of the CNVs for each samples provided
5) The generation of segmented VCFs
6) The creation of a file counting each CNVs found in your samples ("database")
7) BETA (to be updated): The creation of a tabulated file with annotated informations: frequency in your cohort, refseq genes impacted, OMIM

## Prerequisite
Have an activated conda environment with a recent version of snakemake

## Installation process

1) Download zip or clone the git repo
    
    ~~~
    git clone https://gitlab.univ-nantes.fr/benjamin.co/sv-exome.git
    cd sv-exome
    ~~~

2) Unzip the data reference files 
    ~~~
    gunzip ref/*.gz
    ~~~

3) Edit the environment file "envs/gatkcondaenv.yml": Add the absolute path of the last line: /path/to/sv-exome/envs/gatkPythonPackageArchive.zip

## Launch of the pipeline

1) First, edit the configuration file **conf/config.yaml**

2) Then launch pipeline with snakemake:
    ~~~
    snakemake -s Snakefile_exome -j /Number of threads you wish/ --use-conda --configfile conf/config.yaml --resources cnv_caller=1
    ~~~
    The parameter "--resources cnv_caller=1" is not required but allows to specifically attributes the number of germline cnv calling jobs run in parallel as it is really consuming a lot of RAM, which increase proportionally to the number of samples processed. Exemple: for 150 exomes, if you have more than 100Go RAM and at least 20 CPUs available (default value of 5 CPUs per job in the config file), you can consider using the maximum: cnv_caller=4 (default value if you choose -j options of at least 4CPUs)

## DEBUG

If during conda environment installation you have the following message:

~~~
ResolvePackageNotFound:
  - certifi==2016.2.28=py36_0
  - xz==5.2.3=0
  - readline==6.2=2
  - openssl==1.0.2l=0
  - tk==8.5.18=0
  - pip==9.0.1=py36_1
  - python==3.6.2=0
  - zlib==1.2.11=0
  - sqlite==3.13.0=0
  - setuptools==36.4.0=py36_1
  - wheel==0.29.0=py36_0
~~~

Just do before running snakemake: 
~~~
conda config --set restore_free_channel true
~~~