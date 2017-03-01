# CMIT Metagenomic Pipeline Common Use Standards

## Launch an instance of the Pipeline
* Launch the instance from AMI `almlab_cluster_09162016 (ami-fcbbcdeb)`. This will provide the environment for the pipeline. Remember to leave port 8888 open.
* Create a volume from Snapshot `snap-8332db98`. The volume will contain the tools that required for the pipeline to run. Mount this volume to the folder `/home/ubuntu/tools`.

## Specify the absolute path to your data
The pipeline also requires the path to certain data in order to perform certain task. You should specify the path to the data in file `luigi.cfg`

```{bash}
[core]
default-scheduler-port:8888 

[GlobalParameter] # Remember that all the path should be absolute
ref_genome=/home/ubuntu/Data/Amphora/all.649.amphora.fst #path to amphora gene fasta file
human_ref=/home/ubuntu/Data/human/hg18.fa # path to human genomes (bowtie2 indexed)
amphorafolder=/home/ubuntu/Data/Amphora/ # path to Amphora folder
basefolder=/home/ubuntu/work # the absolute path to your working folder where you want to run the pipeline and where your raw fastq data is stored.
```
You can download all this data from `s3://almlab.bucket/xiaofang/Luigi_Workflow_DB`

If you need to run Kraken, specify the path to your database by run command `export KRAKEN_DEFAULT_DB=path_to_kraken_db`. 

## Prepare your data for the pipeline
1. Create a folder as your working folder (Specified by `basefolder` in your `luigi.cfg` file )
2. Put the raw metagenomic fastq file in folder `Raw` inside your working folder
3. Keep your raw fastq files in the format `samplename_1.fastq` and  `samplename_2.fastq` (paired end).
4. Create a file `sample.list` containg the names of all samples, each sample in one line without any spaces.
5. Run command line `luigid --background --port 8888`.
6. Put `Pipeline_MG.py`, `luigi.cfg` in your working folder.

## Run your task
Task examples:
### The following command line will run three tasks subsequently for each sample
+ TrimTrimmomatic: trim low QC fastq files
+ DereplicateFastuniq: remove replicates from fastq files
+ ContaimRemoveBwa: remove human contamination

```{bash}
python Pipeline_MG.py  ContaimRemoveBwaList --samplelistfile sample.list --workers 2 1>ContaimRemoveBwaList.log 2>ContaimRemoveBwaList.err 
```
The results will be in the folder: `TrimTrimmomatic` ,`DereplicateFastuniq` and `ContaimRemoveBwa` of your working folder.

### Run Kraken
```{bash}
python Pipeline_MG.py  TaxonProfileKrakenList --samplelistfile sample.list --workers 2 1>TaxonProfileKrakenList.log 2>TaxonProfileKrakenList.err 
```
This step will automatically run task `TrimTrimmomatic` ,`DereplicateFastuniq` and `ContaimRemoveBwa` if you have not done so. 

### Run funcational annotation
```{bash}
python Pipeline_MG.py  CDSRefCogAnnotation --samplelistfile sample.list --workers 2 1>CDSRefCogAnnotation.log 2>CDSRefCogAnnotation.err 
python Pipeline_MG.py  AbundanceEstimateList --samplelistfile sample.list --workers 2 1>AbundanceEstimateList.log 2>AbundanceEstimateList.err 
```
This will take a long time.

1. CDSRefCogAnnotation: 
  * Assemble and predict protein coding genes for each sample
  * Combine all protein coding genes and create a non-redundant fasta file
  * Annotate the non-redundant fasta file with COG terms
2. AbundanceEstimateList: (Run `CDSRefCogAnnotation` task first before run this step)
  * Align each sample to the non-redundant fasta file to estimate the abundance of genes in each sample

## Tips:
1. The pipeline is based on the spotify workflow [luigi](https://github.com/spotify/luigi). You can modify the pipeline to your desire based on the [document](http://luigi.readthedocs.io/en/stable/workflows.html).
2. Choose the number of workers based on the number of cores in your instance. Most tasks in the pipeline use 16 cores. So if your instance has 48 cores, your worker should be `48/16=3`.
3. You can check your task status and get better understanding of dependence of individual task from browser by entering `http://youripaddress:8888`.
