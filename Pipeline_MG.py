#!/usr/bin/env python
import time
import json
import re
import luigi
import os
import subprocess
import logging


#############################################################################################################################
class GlobalParameter(luigi.Config):
    human_ref = luigi.Parameter(default='/home/ubuntu/Data/human/hg18.fa')
    ref_genome= luigi.Parameter(default='/home/ubuntu/Data/Amphora/all.649.amphora.fst')
    basefolder= luigi.Parameter(default='/home/ubuntu/work/')
    amphorafolder = luigi.Parameter(default='/home/ubuntu/Data/Amphora/')
    cogfolder = luigi.Parameter(default='/home/ubuntu/Data/COG_annotation/')

def run_cmd(cmd):
    cmd = "t1=`date +%s`;"+cmd+";t2=`date +%s`;tdiff=`echo 'scale=3;('$t2'-'$t1')/60' | bc`;echo '##### Total time:  '$tdiff' mins'"
    p = subprocess.Popen(cmd,bufsize=-1, shell=True, universal_newlines=True, stdout=subprocess.PIPE,executable='/bin/bash')
    output = p.communicate()[0]
    return output


#############################################################################################################################
class TrimTrimmomatic(luigi.Task):
    sample = luigi.Parameter()
    n_cpu = luigi.Parameter(default="16")

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'p1':luigi.LocalTarget(folder+self.sample+"_P_1.fastq"),
                'u1':luigi.LocalTarget(folder+self.sample+"_U_1.fastq"),
                'p2':luigi.LocalTarget(folder+self.sample+"_P_2.fastq"),
                'u2':luigi.LocalTarget(folder+self.sample+"_U_2.fastq")}

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = "[ -d  {folder} ] || mkdir {folder};trimmomatic PE -threads {n_cpu}  {inputfolder}{sample}_1.fastq {inputfolder}{sample}_2.fastq {folder}{sample}_P_1.fastq {folder}{sample}_U_1.fastq {folder}{sample}_P_2.fastq {folder}{sample}_U_2.fastq LEADING:20 TRAILING:20 MINLEN:50".format(sample=self.sample,folder=folder,n_cpu=self.n_cpu,inputfolder=os.path.abspath(GlobalParameter().basefolder)+"/raw/")
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class DereplicateFastuniq(luigi.Task):
    sample = luigi.Parameter()
    def requires(self):
        return TrimTrimmomatic(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'p1':luigi.LocalTarget(folder+self.sample+"_P_1.fastq"),
                'p2':luigi.LocalTarget(folder+self.sample+"_P_2.fastq")}

    def run(self):
        folder=self.__class__.__name__+"/"
        cmd = ''' [ -d  {folder} ] || mkdir {folder};
        fastuniq -i <(printf "{p1}\n{p2}\n") -t q -o {p1m} -p {p2m} '''.format( p1=self.input()['p1'].path,p2=self.input()['p2'].path,p1m=self.output()['p1'].path,p2m=self.output()['p2'].path,folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class ContaimRemoveBmtagger(luigi.Task):
    sample = luigi.Parameter()
    human_ref = luigi.Parameter(default=GlobalParameter().human_ref) 
    def requires(self):
        return DereplicateFastuniq(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'p1m':luigi.LocalTarget(folder+self.sample+"_P_1.fastq"),
                'p2m':luigi.LocalTarget(folder+self.sample+"_P_2.fastq")}

    def run(self):
        folder=self.__class__.__name__+"/"
        cmd = ''' [ -d  {folder} ] || mkdir {folder};bmtagger.sh -b {human_ref}.w18.bitmask -x {human_ref}.srprism -q 1 -1 {p1} -2 {p2} -o {folder}{sample}_bmtagger.out
        cat {p1} |post_bmtagger_fqextract.py {folder}{sample}_bmtagger.out >{p1m}
        cat {p2} |post_bmtagger_fqextract.py {folder}{sample}_bmtagger.out >{p2m} '''.format(human_ref=self.human_ref, p1=self.input()['p1'].path,p2=self.input()['p2'].path,p1m=self.output()['p1m'].path,p2m=self.output()['p2m'].path,folder=folder,sample=self.sample)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class ContaimRemoveBwa(luigi.Task):
    sample = luigi.Parameter()
    human_ref = luigi.Parameter(default=GlobalParameter().human_ref) 
    n_cpu = luigi.Parameter(default="16")
    def requires(self):
        return DereplicateFastuniq(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'p1':luigi.LocalTarget(folder+self.sample+"_P_1.fastq"),
                'p2':luigi.LocalTarget(folder+self.sample+"_P_2.fastq")}

    def run(self):
        folder=self.__class__.__name__+"/"
        cmd = ''' [ -d  {folder} ] || mkdir {folder};
        bwa mem -t {n_cpu} {human_ref} {p1} {p2}|samtools fastq -f 13 -1 {p1m}.tmp  -2 {p2m}.tmp -
        sed '1~4 s/$/.1/g' {p1m}.tmp  >{p1m}
        sed '1~4 s/$/.2/g' {p2m}.tmp  >{p2m}
        rm -rf {p1m}.tmp {p2m}.tmp '''.format(human_ref=self.human_ref, p1=self.input()['p1'].path,p2=self.input()['p2'].path,p1m=self.output()['p1'].path,p2m=self.output()['p2'].path,folder=folder,sample=self.sample,n_cpu=self.n_cpu)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class PairReadsMergeFlash(luigi.Task):
    sample = luigi.Parameter()
    def requires(self):
        return ContaimRemoveBwa(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'p1':luigi.LocalTarget(folder+self.sample+".notCombined_1.fastq"),
                'p2':luigi.LocalTarget(folder+self.sample+".notCombined_2.fastq"),
                'm':luigi.LocalTarget(folder+self.sample+".extendedFrags.fastq")}

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = "[ -d  {folder} ] || mkdir {folder};flash  {p1} {p2} -d {folder} -o {sample}".format(p1=self.input()['p1'].path,p2=self.input()['p2'].path,folder=folder,sample=self.sample)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class AlignmentBowtie2(luigi.Task):
    sample = luigi.Parameter()
    ref_genome = luigi.Parameter(default=GlobalParameter().ref_genome)
    n_cpu = luigi.Parameter(default="16")
    def requires(self):
        return ContaimRemoveBwa(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'bam':luigi.LocalTarget(folder+self.sample+".bam"), 'sam':luigi.LocalTarget(folder+self.sample+".sam")}

    def run(self):
        folder= GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''
        [ -d  {folder} ] || mkdir {folder}
        bowtie2 --no-unal -p {n_cpu} -x {ref_genome}  -1 {p1} -2 {p2} -S {folder}{sample}.sam; 
        samtools view -hubS -@ {n_cpu} {folder}{sample}.sam |samtools sort -@ {n_cpu} 1>{folder}{sample}.bam '''.format(ref_genome=self.ref_genome,p1=self.input()['p1'].path,p2=self.input()['p2'].path,folder=folder,sample=self.sample,n_cpu=self.n_cpu)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class TaxonProfileKraken(luigi.Task):
    sample = luigi.Parameter()
    n_cpu = luigi.Parameter(default="16")
    def requires(self):
        return ContaimRemoveBwa(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'kraken':luigi.LocalTarget(folder+self.sample+".report")}

    def run(self):
	folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''
        [ -d  {folder} ] || mkdir -p {folder};
        kraken --threads {n_cpu}  --fastq-input --check-names --paired {p1} {p2} --output {folder}{sample}.kraken
        kraken-report {folder}{sample}.kraken > {folder}{sample}.report '''.format(p1=self.input()['p1'].path,p2=self.input()['p2'].path,sample=self.sample,folder=folder,n_cpu=self.n_cpu)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)


class AssemblySpades(luigi.Task):
    sample = luigi.Parameter()
    n_cpu = luigi.Parameter(default="16")
    def requires(self):
        return ContaimRemoveBwa(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"+self.sample+"/"
        return luigi.LocalTarget(folder+"scaffolds.fasta")

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = " [ -d  {folder} ] || mkdir -p {folder}".format(folder=folder)
        run_cmd(cmd)
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"+self.sample+"/"
        cmd = '''
        spades.py --meta -t {n_cpu} -1 {p1} -2 {p2} -o {folder}'''.format(
                p1=self.input()['p1'].path,
                p2=self.input()['p2'].path,
                n_cpu=self.n_cpu,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

#############################################################################################################################

class GeneFinderProdigal(luigi.Task):
    sample = luigi.Parameter()
    n_cpu = luigi.Parameter(default="16")
    assembler = luigi.Parameter(default="spades")

    def requires(self):
        if self.assembler == "velvet":
            return AssemblyVelvet(sample=self.sample) 
        elif self.assembler == "ray":
            return AssemblyRay(sample=self.sample) 
        else:
            return AssemblySpades(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'protein':luigi.LocalTarget(folder+self.sample+".genes.faa"),'transcript':luigi.LocalTarget(folder+self.sample+".genes.fna"),'gff':luigi.LocalTarget(folder+self.sample+".gff")} 

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''[ -d  {folder} ] || mkdir {folder};
        prodigal -q -m -p meta -a {folder}{sample}.genes.faa -d {folder}{sample}.genes.fna -f gff -o {folder}{sample}.gff -i {contigs}'''.format(
                contigs=self.input().path,
                sample=self.sample,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class CogAnnotation(luigi.Task):
    sample = luigi.Parameter()
    n_cpu = luigi.Parameter(default="16")
    cogfolder = luigi.Parameter(GlobalParameter().cogfolder)

    def requires(self):
        return GeneFinderProdigal(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'cogstats':luigi.LocalTarget(folder+self.sample+"/cog_stats.txt"), 'funcstats':luigi.LocalTarget(folder+self.sample+"/func_stats.txt"), 'proteinidcog':luigi.LocalTarget(folder+self.sample+"/protein-id_cog.txt"), 'rpsblastcog':luigi.LocalTarget(folder+self.sample+"/rps-blast_cog.txt") } 

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''[ -d  {folder} ] || mkdir {folder};
        cp {protein} {sample}
        fastasplit  -f {sample}.aa -o {folder} -c {n_cpu}
        parallel rpsblast -query {sample}.aa_chunk_0000{{}} -db {cogfolder}Cog  -out {folder}/{sample}.out_{{}} -evalue 1e-3 -outfmt 6 ::: $(eval echo {{000..{n_cpu_1}}})
        eval echo {folder}{sample}.out_{{000..{n_cpu_1}}} | xargs  cat  >{folder}{sample}.out
        cdd2cog.pl -r {folder}{sample}.out -c {cogfolder}/cddid.tbl -f {cogfolder}/fun.txt -w {cogfolder}/whog -o {folder}/{sample}/ '''.format(
                protein=self.input()['protein'].path,
                sample=self.sample,
                n_cpu=self.n_cpu,
                n_cpu_1=str(int(self.n_cpu)-1),
                cogfolder=self.cogfolder,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class CDSRefCogAnnotation(luigi.Task):
    samplelistfile = luigi.Parameter()
    sample="CDSRef"
    n_cpu = luigi.Parameter(default="16")
    cogfolder = luigi.Parameter(GlobalParameter().cogfolder)

    def requires(self):
        return CDSRefCreate(samplelistfile=self.samplelistfile)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'cogstats':luigi.LocalTarget(folder+self.sample+"/cog_stats.txt"), 'funcstats':luigi.LocalTarget(folder+self.sample+"/func_stats.txt"), 'proteinidcog':luigi.LocalTarget(folder+self.sample+"/protein-id_cog.txt"), 'rpsblastcog':luigi.LocalTarget(folder+self.sample+"/rps-blast_cog.txt") } 

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''[ -d  {folder} ] || mkdir {folder};
        cd {folder};
        cp {protein} {sample}.aa
        fastasplit  -f {sample}.aa -o {folder} -c {n_cpu}
        parallel rpsblast -query {sample}.aa_chunk_0000{{}} -db {cogfolder}Cog  -out {folder}/{sample}.out_{{}} -evalue 1e-3 -outfmt 6 ::: $(eval echo {{000..{n_cpu_1}}})
        eval echo {folder}{sample}.out_{{000..{n_cpu_1}}} | xargs  cat  >{folder}{sample}.out
        cdd2cog.pl -r {folder}{sample}.out -c {cogfolder}/cddid.tbl -f {cogfolder}/fun.txt -w {cogfolder}/whog -o {folder}/{sample}/ '''.format(
                protein=self.input()['aa'].path,
                sample=self.sample,
                n_cpu=self.n_cpu,
                n_cpu_1=str(int(self.n_cpu)-1),
                cogfolder=self.cogfolder,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class InterproscanAnnotation(luigi.Task):
    sample = luigi.Parameter()

    def requires(self):
        return GeneFinderProdigal(self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'tsv':luigi.LocalTarget(folder+self.sample+".tsv"), 'tsv':luigi.LocalTarget(folder+self.sample+".tsv"), 'gff3':luigi.LocalTarget(folder+self.sample+".gff3"), 'xml':luigi.LocalTarget(folder+self.sample+".xml")} 

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''[ -d  {folder} ] || mkdir {folder};
        sed -i 's/\*//g' {protein}
        interproscan.sh -i {protein} -f TSV,XML,GFF3,HTML,SVG -iprlookup -goterms -pa -b {folder}/{sample} '''.format(
                protein=self.input()['protein'].path,
                sample=self.sample,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

############################################################################################################################
class Sam2Amphora(luigi.Task):
    sample = luigi.Parameter()
    Pct = luigi.Parameter(default=90)
    amphorafolder = luigi.Parameter(GlobalParameter().amphorafolder)
    def requires(self):
        return AlignmentBowtie2(sample=self.sample)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'amphora_sample':luigi.LocalTarget(folder+self.sample+".txt")}

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''
        [ -d  {folder} ] || mkdir {folder};
        python {amphorafolder}/scripts/1.filter_sam.py {sam} {Pct} 40 |samtools view -hubS - |samtools sort > {folder}/{sample}.sort.bam;samtools index {folder}/{sample}.sort.bam
        perl {amphorafolder}/scripts/4.kpileup.generic.pl {sample} {folder}/{sample}.sort.bam  {amphorafolder}/amphora_gene_file.txt 0 0 10 > {folder}/{sample}.txt'''.format(
                amphorafolder=self.amphorafolder,
                sam=self.input()['sam'].path,
                Pct=self.Pct,
                sample=self.sample,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class AmphoraPerGene(luigi.Task):
    samplelistfile = luigi.Parameter()
    Pct = luigi.Parameter(default="90")
    amphorafolder = luigi.Parameter(GlobalParameter().amphorafolder)
    gene = luigi.Parameter()
    def requires(self):
        return [Sam2Amphora(sample=i,Pct=self.Pct,amphorafolder=self.amphorafolder) for i in [line.strip() for line in open(self.samplelistfile,"r")]]

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'amphora':luigi.LocalTarget(folder+self.gene+".amphora."+str(self.Pct)+".pickle")}

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''
        [ -d  {folder} ] || mkdir {folder};
        python {amphorafolder}/scripts/5.kp2np.py {Pct} {gene} {amphorafolder}/amphora_gene_file.txt {samplelistfile} {inputfolder} {folder}'''.format(
                amphorafolder=self.amphorafolder,
                inputfolder=os.path.dirname(self.input()[0]['amphora_sample'].path),
                Pct=self.Pct,
                samplelistfile = self.samplelistfile,
                gene=self.gene,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class Amphora2OTU(luigi.Task):
    samplelistfile = luigi.Parameter()
    Pct = luigi.Parameter(default=90)
    amphorafolder = luigi.Parameter(GlobalParameter().amphorafolder)
    genelistfile = luigi.Parameter(GlobalParameter().amphorafolder+"/bacterial_gene_list.txt")
    def requires(self):
        return [AmphoraPerGene(gene=i,samplelistfile=self.samplelistfile,Pct=self.Pct,amphorafolder=self.amphorafolder) for i in [line.strip() for line in open(self.genelistfile,"r")]]

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'amphora_otu':luigi.LocalTarget(folder+"cov.table.txt")}

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''
        [ -d  {folder} ] || mkdir {folder}; cd {inputfolder};
        for i in $(ls *.amphora.*.pickle);do python {amphorafolder}/scripts/sum_mean.py  $i ;done >{folder}cov.table.txt '''.format(
                amphorafolder=self.amphorafolder,
                inputfolder=os.path.dirname(self.input()[0]['amphora'].path),
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class StrainAmphorePerGene(luigi.Task):
    samplelistfile = luigi.Parameter()
    Pct = luigi.Parameter(default="90")
    amphorafolder = luigi.Parameter(GlobalParameter().amphorafolder)
    gene = luigi.Parameter()
    def requires(self):
        return AmphoraPerGene(samplelistfile=self.samplelistfile,Pct=self.Pct, amphorafolder=self.amphorafolder,gene=self.gene)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return [ luigi.LocalTarget(folder+self.gene+"_"+str(k)+".em") for k in range(2,10) ]

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = '''
        [ -d  {folder} ] || mkdir {folder};
        python {amphorafolder}/scripts/6.proc.py {out} {folder}/{gene}.out  {samplelistfile} 
        parallel "python {amphorafolder}/scripts/StrainFinder.py -N {{}} --aln {folder}/{gene}.out.aln.pickle --s_reps 1000 --s_iter 3 --d_reps 1 --d_iter 1000 --n_keep 3 --dtol 1 --ntol 2 --em_out {folder}{gene}_{{}}.em --log {folder}{gene}_{{}}.log --force_update --min_fdist 0.1 --min_gdist 0.01" ::: $(seq 2 10) '''.format(
                amphorafolder=self.amphorafolder,
                out=self.input()['amphora'].path,
                samplelistfile = self.samplelistfile,
                Pct=self.Pct,
                gene=self.gene,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class StrainAmphore(luigi.Task):
    samplelistfile = luigi.Parameter()
    Pct = luigi.Parameter(default="90")
    amphorafolder = luigi.Parameter(GlobalParameter().amphorafolder)
    genelistfile = luigi.Parameter(GlobalParameter().amphorafolder+"/bacterial_gene_list.txt")
    def requires(self):
        return [StrainAmphorePerGene(samplelistfile=self.samplelistfile,Pct=self.Pct, amphorafolder=self.amphorafolder,gene=i) for i in [line.strip() for line in open(self.genelistfile,"r")]] 

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return luigi.LocalTarget(folder+self.Pct+"SF_results.pkl")

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd =''' 
        [ -d  {folder} ] || mkdir {folder};
        python {amphorafolder}/scripts/get_strain_estimates.py {genelistfile} {folder}{Pct}. {inputfolder} {Pct}'''.format(
                amphorafolder=self.amphorafolder,
                genelistfile=self.genelistfile,
                inputfolder= os.path.dirname(self.input()[0][0].path),
                Pct=self.Pct,
                folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

############################################################################################################################

class TaxonProfileKrakenList(luigi.Task):
    samplelistfile = luigi.Parameter()
    def requires(self):
        return [ TaxonProfileKraken(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]

    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))

class ContaimRemoveBwaList(luigi.Task):
    samplelistfile = luigi.Parameter()
    def requires(self):
        return [ContaimRemoveBwa(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]


    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))

class GeneFinderProdigalList(luigi.Task):
    samplelistfile = luigi.Parameter()
    def requires(self):
        return [GeneFinderProdigal(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]


    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())

        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))

class CogAnnotationList(luigi.Task):
    samplelistfile = luigi.Parameter()
    def requires(self):
        return [ CogAnnotation(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]


    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))

class InterproscanAnnotationList(luigi.Task):
    samplelistfile = luigi.Parameter()
    def requires(self):
        return [ InterproscanAnnotation(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]

    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))

class PairReadsMergeFlashList(luigi.Task):
    samplelistfile = luigi.Parameter()
    def requires(self):
        return [ PairReadsMergeFlash(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]

    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))

class DereplicateFastuniqList(luigi.Task):
    samplelistfile = luigi.Parameter()
    def requires(self):
        return [ DereplicateFastuniq(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]

    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))


class TrimTrimmomaticList(luigi.Task):
    samplelistfile = luigi.Parameter()
    def requires(self):
        return [ TrimTrimmomatic(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]

    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))


class CDSRefCreate(luigi.Task):
    samplelistfile = luigi.Parameter()

    def requires(self):
        return [ GeneFinderProdigal(sample=i) for i in [line.strip() for line in open(self.samplelistfile,"r")]]

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return {'nt':luigi.LocalTarget(folder+"/cds.nr.nt"), 'aa':luigi.LocalTarget(folder+"/cds.nr.aa")}

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = ''' [ -d  {folder} ] || mkdir {folder};
        for i in {nt_list};do tmp=$(basename $i); tmps=${{tmp/.genes.fna}};sed  "s/>/>${{tmps}}_/g" $i;done > {folder}/cds.nt ;
        for i in {nt_list};do tmp=$(basename $i); tmps=${{tmp/.genes.fna}};sed  "s/>/>${{tmps}}_/g" ${{i/.genes.fna}}.genes.faa;done > {folder}/cds.aa;
        cd-hit-est -d 0 -n 10 -l 100 -p 1 -G 0 -c 0.95 -aS 0.8 -M 0 -T 0 -i {folder}/cds.nt -o {folder}/cds.nr.nt
        awk '/^>/{{print "seq" ++i, $0; next}}' {folder}/cds.nr.nt|sed 's/ >/\t/g' >{folder}/rename.table
        awk '/^>/{{print $0; next}}' {folder}/cds.nr.nt|sed 's/ >/\t/g' >{folder}/name.list
        awk '/^>/{{print ">seq" ++i; next}}{{print}}' {folder}/cds.nr.nt >tmp;mv tmp {folder}/cds.nr.nt 
        seqtk subseq {folder}/cds.aa  {folder}/name.list >{folder}/cds.aa.nr
        extract_rename_seq.py {folder}/cds.aa.nr {folder}/name.list {folder}/cds.nr.aa
        bowtie2-build {folder}/cds.nr.nt {folder}/cds.nr.nt '''.format(nt_list=" ".join([ i['transcript'].path for i in self.input()]),folder=folder)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class AbundanceEstimate(luigi.Task):
    sample = luigi.Parameter()
    ref_genome = luigi.Parameter()

    def requires(self):
        return AlignmentBowtie2(sample=self.sample,ref_genome=self.ref_genome)

    def output(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        return luigi.LocalTarget('{folder}/{sample}.table'.format(folder=folder,sample=self.sample))

    def run(self):
        folder=GlobalParameter().basefolder+"/"+self.__class__.__name__+"/"
        cmd = ''' [ -d  {folder} ] || mkdir {folder};
        gen_contig_cov_per_bam_table.py {ref_genome} {bam} >{folder}/{sample}.table '''.format(ref_genome=self.ref_genome,bam=self.input()['bam'].path,folder=folder,sample=self.sample)
        print("****** NOW RUNNING COMMAND ******: " + cmd)
        print run_cmd(cmd)

class AbundanceEstimateList(luigi.Task):
    samplelistfile = luigi.Parameter()
    ref_genome = luigi.Parameter(default="CDSRefCreate/cds.nr.nt")

    def requires(self):
        return [ AbundanceEstimate(sample=i,ref_genome=self.ref_genome) for i in [line.strip() for line in open(self.samplelistfile,"r")]]

    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('workflow finished at {t}'.format(t=timestamp))

#########################################################################################################################
if __name__=='__main__':
    luigi.run()
