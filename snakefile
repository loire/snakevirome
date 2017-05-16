import glob
import re
import sys
from os.path import join 

def message(mes):
	sys.stderr.write("|---- " + mes + "\n")

def errormes(mes):
	sys.stderr.write("| ERROR ----" + mes + "\n")

configfile : "config.yaml"

datadir = config["fastq"] 

SAMPLES, = glob_wildcards(datadir+"{sample}.fastq.gz")
READS, =  glob_wildcards(datadir+"{readfile}.fastq.gz")

NBSAMPLES = len(SAMPLES)
NBREADS = len(READS)

RUN = config["run"]
THREADS = config["threads"]

message(str(NBSAMPLES)+"samples  will be analysed")
message(str(len(READS))+"fastq files  will be processed")
message("Run name: "+RUN)


 
rule All:
	input:
		config["host"]+".bwt",
		expand("MappingOnAssembly/{smp}_on_"+RUN+".bam",smp=SAMPLES),
		"Megahit_results/"+RUN+".contigs.fa.bwt",
		"Blast_results/"+RUN+".blastnresults.tsv",
		RUN+".taxo.txt",
		RUN+"_counts_contigs.mat",
		RUN+"_DE_list.txt"
##		RUN+"_resultat.mk"
	shell:
		"""
		mv snakejob* snakejoblog/
		"""

rule Cutadapt_1:
	input:
		datadir+"{readfile}.fastq.gz"
	output:
		"cutadaptfiles/{readfile}.trimmed.fastq"
	log:
		"logscutadapt/{readfile}_cut1.log"
	params:
		A3 = config["A3"],
		A5 = config["A5"]
	shell:
		"""
		module load bioinfo/cutadapt/default/ ;
		cutadapt  -n 10 -g {params.A5} -a {params.A3} --overlap 15 {input} -o {output} &> {log}
		"""

rule Cutadapt_2:
	input:
		"cutadaptfiles/{readfile}.trimmed.fastq"	
	output:
		"cutadaptfiles/{readfile}.clean.fastq"	
	log:
		"logscutadapt/{readfile}_cut2.log"
	shell:
		"""
		module load bioinfo/cutadapt/default/ ;
		cutadapt -q 30,30 -m 40 -o {output} {input} &> {log}
		"""

####### All this part is on host ribosomal dna sequence ######

####### I think you can add bacterial sequences and filters the BAM files to get whatever you want from it ####
####### To discuss with serafin ######

rule Index_Host_Sequences:
	input:
		config["host"]
	output:
		config["host"]+".bwt"
	shell:
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/bwa/0.7.15;
		bwa index {input}
		"""

rule Map_Reads_On_Host:
	input:
		host = config["host"],
		hostindex = config["host"]+".bwt",
		Wi = "cutadaptfiles/{readfile}.clean.fastq"	
	output:
		"HostMapping/{smp}.bam"
	log:	
		"logsMapHost/{smp}_bwa__on_host.log"
	threads: THREADS
	shell:		
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/samtools/1.3;
		module load bioinfo/bwa/0.7.15;
		bwa mem  -t {threads} {input.host} {input.Wi} 2> {log} | samtools view -b - > {output}
		"""
rule Extract_Unmapped:
	input:
		"HostMapping/{smp}.bam"
	output:
		"HostMapping/unmapped_{smp}.bam"
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools view -b -hf 0x4 {input} > {output} 	
		"""	
rule Get_Widows_fastq_Filtered:
	input:
		"HostMapping/unmapped_{smp}.bam"
	output: 
		"FilteredFastq/filtered_{smp}.fastq"
	log:
		"logsGETWI/{smp}_getfiltered_fastq_picard.log"	
	shell:
		"""
		java -jar /usr/local/bioinfo/picard-tools/1.130/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT I={input} F={output} 2> {log}
		"""
# Files lists for assembly, needed to construct the string of comma separated arguments 
SElist=expand("FilteredFastq/filtered_{smp}.fastq",smp=SAMPLES)   

rule Megahit:
	input:
		SEs = SElist
	params:
		prefix="Megahit_results",
		commaSEs = ",".join(SElist)
	output:
		"Megahit_results/{RUN}.contigs.fa"
	log:	"logsMegahit/Megahit_{RUN}.log"
	threads: THREADS
	shell:
		"""
		megahit -t {threads}  -r {params.commaSEs} -o {params.prefix} --out-prefix {RUN} --continue  2> {log}
		"""	

rule Index_Assembly:
	input:
		"Megahit_results/{RUN}.contigs.fa"
	output:
		"Megahit_results/{RUN}.contigs.fa.bwt"
	shell:	
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/bwa/0.7.15;
		bwa index {input}
		"""
rule Map_SE_On_Assembly:
	input:
		SE="FilteredFastq/filtered_{smp}.fastq",
		INDEX =  "Megahit_results/{RUN}.contigs.fa.bwt",
		CONTIGS = "Megahit_results/{RUN}.contigs.fa"
	output:
		"MappingOnAssembly/{smp}_on_{RUN}.bam"
	threads:	THREADS
	shell:
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/bwa/0.7.15;
		bwa mem -t {threads} {input.CONTIGS} {input.SE} | samtools view -b - > {output}
		"""

rule Blast_contigs:
	input: 
		CONTIGS = "Megahit_results/{RUN}.contigs.fa",
		VIRAL_DATA_BASE = "ntvir"
	output:
		"Blast_results/{RUN}.blastnresults.tsv"
	threads:	THREADS
	shell:
		"""
		blastn {input.CONTIGS} {input.VIRAL_DATA_BASE} > {output}
		"""

SElist=expand("MappingOnAssembly/{smp}_SE_on_{RUN}.bam",smp=SAMPLES,RUN=RUN)			

rule Quantify_contigs:
	input:
		SE = SElist,
	output:
		"{RUN}_counts_contigs.mat"
	threads: 1
	shell:
		"""
		./count {input.SE} > {output}
		"""


rule DE_expressed_contigs:
	input:
		Sample_desc = "Sample_desc.mat",
		counts = "{RUN}_counts_contigs.mat"
	output:
		"{RUN}_DE_list.txt"
	threads: 1
	shell:
		"""
		./count {input.Sample_desc} {input.counts} > {output}
		"""


rule Taxonomic_assignation:
	input:
		BLASTRESULTS="Blast_results/{RUN}.blastnresults.tsv",
		taxdb  = "taxo.dump"
	output:
		"{RUN}.taxo.txt"
	threads: 1
	shell:
		"""
		./count {input.BLASTRESULTS} {input.taxdb} > {output}
		""



#
#
#
#rule Pandoc:
#	input:
#		taxo = "{RUN}.taxo.txt",
#		DE = "{RUN}_DE_list.txt",
#		counts = "{RUN}_counts_contigs.mat"
#	output:
#		"{RUN}_resultat.mk"
#	threads: 1
#	shell:
#		"""
#		./count {input.SE} {input.PE} > {output}
#		""
#

#
#rule Taxonomic_assignation:
#	input:
#		BLASTRESULTS="Blast_results/{RUN}.blastnresults.tsv"
#	ouput:
#		"{RUN}.taxo.txt"
#	threads: 1
#	shell:
#		"""
#		"""
#

##
#
#threads: 1 
#	shell:
#		"""
#		kronatools {input.BLASTRESULT} {input.TAXODB} > {output}
#		"""
##
#rule DE_analysis:
#	input:
#		SE = "MappingOnAssembly/{smp}_pairs_on_{RUN}.bam",
#		PE = "MappingOnAssembly/{smp}_SE_on_{RUN}.bam",
#		#TAXO = "{RUN}.taxo.txt"
#	output:
#		"{RUN}_DE_contigs.tsv"
#	threads: 1
#	shell:
#		"""
#		EdgeR {input.SE} {input.PE} > {output}
#		"""
#

rule clean:
	shell:
		"""
		rm -rf FilteredFastq HostMapping Megahit_results logs cutadaptfiles MappingOnAssembly snakejob*
		"""


