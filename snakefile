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

SAMPLES, = glob_wildcards(datadir+"{sample}_1.fastq")
READS, =  glob_wildcards(datadir+"{readfile}.fastq")

NBSAMPLES = len(SAMPLES)
NBREADS = len(READS)

RUN = config["run"]
THREADS = config["threads"]

message(str(NBSAMPLES)+"samples  will be analysed")
message(str(len(READS))+"fastq files  will be processed")

for i in READS:
	message("read file: "+i)
 

REPAIRPATH = config["RepairPATH"]



message("Run name: "+RUN)



if NBREADS != 2*NBSAMPLES:
	errormes("Please provide two reads file per sample")
	sys.exit()


 
rule All:
	input:
		config["host"]+".bwt",
		expand("MappingOnAssembly/{smp}_pairs_on_"+RUN+".bam",smp=SAMPLES),
		"Megahit_results/"+RUN+".contigs.fa.bwt",
		expand("MappingOnAssembly/{smp}_SE_on_"+RUN+".bam",smp=SAMPLES),
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
		datadir+"{readfile}.fastq"
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
rule Repair_Pairs:
	input:
		R1="cutadaptfiles/{smp}_1.clean.fastq",
		R2="cutadaptfiles/{smp}_2.clean.fastq"
	params: exec = config["RepairPATH"]
	output:	
		"cutadaptfiles/{smp}_1.clean.fastq_pairs_R1.fastq",
		"cutadaptfiles/{smp}_2.clean.fastq_pairs_R2.fastq",
		"cutadaptfiles/{smp}_1.clean.fastq_widows.fastq"
	log: 	"logsRepairspairs/{smp}_repair.log"
	shell:
		"""
		python2 {params.exec} {input.R1} {input.R2} 2> log
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


rule Map_Pairs_On_Host:
	input:
		host = config["host"],
		hostindex = config["host"]+".bwt",
		R1="cutadaptfiles/{smp}_1.clean.fastq_pairs_R1.fastq",
		R2="cutadaptfiles/{smp}_2.clean.fastq_pairs_R2.fastq"
	output:
		"HostMapping/{smp}_pairs.bam"
	log:	
		"logsMapHost/{smp}_bwa_pairs_on_hist.log"
	
	threads: THREADS
	shell:		
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/samtools/1.3;
		module load bioinfo/bwa/0.7.15;
		bwa mem -t {threads} {input.host} {input.R1} {input.R2} 2> {log} | samtools view -b - > {output}
		"""

rule Map_Widows_On_Host:
	input:
		host = config["host"],
		hostindex = config["host"]+".bwt",
		Wi="cutadaptfiles/{smp}_1.clean.fastq_widows.fastq",
	output:
		"HostMapping/{smp}_widows.bam"
	log:	
		"logsMapHost/{smp}_bwa_widows_on_hist.log"
	threads: THREADS
	shell:		
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/samtools/1.3;
		module load bioinfo/bwa/0.7.15;
		bwa mem  -t {threads} {input.host} {input.Wi} 2> {log} | samtools view -b - > {output}
		"""

rule Extract_Unmapped_PE_Reads:
	input:
		"HostMapping/{smp}_pairs.bam"
	output:
		"HostMapping/unmapped_{smp}_pairs.bam"
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools view -b -hf 0x4 {input} > {output} 	
		"""	

rule Extract_Unmapped_Widows:
	input:
		"HostMapping/{smp}_widows.bam"
	output:
		"HostMapping/unmapped_{smp}_widows.bam"
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools view -b -hf 0x4 {input} > {output} 	
		"""	
rule Get_PE_Fastq_Filtered:
	input:
		"HostMapping/unmapped_{smp}_pairs.bam"
	output: 
		R1="FilteredFastq/filtered_{smp}_R1.fastq",
		R2="FilteredFastq/filtered_{smp}_R2.fastq"
	log:
		"logsGETPE/{smp}_getfilteredfastq_picard.log"	
	shell:
		"""
		java -jar /usr/local/bioinfo/picard-tools/1.130/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT I={input} F={output.R1} F2={output.R2} 2> {log}
		"""

rule Get_Widows_fastq_Filtered:
	input:
		"HostMapping/unmapped_{smp}_widows.bam"
	output: 
		"FilteredFastq/filtered_{smp}_widows.fastq"
	log:
		"logsGETWI/{smp}_getfiltered_widows_fastq_picard.log"	
	shell:
		"""
		java -jar /usr/local/bioinfo/picard-tools/1.130/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT I={input} F={output} 2> {log}
		"""

rule Merge_Pairs_With_Flash:
	input:
		R1="FilteredFastq/filtered_{smp}_R1.fastq",
		R2="FilteredFastq/filtered_{smp}_R2.fastq"
	params:
		prefix="FilteredFastq/{smp}"
	output:
		ext="FilteredFastq/{smp}.extendedFrags.fastq",
		R1="FilteredFastq/{smp}.notCombined_1.fastq",
		R2="FilteredFastq/{smp}.notCombined_2.fastq"
	log:
		"logsFLASH/{smp}_flash.log"
	shell:
		"""
		flash -M 250 {input.R1} {input.R2} -o {params.prefix} &> {log}	
		"""
rule Concatenate_Widows_And_Merged:
	input: 
		F="FilteredFastq/{smp}.extendedFrags.fastq",
		W="FilteredFastq/filtered_{smp}_widows.fastq"
	output:
		"FilteredFastq/filtered_{smp}_SE.fastq"
	shell:
		"""
		cat {input.F} {input.W} > {output}
		"""

# Files lists for assembly, needed to construct the string of comma separated arguments 
R1list=expand("FilteredFastq/{smp}.notCombined_1.fastq",smp=SAMPLES)			
R2list=expand("FilteredFastq/{smp}.notCombined_2.fastq",smp=SAMPLES)			
SElist=expand("FilteredFastq/filtered_{smp}_SE.fastq",smp=SAMPLES)   

rule Megahit:
	input:
		R1s = R1list,
		R2s = R2list,
		SEs = SElist
	params:
		prefix="Megahit_results",
		commaR1s = ",".join(R1list),
		commaR2s = ",".join(R2list),
		commaSEs = ",".join(SElist)
	output:
		"Megahit_results/{RUN}.contigs.fa"
	log:	"logsMegahit/Megahit_{RUN}.log"
	threads: THREADS
	shell:
		"""
		megahit -t {threads}  -m 180e9  -1 {params.commaR1s} -2 {params.commaR2s} -r {params.commaSEs} -o {params.prefix} --out-prefix {RUN} --continue  2> {log}
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

rule Map_PE_On_Assembly:
	input:
		R1="FilteredFastq/{smp}.notCombined_1.fastq",		
		R2="FilteredFastq/{smp}.notCombined_2.fastq",			
		CONTIGS = "Megahit_results/{RUN}.contigs.fa",
		INDEX =  "Megahit_results/{RUN}.contigs.fa.bwt"
	output:
		"MappingOnAssembly/{smp}_pairs_on_{RUN}.bam"
	threads:  THREADS
	shell:
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/bwa/0.7.15;
		bwa mem -t {threads} {input.CONTIGS} {input.R1} {input.R2} | samtools view -b - > {output}
		"""

rule Map_SE_On_Assembly:
	input:
		SE="FilteredFastq/filtered_{smp}_SE.fastq",
		INDEX =  "Megahit_results/{RUN}.contigs.fa.bwt",
		CONTIGS = "Megahit_results/{RUN}.contigs.fa"
	output:
		"MappingOnAssembly/{smp}_SE_on_{RUN}.bam"
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
PElist=expand("MappingOnAssembly/{smp}_pairs_on_{RUN}.bam",smp=SAMPLES,RUN=RUN)			

rule Quantify_contigs:
	input:
		SE = SElist,
		PE  = PElist
	output:
		"{RUN}_counts_contigs.mat"
	threads: 1
	shell:
		"""
		./count {input.SE} {input.PE} > {output}
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


