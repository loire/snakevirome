import re
import sys
from os.path import join

def message(mes):
	sys.stderr.write("|---- " + mes + "\n")

def errormes(mes):
	sys.stderr.write("| ERROR ----" + mes + "\n")

configfile : "config.yaml"

datadir = config["fastq"]
scriptdir = config["Scripts"]
KronaDir  = config["KronaDir"]
basse_virale = config["basse_virale"]
basse_nt = config["basse_nt"]
SAMPLES, = glob_wildcards(datadir+"{sample}_1.fastq.gz")
READS, =  glob_wildcards(datadir+"{readfile}.fastq.gz")

NBSAMPLES = len(SAMPLES)
NBREADS = len(READS)

RUN = config["run"]
THREADS = config["threads"]

message(str(NBSAMPLES)+" samples  will be analysed")
message(str(len(READS))+" fastq files  will be processed")
message("Run name: "+RUN)

if NBREADS != 2*NBSAMPLES:
	errormes("Please provide two reads file per sample")
	sys.exit()


rule All:
	input:
		config["host"]+".bwt",
		expand("cutadaptfiles/{smp}_1.clean.fastq_pairs_R1.fastq",smp=SAMPLES),
		expand("cutadaptfiles/{smp}_2.clean.fastq_pairs_R2.fastq",smp=SAMPLES),
		expand("cutadaptfiles/{smp}_1.clean.fastq_widows.fastq",smp=SAMPLES),
		expand("logs_contaminent/Stats_mapping_contaminent_widows_{smp}.txt",smp=SAMPLES),
		"Megahit_results/"+RUN+"_Assembly_results",
		"Megahit_results/"+RUN+"_.contigs.cap.fa",
		"Megahit_results/"+RUN+"_.contigs.cap.fa.bwt",
		"logs/"+RUN+"mega_assembly_inf.txt",
		"logs/"+RUN+"cap_assembly_inf.txt",
		expand("MappingOnAssembly/{smp}_R_1_2_on_"+RUN+".bam",smp=SAMPLES),
		expand("MappingOnAssembly/{smp}_SE_on_"+RUN+".bam",smp=SAMPLES),
		expand("MappingOnAssembly/{smp}_PE_on_"+RUN+".bam",smp=SAMPLES),
		expand("logs/{smp}_"+RUN+"_stats_mapping_assembly.txt",smp=SAMPLES),
		expand("Unmapped/{smp}_unmapped_"+RUN+".fastq",smp=SAMPLES),
		expand("Unmapped/{smp}_PE_unmapped_"+RUN+".fastq",smp=SAMPLES),
		expand("Unmapped/{smp}_SE_unmapped_"+RUN+".fastq",smp=SAMPLES),
		expand("Unmapped/{smp}_PE_unmapped_"+RUN+".fa",smp=SAMPLES),
		expand("Unmapped/{smp}_SE_unmapped_"+RUN+".fa",smp=SAMPLES),
		expand("CountsMapping/{smp}_SE_counts_contigs_"+RUN+".mat",smp=SAMPLES),
		expand("CountsMapping/{smp}_PE_counts_contigs_"+RUN+".mat",smp=SAMPLES),
		expand("CountsMapping_all/{smp}_R_1_2_counts_contigs_"+RUN+".mat",smp=SAMPLES),
		expand("CountsMapping_all/{smp}_R_1_2_counts_unmapped_"+RUN+".mat",smp=SAMPLES),
		"Blast_nt_vir_results/contigs_"+RUN+".blastnresults.tsv",
		expand("Blast_nt_vir_results/{smp}_unmapped_PE_"+RUN+".blastnresults.tsv",smp=SAMPLES),
		expand("Blast_nt_vir_results/{smp}_unmapped_SE_"+RUN+".blastnresults.tsv",smp=SAMPLES),
		expand("logs/{smp}_"+RUN+"_blast_ntvir.txt",smp=SAMPLES),
		"Blast_ntv_hits_seq/hits_contigs_"+RUN+".fa",
		expand("Blast_ntv_hits_seq/{smp}_PE_hits_unmapped_"+RUN+".fa",smp=SAMPLES),
		expand("Blast_ntv_hits_seq/{smp}_SE_hits_unmapped_"+RUN+".fa",smp=SAMPLES),
		"Blast_nt_results/"+RUN+".blast_nt_nresults.tsv",
		expand("Blast_nt_results/{smp}_unmapped_PE_"+RUN+".blastnresults.tsv",smp=SAMPLES),
		expand("Blast_nt_results/{smp}_unmapped_SE_"+RUN+".blastnresults.tsv",smp=SAMPLES),
		dynamic("tmp/id_for_tax_{num}"),
		"Taxonomy/"+RUN+"_Seq_hits_info.csv",
		"Coverage/"+RUN+"_coverage_contigs_by_sample.csv",
		"Coverage/"+RUN+"_coverage_lineage_by_sample.csv",
		expand("logs/{smp}_report.md",smp=SAMPLES),
		"logs/report_"+RUN+"_Assembly.md",
		expand("logs/{smp}_results.html",smp=SAMPLES),


rule Remove_sequencing_adapters:
	input:
		datadir+"{readfile}.fastq.gz",
	output:
		"cutadaptfiles/{readfile}.trimmed.fastq",
	log:
		"logscutadapt/{readfile}_cut1.log",
	benchmark:
		"benchmarks/{readfile}.cut1.benchmark.txt"

	params:
		A3 = config["A3"],
		A5 = config["A5"],
	shell:
		"""
		module load bioinfo/cutadapt/1.8.1 ;
		cutadapt  -n 10 -g {params.A5} -a {params.A3} --overlap 15 {input} -o {output} &> {log}
		"""

rule Quality_trimming:
	input:
		"cutadaptfiles/{readfile}.trimmed.fastq"
	output:
		"cutadaptfiles/{readfile}.clean.fastq"
	log:
		"logscutadapt/{readfile}_cut2.log"
	benchmark:
		"benchmarks/{readfile}.cut2.benchmark.txt"
	shell:
		"""
		module load bioinfo/cutadapt/1.8.1 ;
		cutadapt -q 30,30 -m 40 -o {output} {input} &> {log}
		"""

rule Repair_Pairs:
	input:
		R1="cutadaptfiles/{smp}_1.clean.fastq",
		R2="cutadaptfiles/{smp}_2.clean.fastq"
	output:
		"cutadaptfiles/{smp}_1.clean.fastq_pairs_R1.fastq",
		"cutadaptfiles/{smp}_2.clean.fastq_pairs_R2.fastq",
		"cutadaptfiles/{smp}_1.clean.fastq_widows.fastq",
	params:
		scriptdir+"/fastqCombinePairedEnd.py"
	log:
		"logsRepairspairs/{smp}_repair.log"
	benchmark:
		"benchmarks/{smp}.cut2.benchmark.txt"
	shell:
		"""
		python2 {params} {input.R1} {input.R2} 2> log
		"""

###### All this part is on host ribosomal dna sequence ######

###### I think you can add bacterial sequences and filters the BAM files to get whatever you want from it ####
###### To discuss with serafin ######

rule Index_Host_Sequences:
	input:
		config["host"]
	output:
		config["host"]+".bwt"
	benchmark:
		"benchmarks/index_host_seq.benchmark.txt"
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
	benchmark:
		"benchmarks/{smp}.Map_Pairs_On_Host.benchmark.txt"
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
	benchmark:
		"benchmarks/{smp}.Map_Widows_On_Host.benchmark.txt"
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
		"HostMapping/unmapped_{smp}_pairs.bam",
	benchmark:
		"benchmarks/{smp}.Extract_Unmapped_PE_Reads.benchmark.txt"
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
	benchmark:
		"benchmarks/{smp}.Extract_Unmapped_Widows.benchmark.txt"
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools view -b -hf 0x4 {input} > {output}
		"""

rule log_map_host:
	input:
		p="HostMapping/{smp}_pairs.bam",
		w="HostMapping/{smp}_widows.bam",
	output:
		p="logs_contaminent/Stats_mapping_contaminent_pair_{smp}.txt",
		w="logs_contaminent/Stats_mapping_contaminent_widows_{smp}.txt",
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools flagstat {input.p} > {output.p}
		samtools flagstat {input.w} > {output.w}
		"""

rule Get_PE_Fastq_Filtered:
	input:
		"HostMapping/unmapped_{smp}_pairs.bam"
	output:
		R1="FilteredFastq/filtered_{smp}_R1.fastq",
		R2="FilteredFastq/filtered_{smp}_R2.fastq"
	benchmark:
		"benchmarks/{smp}.Get_PE_Fastq_Filtered.benchmark.txt"
	log:
		"logsGETPE/{smp}_getfilteredfastq_picard.log"
	shell:
		"""
		module load system/java/jre8
		java -jar /usr/local/bioinfo/picard-tools/1.130/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT I={input} F={output.R1} F2={output.R2} 2> {log}
		"""

rule Get_Widows_fastq_Filtered:
	input:
		"HostMapping/unmapped_{smp}_widows.bam"
	output:
		"FilteredFastq/filtered_{smp}_widows.fastq"
	benchmark:
		"benchmarks/{smp}.Get_Widows_fastq_Filtered.benchmark.txt"
	log:
		"logsGETWI/{smp}_getfiltered_widows_fastq_picard.log"
	shell:
		"""
		module load system/java/jre8
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
	benchmark:
		"benchmarks/{smp}.Merge_Pairs_With_Flash.benchmark.txt"
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
		"FilteredFastq/filtered_{smp}_PE.fastq"
	benchmark:
		"benchmarks/{smp}.Concatenate_Widows_And_Merged.benchmark.txt"
	shell:
		"""
		cat {input.F} {input.W} > {output}
		"""


# Files lists for assembly, needed to construct the string of comma separated arguments
R1list=expand("FilteredFastq/{smp}.notCombined_1.fastq",smp=SAMPLES)
R2list=expand("FilteredFastq/{smp}.notCombined_2.fastq",smp=SAMPLES)
PElist=expand("FilteredFastq/filtered_{smp}_PE.fastq",smp=SAMPLES)

rule Megahit_Assembly:
	input:
		R1s = R1list,
		R2s = R2list,
		PEs = PElist
	params:
		prefix="Megahit_results",
		commaR1s = ",".join(R1list),
		commaR2s = ",".join(R2list),
		commaPEs = ",".join(PElist)
	output:
		"Megahit_results/{RUN}.contigs.fa"
	log:
		"logsMegahit/Megahit_{RUN}.log"
	benchmark:
		"benchmarks/Megahit_Assembly.benchmark.txt"
	threads: THREADS
	shell:
		"""
		megahit -t {threads}  -m 180e9  -1 {params.commaR1s} -2 {params.commaR2s} -r {params.commaPEs} -o {params.prefix} --out-prefix {RUN} --continue  2> {log}
		"""


rule Cap3_Assembly:
	input:
		"Megahit_results/{RUN}.contigs.fa"
	output:
		ass="Megahit_results/{RUN}_Assembly_results",
		cont="Megahit_results/{RUN}.contigs.fa.cap.contigs",
		sig="Megahit_results/{RUN}.contigs.fa.cap.singlets",
	threads: THREADS
	benchmark:
		"benchmarks/Cap3_Assembly.benchmark.txt"
	shell:
		"""
		module load bioinfo/CAP3/20150317
		cap3 {input}>{output.ass}
		"""

rule Merge_contigs:
	input:
		cont="Megahit_results/{RUN}.contigs.fa.cap.contigs",
		sig="Megahit_results/{RUN}.contigs.fa.cap.singlets",
	benchmark:
		"benchmarks/Merge_contigs.benchmark.txt"
	output:
		"Megahit_results/{RUN}_.contigs.cap.fa"
	shell:
		"""
		cat {input.cont} {input.sig} > {output}
		"""

rule Assembly_informations:
	input:
		Mega="Megahit_results/{RUN}.contigs.fa",
		Cap="Megahit_results/{RUN}_.contigs.cap.fa",
	output:
		Mega="logs/{RUN}mega_assembly_inf.txt",
		Cap="logs/{RUN}cap_assembly_inf.txt",

	shell:
		"""
		cat {input.Mega} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}'> {output.Mega}
		cat {input.Cap} | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}'> {output.Cap}
		"""

rule Index_Assembly:
	input:
		"Megahit_results/{RUN}_.contigs.cap.fa"
	output:
		"Megahit_results/{RUN}_.contigs.cap.fa.bwt"
	benchmark:
		"benchmarks/Index_Assembly.benchmark.txt"
	shell:
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/bwa/0.7.15;
		bwa index {input} > {output}
		"""

rule Map_SE_On_Assembly:
	input:
		R1= "FilteredFastq/{smp}.notCombined_1.fastq",
		R2= "FilteredFastq/{smp}.notCombined_2.fastq",
		CONTIGS = "Megahit_results/{RUN}_.contigs.cap.fa",
		INDEX =  "Megahit_results/{RUN}_.contigs.cap.fa.bwt"
	output:
		"MappingOnAssembly/{smp}_SE_on_{RUN}.bam"
	benchmark:
		"benchmarks/Map_SE_On_Assembly.benchmark.txt"
	threads:  THREADS
	shell:
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/bwa/0.7.15;
		module load bioinfo/samtools/1.3;
		bwa mem -t {threads} {input.CONTIGS} {input.R1} {input.R2} | samtools view -b -> {output}
		"""

rule Map_PE_On_Assembly:
	input:
		PE="FilteredFastq/filtered_{smp}_PE.fastq",
		CONTIGS = "Megahit_results/{RUN}_.contigs.cap.fa",
		INDEX =  "Megahit_results/{RUN}_.contigs.cap.fa.bwt"
	output:
		"MappingOnAssembly/{smp}_PE_on_{RUN}.bam"
	benchmark:
		"benchmarks/Map_PE_On_Assembly.benchmark.txt"
	threads:	THREADS
	shell:
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/bwa/0.7.15;
		module load bioinfo/samtools/1.3;
		bwa mem -t {threads} {input.CONTIGS} {input.PE} | samtools view -b - > {output}
		"""

rule Map_On_Assembly:
	input:
		R1="cutadaptfiles/{smp}_1.clean.fastq",
		R2="cutadaptfiles/{smp}_2.clean.fastq",
		CONTIGS = "Megahit_results/{RUN}_.contigs.cap.fa",
		INDEX =  "Megahit_results/{RUN}_.contigs.cap.fa.bwt"
	output:
		R1="MappingOnAssembly/{smp}_R1_on_{RUN}.bam",
		R2="MappingOnAssembly/{smp}_R2_on_{RUN}.bam",
		R="MappingOnAssembly/{smp}_R_1_2_on_{RUN}.bam"
	benchmark:
		"benchmarks/Map_ALL_On_Assembly.benchmark.txt"

	threads:  THREADS
	shell:
		"""
		module load compiler/gcc/4.9.2;
		module load bioinfo/bwa/0.7.15;
		module load bioinfo/samtools/1.3;
		bwa mem -t {threads} {input.CONTIGS} {input.R1} | samtools view -b -> {output.R1};
		bwa mem -t {threads} {input.CONTIGS} {input.R2} | samtools view -b -> {output.R2};
		samtools merge {output.R} {output.R1} {output.R2}
		"""


rule Extract_All_Umapped_on_contigs:
	input:
		"MappingOnAssembly/{smp}_R_1_2_on_{RUN}.bam",
	output:
		"Unmapped/{smp}_unmapped_{RUN}.fastq",
	benchmark:
		"benchmarks/Extract_All_Umapped_on_contigs.benchmark.txt"

	threads:  THREADS
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools view -b -hf 0x4 {input} | samtools bam2fq -> {output};
		"""

rule Mapping_information:
	input:
		"MappingOnAssembly/{smp}_R_1_2_on_{RUN}.bam",
	output:
		"logs/{smp}_{RUN}_stats_mapping_assembly.txt",
	threads: 1
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools flagstat {input} > {output};
		"""

rule Extract_filtered_Umapped_on_contigs:
	input:
		PE="MappingOnAssembly/{smp}_PE_on_{RUN}.bam",
		SE="MappingOnAssembly/{smp}_SE_on_{RUN}.bam",
	output:
		PE="Unmapped/{smp}_PE_unmapped_{RUN}.fastq",
		SE="Unmapped/{smp}_SE_unmapped_{RUN}.fastq",
	benchmark:
		"benchmarks/Extract_filtered_Umapped_on_contigs.benchmark.txt"
	threads: 1
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools view -b -hf 0x4 {input.PE} | samtools bam2fq -> {output.PE};
		samtools view -b -hf 0x4 {input.SE} | samtools bam2fq -> {output.SE};
		"""

rule Unmapped_fq_to_fa:
	input:
		PE="Unmapped/{smp}_PE_unmapped_{RUN}.fastq",
		SE="Unmapped/{smp}_SE_unmapped_{RUN}.fastq",
	output:
		PE="Unmapped/{smp}_PE_unmapped_{RUN}.fa",
		SE="Unmapped/{smp}_SE_unmapped_{RUN}.fa",
	benchmark:
		"benchmarks/Unmapped_fq_to_fa.benchmark.txt"

	threads: 1
	shell:
		"""
		seqtk seq -A {input.PE} > {output.PE};
		seqtk seq -A {input.SE} > {output.SE};
		"""


rule Quantify_contigs_coverage:
	input:
		SE = "MappingOnAssembly/{smp}_SE_on_{RUN}.bam",
		PE = "MappingOnAssembly/{smp}_PE_on_{RUN}.bam",
		R_1_2= "MappingOnAssembly/{smp}_R_1_2_on_{RUN}.bam"
	output:
		SE= "CountsMapping/{smp}_SE_counts_contigs_{RUN}.mat",
		PE= "CountsMapping/{smp}_PE_counts_contigs_{RUN}.mat",
		R_1_2= "CountsMapping_all/{smp}_R_1_2_counts_contigs_{RUN}.mat",
		Unmapped="CountsMapping_all/{smp}_R_1_2_counts_unmapped_{RUN}.mat"
	benchmark:
		"benchmarks/Quantify_contigs_coverage.benchmark.txt"
	threads: 1
	shell:
		"""
		module load bioinfo/samtools/1.3;
		samtools view -F 0x4 {input.SE}| cut -f 3 | sort | uniq -c - > {output.SE};
		samtools view -F 0x4 {input.PE}| cut -f 3 | sort | uniq -c - > {output.PE};
		samtools view -F 0x4 {input.R_1_2}| cut -f 3 | sort | uniq -c - > {output.R_1_2};
		samtools view -f 0x4 {input.R_1_2}| cut -f 1 | sort | uniq -c - > {output.Unmapped};
		"""


rule Blast_contigs_on_nt_vir:
	input:
		"Megahit_results/{RUN}_.contigs.cap.fa",

	output:
		"Blast_nt_vir_results/contigs_{RUN}.blastnresults.tsv",
	params:
		blastDBpath=basse_virale
	benchmark:
		"benchmarks/Blast_contigs_on_nt_vir.benchmark.txt"
	threads: 1
	shell:
		"""
		module load bioinfo/ncbi-blast/2.2.30;
		blastn -task blastn -db {params.blastDBpath} -query {input} -num_threads {threads} -evalue 0.001 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output};
		"""



rule Blast_unmapped_on_nt_vir:
	input:
		PE= "Unmapped/{smp}_PE_unmapped_{RUN}.fa",
		SE= "Unmapped/{smp}_SE_unmapped_{RUN}.fa",

	output:
		BLAST_PE= "Blast_nt_vir_results/{smp}_unmapped_PE_{RUN}.blastnresults.tsv",
		BLAST_SE= "Blast_nt_vir_results/{smp}_unmapped_SE_{RUN}.blastnresults.tsv",
	params:
		blastDBpath=basse_virale
	benchmark:
		"benchmarks/Blast_unmapped_on_nt_vir.benchmark.txt"
	threads: 1
	shell:
		"""
		module load bioinfo/ncbi-blast/2.2.30;
	 	blastn -task blastn -db {params.blastDBpath} -query {input.PE} -num_threads {threads} -evalue 0.001 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output.BLAST_PE};
	  	blastn -task blastn -db {params.blastDBpath} -query {input.SE} -num_threads {threads} -evalue 0.001 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output.BLAST_SE};
		"""

rule Extract_hits_on_nt_vir_contigs_fa:
	input:
		BLAST= "Blast_nt_vir_results/contigs_{RUN}.blastnresults.tsv",
		CONTIGS= "Megahit_results/{RUN}_.contigs.cap.fa",
	output:
		"Blast_ntv_hits_seq/hits_contigs_{RUN}.fa"
	params:
		scriptdir+"extract_blast_hits.py"
	benchmark:
		"benchmarks/Extract_hits_on_nt_vir_contigs_fa.benchmark.txt"
	shell:
		"""
		python {params} {input.BLAST} {input.CONTIGS} {output}
		"""

rule Extract_hits_on_nt_vir_unmapped_fa:
	input:
		BLAST_PE= "Blast_nt_vir_results/{smp}_unmapped_PE_{RUN}.blastnresults.tsv",
		BLAST_SE= "Blast_nt_vir_results/{smp}_unmapped_SE_{RUN}.blastnresults.tsv",
		PE= "Unmapped/{smp}_PE_unmapped_{RUN}.fa",
		SE= "Unmapped/{smp}_SE_unmapped_{RUN}.fa"
	output:
		SEQ_PE="Blast_ntv_hits_seq/{smp}_PE_hits_unmapped_{RUN}.fa",
		SEQ_SE="Blast_ntv_hits_seq/{smp}_SE_hits_unmapped_{RUN}.fa",
	params:
		scriptdir+"extract_blast_hits.py"
	benchmark:
		"benchmarks/Extract_hits_on_nt_vir_unmapped_fa.benchmark.txt"
	shell:
		"""
		python {params} {input.BLAST_PE} {input.PE} {output.SEQ_PE};
		python {params} {input.BLAST_SE} {input.SE} {output.SEQ_SE};
		"""
rule log_ntvir:
	input:
		SEQ_PE="Unmapped/{smp}_PE_unmapped_{RUN}.fa",
		SEQ_SE="Unmapped/{smp}_SE_unmapped_{RUN}.fa",
		contigs="Megahit_results/{RUN}_.contigs.cap.fa",
		BLAST_PE="Blast_nt_vir_results/{smp}_unmapped_PE_{RUN}.blastnresults.tsv",
		BLAST_SE="Blast_nt_vir_results/{smp}_unmapped_SE_{RUN}.blastnresults.tsv",
		BLAST_contigs="Blast_nt_vir_results/contigs_{RUN}.blastnresults.tsv",
	output:
		"logs/{smp}_{RUN}_blast_ntvir.txt"
	shell:
		"""
		grep -c '>' {input.SEQ_SE} {input.SEQ_PE} >> {output}
		wc -l {input.BLAST_SE} {input.BLAST_PE} >> {output}
		grep -c '>' {input.contigs} >> {output}
		wc -l {input.BLAST_contigs} >> {output}
		"""


rule Blast_contigs_on_nt:
	input:
		"Blast_ntv_hits_seq/hits_contigs_{RUN}.fa",
	output:
		"Blast_nt_results/{RUN}.blast_nt_nresults.tsv"
	params:
		blastDBpath = basse_nt
	benchmark:
		"benchmarks/Blast_contigs_on_nt.benchmark.txt"
	threads: 1
	shell:
		"""
		module load bioinfo/ncbi-blast/2.2.30;
		blastn -task blastn -db {params.blastDBpath} -query {input} -num_threads {threads} -evalue 0.001  -max_hsps 1 -max_target_seqs 10 -outfmt 6  -out {output}
		"""


rule Blast_unmapped_on_nt:
	input:
		SEQ_PE="Blast_ntv_hits_seq/{smp}_PE_hits_unmapped_{RUN}.fa",
		SEQ_SE="Blast_ntv_hits_seq/{smp}_SE_hits_unmapped_{RUN}.fa",
	output:
		BLAST_PE="Blast_nt_results/{smp}_unmapped_PE_{RUN}.blastnresults.tsv",
		BLAST_SE="Blast_nt_results/{smp}_unmapped_SE_{RUN}.blastnresults.tsv",
	params:
		blastDBpath = basse_nt
	benchmark:
		"benchmarks/Blast_unmapped_on_nt.benchmark.txt"
	threads: 1
	shell:
		"""
		module load bioinfo/ncbi-blast/2.2.30;
		blastn -task blastn -db {params.blastDBpath} -query {input.SEQ_PE} -num_threads {threads} -evalue 0.001 -max_hsps 1 -max_target_seqs 10 -outfmt 6 -out {output.BLAST_PE}
		blastn -task blastn -db {params.blastDBpath} -query {input.SEQ_SE} -num_threads {threads} -evalue 0.001 -max_hsps 1 -max_target_seqs 10 -outfmt 6 -out {output.BLAST_SE}
		"""



rule Join_best_hits:
	input: 
	output:	"tmp/join_hits",
	benchmark:
		"benchmarks/Join_best_hits.benchmark.txt"
	shell:
		"""
		cat Blast_nt_results/* > {output}
		"""

rule extract_seq_acc:
	input:"tmp/join_hits",
	output: "tmp/seq_acc_ids.txt",
	benchmark:
		"benchmarks/extract_seq_acc.benchmark.txt"
	shell:
		"""
		awk  '{{print $1, $2 ,$3, $11, $12}}' {input} | awk -F'|' '{{print $1, $2, $3, $4, $5}}'| awk  '{{print $1 ,substr($5,1, length($5)-2) , $(NF-2),$(NF-1),$(NF)}}' >> {output}
		"""

rule extract_acc_and_split:
	input:
		"tmp/seq_acc_ids.txt",
	output:
		dynamic("tmp/id_for_tax_{num}"),
	benchmark:
		"benchmarks/extract_acc_and_split.benchmark.txt",
	shell:
		"""
		awk '{{print $2}}' {input} | split -d -l 1000 - {output}
		"""

rule extract_tax_ids:
	input:
		taxo_dmp= KronaDir+"all.accession2taxid.sorted",
		acc_ids= "tmp/id_for_tax_{num}"
	output:
		"tmp/tax_ids_{num}",
	benchmark:
		"benchmarks/extract_tax_ids.benchmark.txt"
	shell:
		"""
		for i in `cat {input.acc_ids}` ; do  look $i {input.taxo_dmp}  >> {output}; done || true
		"""

rule join_acc_tax:
	input:  tax_ids=dynamic("tmp/tax_ids_{num}"),
	output: tax_ids="tmp/acc_tax_id.txt"
	benchmark:
		"benchmarks/join_acc_tax.benchmark.txt"
	shell:
		"""
		cat {input.tax_ids}  > {output.tax_ids}
		"""

rule Join_seq_acc_taxo :
	input: tax_ids="tmp/acc_tax_id.txt", ids="tmp/seq_acc_ids.txt"
	output: "Taxonomy/{RUN}_Seq_hits_info.csv",
	benchmark:
		"benchmarks/Join_seq_acc_taxo.benchmark.txt"
	shell:
		"""
		awk 'NR==FNR {{h[$1] = $2; next}} {{print $1,$2,$3,$4,$5,h[$2]}}' {input.tax_ids} {input.ids} | awk 'NF==6' | sort -k1,1 -k5,5gr -k4,4g -k3,3gr> {output}
		"""

rule get_lineage_from_taxids:
	input: "Taxonomy/{RUN}_Seq_hits_info.csv"
	output: "Taxonomy/{RUN}_lineage.csv"
	params: scriptdir+"get_rank.py"
	benchmark:
		"benchmarks/get_lineage_from_taxids.benchmark.txt"
	shell:
		"""
		export PATH=~/anaconda_ete/bin:$PATH
		awk '{{print $6}}' < {input} |sort -u | paste -s -d, | tr -d '\n'| python {params} {output}
		"""


rule Build_array_coverage:
	input:
		blast_info="Taxonomy/Aedesdataset_Seq_hits_info.csv",
		lineage="Taxonomy/Aedesdataset__lineage.csv",
	output:
		contigs_coverage="Coverage/{RUN}_coverage_contigs_by_sample.csv",
		lineage_coverage="Coverage/{RUN}_coverage_lineage_by_sample.csv"
	params:
		scriptdir+"build_tables_coverage.py",
		countdir = "CountsMapping_all/"
	benchmark:
		"benchmarks/Build_array_coverage.benchmark.txt"
	shell:
		"""
		python {params} {input.blast_info} {params.countdir} {input.lineage} {output.contigs_coverage} {output.lineage_coverage}
		"""

rule Create_logs_report:
	input:
	output:
		"logs/{smp}_report.md"
	params:
		files=datadir+"*_1.fastq.gz",
		script=scriptdir+"create_results_doc.py"
	benchmark:
		"benchmark:"
	benchmark:
		"benchmarks/Create_logs_report.benchmark.txt"
	shell:
		"""
		for entry in {params.files} ; do parts=(${{entry//_1./ }}); parts1=(${{parts//// }}) ; s=(${{parts1[5]}});python {params.script} ${{s}} ; done
		"""

rule Create_assembly_report:
	input:
		Mega="logs/{RUN}mega_assembly_inf.txt",
		Cap="logs/{RUN}cap_assembly_inf.txt",
	output:
		"logs/report_{RUN}_Assembly.md"
	params:
		script=scriptdir+"create_results_assembly_doc.py",
	benchmark:
		"benchmarks/Create_assembly_report.benchmark.txt"
	shell:
		"""
		python {params.script} {input.Mega} {input.Cap} {output}
		"""

rule Report_html:
	input:
		"logs/{smp}_report.md"
	output:
		"logs/{smp}_results.html"
	params:
		scriptdir+"github-pandoc.css"
	benchmark:
		"benchmarks/Report_html.benchmark.txt"
	shell:
		"""
		pandoc -f markdown_github -t html5 -c {params} -o {output} {input}
		"""
