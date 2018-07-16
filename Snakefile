GENOMES, = glob_wildcards("fasta/{genome}.fna")
# divisions are: Q1<2.5Kp; 2.5Kbp<Q2<5Kbp; 5Kbp<Q3<10Kbp; 10Kbp<Q4
LENGTH = ["Q1","Q2"]
COVERAGE = [5,10]
SIM_NANO_READS = 2000

rule all:
	input:
		expand("results/illumina/{genome}.illumina.bwa.read1.fastq", genome=GENOMES),
		expand("results/illumina/{genome}.illumina.bwa.read2.fastq", genome=GENOMES),
		expand("results/nanopore/{genome}_{read_len}_cov{cov}_reads.fasta", genome=GENOMES, read_len=LENGTH, cov=COVERAGE)

rule simulate_illumina:
	input:
		"fasta/{genome}.fna"
	output:
		"results/illumina/{genome}.illumina.bwa.read1.fastq",
		"results/illumina/{genome}.illumina.bwa.read2.fastq"
	params:
		R1_error="0.0001-0.001",
		R2_error="0.0001-0.001",
		number_reads=4000,
		read_length=250,
		mut_rate=0,
		indel_rate=0,
		random_read=0.01,
		out_prefix="./results/illumina/{genome}.illumina"
	run:
		shell("dwgsim -e {params.R1_error} -E {params.R2_error} -N {params.number_reads} "
		"-1 {params.read_length} -2 {params.read_length} -R {params.indel_rate} -r {params.mut_rate} "
		"-y {params.random_read} {input} {params.out_prefix}")

rule simulate_nanopore:
	input:
		"fasta/{genome}.fna"
	output:
		"results/nanopore/{genome}_{read_len}_reads.fasta"
	params:
		number_reads=200,
		out_prefix="results/nanopore/{genome}_{read_len}",
		error_profile="ecoli",
	run:
		if wildcards.read_len == 'Q1':
			shell("simulator.py circular -r {input} -o {params.out_prefix} -n {params.number_reads} "
			"-c {params.error_profile} --max_len 2500")
		if wildcards.read_len == 'Q2':
			shell("simulator.py circular -r {input} -o {params.out_prefix} -n {params.number_reads} "
			"-c {params.error_profile} --min_len 2500 --max_len 5000")

rule get_total_genome_bp:
	input:
		"fasta/{genome}.fna"
	output:
		"fasta/{genome}_length.txt"
	shell:
		"grep -v '>' {input} | wc -m > {output}"

rule get_total_nanosim_bp:
	input:
		"results/nanopore/{genome}_{read_len}_reads.fasta"
	output:
		"results/nanopore/{genome}_{read_len}_length.txt"
	shell:
		"grep -v '>' {input} | wc -m > {output}"

rule subsample_nanopore:
	input:
		genome_bp="fasta/{genome}_length.txt",
		nanosim_bp="results/nanopore/{genome}_{read_len}_length.txt",
		nano_reads="results/nanopore/{genome}_{read_len}_reads.fasta"
	output:
		"results/nanopore/{genome}_{read_len}_cov{cov}_reads.fasta"
	run:
		#my_num = {wildcards.cov}*{wildcards.cov}
		shell("seqtk sample -s $RANDOM {input.nano_reads} {wildcards.cov} > {output}")

############
# ignore below
############

rule clean:
	run:
		shell("rm results/nanopore/*error_profile")

rule assembly:
	input:
		ont="results/nanopore/{genome}_{read_len}_reads.fasta",
		r1="results/illumina/{genome}.illumina.bwa.read1.fastq",
		r2="results/illumina/{genome}.illumina.bwa.read2.fastq"
	output:
		"results/assemblies/{genome}/{read_len}_contigs.fasta"
	params:
		out_prefix="results/assemblies/{genome}/"
	shell:
		"unicycler -1 {input.r1} -2 {input.r2} -l {input.ont} -o {params.out_prefix}"

rule makeidx:
	input:
		"results/assemblies/{genome}/{read_len}_contigs.fasta"
	output:
		touch("results/assemblies/{genome}/{read_len}_makeidx.done")
	shell:
		"bwa index {input}"

rule map:
	input:
		r1="results/illumina/{genome}.illumina.bwa.read1.fastq",
		r2="results/illumina/{genome}.illumina.bwa.read2.fastq",
		idx="results/assemblies/{genome}/{read_len}_contigs.fasta"
	output:
		"results/mapping/{genome}/{read_len}_illumina.sorted.bam"
	threads:
		8
	shell:
		"bwa mem -t {threads} {input.idx} {input.r1} {input.r2} | samtools sort -o {output}"

rule variants:
	input:
		bam="results/mapping/{genome}/{read_len}_illumina.sorted.bam",
		fasta="results/assemblies/{genome}/{read_len}_contigs.fasta"
	output:
		"results/mapping/{genome}/{read_len}_contigs.vcf"
	shell:
		"samtools mpileup -f {input.fasta} {input.bam} > {output}"

rule depth:
	input:
		"results/mapping/{genome}/{read_len}_illumina.sorted.bam"
	output:
		"results/mapping/{genome}/{read_len}_depth.txt"
	shell:
		"samtools depth {input} > {output}"
