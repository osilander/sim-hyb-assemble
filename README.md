# sim-hyb-assemble
Simulates illumina and ONT reads to test effects on assembly of bacterial genomes.
Consists only of a Snakefile at this point, requires a .fna file containing a genome.
Uses dwgsim and nanosim to simulate paired end Illumina and ONT reads, respectively.
The ONT reads are simulated at different reads lengths, and downsampled to different coverage values.
The two read sets are then assembled using a hybrid approach implemented by Unicycler.
The assembly is compared to the original fasta using quast.
Depends:
	- dwgsim
	- nanosim
	- quast
	- unicycler
	- seqtk