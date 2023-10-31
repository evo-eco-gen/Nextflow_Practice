// Commented nextflow script for evaluation of RNAseq input

/*
 * pipeline input parameters
 */
 
// DEFINE PARAMETERS
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
 
// DEFINE TEXT OUTPUT TO LOG
 
log.info """\
 	RNASEQ PIPELINE
 	+++++++++++++++
 	transcriptome	: ${params.transcriptome_file}
 	reads			: ${params.reads}
 	outdir			: ${params.outdir}
	"""
 	.stripIndent()

// DEFINE INDIVIDUAL PROCESSES, with input, output and script


// define the INDEX process that creates the index of a binary file (salmon Step 1)'

process INDEX {
	input:
	path transcriptome
	
	output:
	path 'salmon_index'
	
	script:
	"""
	salmon index --threads $task.cpus -t $transcriptome -i salmon_index
	"""
}

// define the QUANTIFICATION process to assess expression levels (salmon Step 2)

process QUANTIFICATION {
	// Add a tag for diplay
	// ?? Where does sample_id come from?
		// A: it is globbed from 'data/ggal/*_{1,2}.fq'
	// REVIEW the use of variables in bash.
	tag "Salmon on $sample_id"
	publishDir params.outdir, mode:'copy'
	
	input:
	path salmon_index
	// Use a tuple to connect files with paths
	tuple val(sample_id), path(reads)
	
	output:
	path "$sample_id"
	
	script:
	"""
	salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
	"""
}

// process FASTQC to evaluate the quality of individual reads

process FASTQC {
	tag "FASTQC on $sample_id"
	
	input:
		// Notice sample_id without "$"
	tuple val(sample_id), path(reads)
	
	output:
	path "fastqc_${sample_id}_logs"
	
	script:
	"""
	mkdir fastqc_${sample_id}_logs
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""
}

process MULTIQC {
		/*
		 * Use 'publishDir' to point where to store outputs.
		 * Downstream should access outputs through channels.
		 */
	publishDir params.outdir, mode:'copy'
	
		// ?? What's the default path here?
	input:
	path '*'
	
	output:
	path 'multiqc_report.html'
	
	script:
	"""
	multiqc .
	"""
}

/* 
 * Set up the WORKFLOW 
 * using CHANNELS that run PROCESSES over PARAMETERS:
 */

workflow {
	Channel
		.fromFilePairs(params.reads, checkIfExists: true)
		.set { read_pairs_ch }
		
	index_ch = INDEX(params.transcriptome_file)
	quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
	fastqc_ch = FASTQC(read_pairs_ch)
		// MIX two channels and show outputs together with COLLECT
	MULTIQC(quant_ch.mix(fastqc_ch).collect())
}
