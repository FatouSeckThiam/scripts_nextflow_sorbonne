#!/usr/bin/env nextflow
workflow.onComplete {
    //any worlflow property can be used here
    if ( workflow.success ) {
        println "Pipeline Complete"
    }
    println "Command line: $workflow.commandLine"
}

workflow.onError {
    println "Oops .. something went wrong"
}

params.inputDir = "/SCVOL02/run_13/rnaseq/Fastp/"
params.outDir = "/SCVOL02/run_13/rnaseq/expression/"
params.cpus = 4
params.index = "/refs/references/GENOME/Homo_sapiens.hg38/KALLISTO/Homo_sapiens.hg38.idx"

//genome_index = file(params.index)

inputChannel = Channel.fromFilePairs("${params.inputDir}/*R{1,2}.fastq.gz").ifEmpty { exit 1, "Cannot find any PE reads file in ${params.inputDir}" }.view() //start from fastq
log.info "-------  Kallisto P I P E L I N E  --------------"

log.info ""
log.info "Current home       : $HOME"
log.info "Current user       : $USER"
log.info "Current path       : $PWD"
log.info "Script dir         : $baseDir"
log.info "Working dir        : $workDir"
log.info "Output dir         : ${params.outDir}"
log.info ""


process kallisto{

	module "kallisto/0.45.0"
	memory "40G"

	input:
	tuple val(id), file(reads) from inputChannel


	output:
	tuple val(id), path("${id}/abundance.h5"), path("${id}/abundance.tsv"), path("${id}/run_info.json") into results

	shell:
	"""
	mkdir -p !{id}
	mkdir -p !{params.outDir}/!{id}
	
	kallisto quant -b 10 -t 4 -i !{params.index} -o ./!{id} !{reads[0]} !{reads[1]}
	cp -r ./!{id}/* !{params.outDir}/!{id}
	"""
}
