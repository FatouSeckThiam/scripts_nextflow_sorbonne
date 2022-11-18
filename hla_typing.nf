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

params.inputDir = "/SCVOL02/run_ideation_Hors/exome/Fastp"
params.outDir = "/SCVOL02/run_ideation_Hors/exome/seq2HLA_Germline"
params.cpus = 8 
params.help= false

// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline est fait pour réaliser des appels de variants somatiques avec l'outil Mutect 1")
        println(   " Arguments réquis :\n")
        println("  --samplelist: liste des échantillons")
        println("  --inputDir : [PATH] chemin d'accès pour les fichiers bam")
        println("  --outDir : [PATH] chemin d'accès pour les fichiers de sortie")
        println(" Arguments optionnels :\n")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if (params.inputDir == " " || params.outDir == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}


inputChannel = Channel.fromFilePairs("${params.inputDir}/*Germline_R{1,2}.fastq.gz").ifEmpty { exit 1, "Cannot find any PE reads file in ${params.inputDir}" }.view() //start from fastq



/*log.info "-------  HLA_TYPING P I P E L I N E  --------------"

log.info ""
log.info "Current home       : $HOME"
log.info "Current user       : $USER"
log.info "Current path       : $PWD"
log.info "Script dir         : $baseDir"
log.info "Working dir        : $workDir"
log.info "Output dir         : ${params.outDir}"
log.info ""

*/

process seq2hla{

	module "opt-python:samtools:R:seq2HLA/2.2"
	memory "80G"
	cpus params.cpus

	input:
	tuple val(id), file(reads) from inputChannel 
	

	shell:
	"""
	mkdir -p !{params.outDir}
	python \$SEQ2HLA/seq2HLA.py -1 !{reads[0]} -2 !{reads[1]} -r !{params.outDir}/!{id} -p !{params.cpus} 
	"""
}


