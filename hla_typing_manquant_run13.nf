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

params.samples = "/home/seck-thiam/samplelist_seq2hla_run13"
params.inputDir = "/SCVOL02/run_13/exome/fastp"
//params.inputDir= "/SCVOL02/run_11/rnaseq/fastp"
params.outDir = "/SCVOL02/run_13/exome/seq2HLA_Germline"
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


samplelist =file("${params.samples}")

inputs=[]

reader= samplelist.newReader()
samplelist.withReader {
	String line
	while (line = reader.readLine() ){
	String id = line.split("_")[0].replaceAll(".fastq.gz","")
	file_fq1 = line.split("\t")[0]
	file_fq2 = line.split("\t")[1]
	fq1 = "${params.inputDir}/"+file_fq1
	fq2 = "${params.inputDir}/"+file_fq2
	inputs.add([id,[fq1,fq2]])
	}
}

inputChannel = Channel.fromList(inputs).view()



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
	tuple val(id), path(fq) from inputChannel 
	

	shell:
	"""
	mkdir -p !{params.outDir}
	python \$SEQ2HLA/seq2HLA.py -1 !{fq[0]} -2 !{fq[1]} -r !{params.outDir}/!{id} -p !{params.cpus} 
	"""
}


 //output:
        //tuple val(id), path("${id}*") into hla_typing_ch
//cp !{id}* !{params.outDir}
