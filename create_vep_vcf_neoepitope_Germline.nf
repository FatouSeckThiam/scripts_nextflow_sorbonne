#!/usr/bin/env nextflow

//params.inputDir = '/SCVOL02/run_10/exome/bam/*_Tumour_bqsr_realign_markdup_rg_sorted.bam'
//samples_ch = Channel.fromPath(params.inputDir).map{
        //tuple(it.name.split('_')[0], it)

//}

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

// 1) Définitions des paramètres
//params.samples = "/home/seck-thiam/samplelist_vcf_neo_run_11"
//params.samples = "/home/seck-thiam/samplelist_vcf_germ_neo_run_10"
params.samples = "/home/seck-thiam/samplelist_create_vcf_germ_run_hors_idea"
params.inputDir_vcf_converted = "/SCVOL02/run_ideation_Hors/exome/neoepitope_vcf"
params.outDir = "/SCVOL02/run_ideation_Hors/neoepitope_vcf/Germline"
params.refGenome = "/home/seck-thiam/Homo_sapiens.hg38.fa"
params.VEP_plugins="/SCVOL02/VEP_plugins"
params.cpus = 8
params.help=false
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
} else if (params.samples == " "|| params.inputDir_vcf_converted == " " || params.outDir == " " || params.vcf_radia == " " || params.input_Bam_rnaseq == " " || params.input_kallisto == " " || params.refGenome == " ") {
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

samplelist=file("${params.samples}")

inputs=[]

reader=samplelist.newReader()
samplelist.withReader() {
	String line
	while( line=reader.readLine() ){
	
	String id_mutect1 = line.split("\t")[1].split("_")[0].replaceAll(".vcf","")
	String id_strelka2 = line.split("\t")[2].split("_")[0].replaceAll(".vcf","")
	String mutect1 = line.split("\t")[1]
	String strelka2 = line.split("\t")[2]
	vcf_mutect1 ="${params.inputDir_vcf_converted}/"+mutect1
	vcf_strelka2 = "${params.inputDir_vcf_converted}/"+strelka2
	inputs.add([id_mutect1,[vcf_mutect1], id_strelka2,[vcf_strelka2]])
  }
 
}

vcf_converted_mutect1_strelka2_channel=Channel.fromList(inputs).view()
vcf_converted_mutect1=Channel.fromList(inputs)


process merge_vcf_converted {

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	memory "40GB"

	input:
	tuple val(id), path(vcf_mutect1), val(id_s), path(vcf_strelka2) from vcf_converted_mutect1_strelka2_channel	
	
	output:
	tuple val(id), path("${id}_mutect1_strelka.pass.converted.vcf") into strelk2_mutect1_variants
	
	shell:
	"""
	mkdir -p !{params.outDir}

	java -jar \$PICARD MergeVcfs INPUT=!{params.inputDir_vcf_converted}/!{vcf_mutect1} INPUT=!{params.inputDir_vcf_converted}/!{vcf_strelka2} OUTPUT=./!{id}_mutect1_strelka.pass.converted.vcf
	
	cp ./!{id}_mutect1_strelka.pass.converted.vcf !{params.outDir}	

	"""
}

process vep {
	module "vep/release-99"
	cpus params.cpus
	
	input:
	tuple val(id), path(mutect1_strelka_pass_converted_vcf) from strelk2_mutect1_variants
	
	output:
	tuple val(id), path("${id}_mutect1_strelka.pass.converted.vep.vcf") into vep_vcf_channel

	shell:
	"""
	vep --input_file !{mutect1_strelka_pass_converted_vcf} --output_file ./!{id}_mutect1_strelka.pass.converted.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta !{params.refGenome} --offline --cache --dir_cache /SCVOL01/Tools/ensembl-vep-99.0/.vep ensembl --assembly GRCh38 --plugin Frameshift --plugin Wildtype --dir_plugins !{params.VEP_plugins}

	cp ./!{id}_mutect1_strelka.pass.converted.vep.vcf !{params.outDir}
	"""
}






