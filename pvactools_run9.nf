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
params.samples = "/home/seck-thiam/samplelist_Pvactools_run_09"

//params.samples = "/home/seck-thiam/sample_run9_FR-03-115-Pvactools-GX"
params.inputDir = "/SCVOL02/run_09/neoepitope_vcf/GX_with_value"
params.outDir_v3 = "/SCVOL02/run_09/neoepitope_pvactools/3.0.3/GX_values"
params.outDir_v2 = "/SCVOL02/run_09/neoepitope_pvactools/2.0.1/GX_values"
params.outDir_v1 = "/SCVOL02/run_09/neoepitope_pvactools/1.5.4/GX_values"
params.VEP_plugins="/home/seck-thiam/VEP_plugins"

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
} else if (params.samples == " "|| params.inputDir == " " || params.outDir_v2 == " " || params.outDir_v1 == " " || params.input_allele_hla == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

samplelist=file("${params.samples}")
inputs=[]

reader=samplelist.newReader()
samplelist.withReader() {
        String line
        while( line=reader.readLine() ){
        String id_vcf = line.split("\t")[1].split("_")[0].replaceAll(".vcf","")
        String vcf = line.split("\t")[1]
	String id_hla_I = line.split("\t")[3].split("_")[0].replaceAll("_classI","")
	hla_al_I = line.split("\t")[3].split("_")[1]
	//hla_al_I = "${params.input_allele_hla}/"+hla_I
	String id_hla_II = line.split("\t")[2].split("_")[0].replaceAll("_classII","")
	hla_al_II= line.split("\t")[2].split("_")[1]
	//hla_al_II = "${params.input_allele_hla}/"+hla_II
        vcf_file ="${params.inputDir}/"+vcf
        inputs.add([id_vcf,[vcf_file],id_hla_I,hla_al_I, id_hla_II,hla_al_II])
  }

}

vcf_channel1=Channel.fromList(inputs).view()
vcf_channel2=Channel.fromList(inputs)


process run_pvactools2{

        module "conda3/4.7.12:netMHC/4.0-4.0-3.2:netMHCstabpan/0.1:opt-python"
        cpus params.cpus
        time "396h"

        input:
	tuple val(id), path(vcf_file), val(id_hla_I), val(hla_al_I), val(id_hla_II), val(hla_al_II) from vcf_channel1

        output:
        tuple val(id), path("${id}/combined/${id}_Tumour.all_epitopes.tsv") into combine_neo_channel1

        shell:
	"""
	mkdir -p !{id}
        mkdir -p !{id}/combined/
        mkdir -p !{id}/MHC_Class_I/
        mkdir -p !{id}/MHC_Class_II/
        mkdir -p !{params.outDir_v2}/!{id}

        conda run -p /SCVOL01/Tools/pVACtools-2.0.1 pvacseq run !{params.inputDir}/!{vcf_file} !{id}_Tumour --normal-sample-name !{id}_Germline !{hla_al_I}!{hla_al_II} NetMHC NetMHCpan NetMHCIIpan ./!{id} -e1 8,9,10 -e2 15 --iedb-install-directory /SCVOL01/Projects/Ideation/refs/new_iebd/ --tdna-vaf 0.05 --expn-val 1
        cp -r ./!{id}/* !{params.outDir_v2}/!{id}
        """
}


process run_pvactools3{

	module "conda3/4.7.12:netMHC/4.0-4.0-3.2:netMHCstabpan/0.1:opt-python"
	time "384h"	

	input:
	tuple val(id), path(vcf_file), val(id_hla_I), val(hla_al_I), val(id_hla_II), val(hla_al_II) from vcf_channel2
	
	output:
	tuple val(id), path("${id}/combined/${id}_Tumour.all_epitopes.tsv") into combine_neo_channel 

	shell:
	"""
	mkdir -p !{id}
	mkdir -p !{id}/combined/
	mkdir -p !{id}/MHC_Class_I/
	mkdir -p !{id}/MHC_Class_II/
	mkdir -p !{params.outDir_v3}/!{id}
		
	conda run -p /SCVOL01/Tools/pVACtools-3.0.3 pvacseq run !{params.inputDir}/!{vcf_file} !{id}_Tumour --normal-sample-name !{id}_Germline !{hla_al_I}!{hla_al_II} NetMHC NetMHCpan NetMHCIIpan ./!{id} -e1 8,9,10 -e2 15 --iedb-install-directory /SCVOL01/Projects/Ideation/refs/new_iebd/ --tdna-vaf 0.05 --expn-val 1
	cp -r ./!{id}/* !{params.outDir_v3}/!{id}
	"""
}


process run_pvactools1{

	module "conda3/4.7.12:netMHC/4.0-4.0-3.2:netMHCstabpan/0.1:opt-python"
	time "384h"

	input:
	tuple val(id), path("${id}/combined/${id}_Tumour.all_epitopes.tsv") from combine_neo_channel

	output:
	tuple val(id), path("${id}_Tumour.filtered.condensed.ranked.tsv") into v1_channel

	shell:
	"""
	mkdir -p !{params.outDir_v1}
	conda run -p /SCVOL01/Tools/pVACtools-1.5.4 pvacseq generate_condensed_ranked_report -m median !{id}/combined/!{id}_Tumour.all_epitopes.tsv ./!{id}_Tumour.filtered.condensed.ranked.tsv
	cp ./!{id}_Tumour.filtered.condensed.ranked.tsv !{params.outDir_v1}
	"""

}

