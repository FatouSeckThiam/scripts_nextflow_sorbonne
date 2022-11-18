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
params.samples = "/home/seck-thiam/samplelist_Pvactools_run_10_second_netMHCpan"
params.inputDir = "/SCVOL02/run_10/neoepitope_pvactools/2.0.1/TX/peptide/CMHI"
params.outDir = "/SCVOL02/run_10/neoepitope_pvactools/2.0.1/TX/peptide/CMHI/netMHCpan"



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
        String id_hlaI = line.split("\t")[2].split("_")[0]
	peptide = "${params.inputDir}/"+ id_hlaI+"_uniq_CMHI_MT_Epitopes_Seq.txt"
	hla_allele = line.split("\t")[2].split("_")[1]
        inputs.add([id_hlaI, hla_allele, id_hlaI,[peptide]])
  }

}

hla_channel=Channel.fromList(inputs).view()


process filtre_ideation{
	module "opt-python"

	input:
	tuple val(id), val(hla_allele), val(id_hla), path(peptide) from hla_channel
	
	output:
	tuple val(id), path("${id}_core_netMHC_out.xls") into core_channel

	shell:
	"""
	mkdir -p !{params.outDir}
	
	/SCVOL01/Projects/Ideation/refs/new_iebd/mhc_i/method/netmhcpan-4.1-executable/netmhcpan_4_1_executable/netMHCpan -a !{hla_allele} -p !{params.inputDir}/!{peptide} -s -xls -xlsfile !{id}_core_netMHC_out.xls	
	
	cp !{id}_core_netMHC_out.xls !{params.outDir}
	"""
}


