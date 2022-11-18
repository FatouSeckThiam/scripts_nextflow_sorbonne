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
//params.samples = "/home/seck-thiam/samplelist_Pvactools_run_11_germline"
params.samples = "/home/seck-thiam/samplelist_filtre_manquant_r11_germ"

params.inputDir = "/SCVOL02/run_11/neoepitope_pvactools/2.0.1/Germline"
params.outDir = "/SCVOL02/run_11/neoepitope_pvactools/2.0.1/Germline/filtre_ideation"
params.outDir_pep ="/SCVOL02/run_11/neoepitope_pvactools/2.0.1/Germline/peptide"

//params.outDir_pep ="/SCVOL02/run_10/neoepitope_pvactools/2.0.1/TX/peptide"



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
        String id_tsv = line.split("\t")[3].split("_")[0].replaceAll(".vcf","")
        //String tsv = line.split("\t")[0]	
	tsv_file = "${params.inputDir}/"+id_tsv+"_wes_refait"+"/combined/"+ id_tsv+"_Tumour.all_epitopes.tsv"
        inputs.add([id_tsv,[tsv_file]])
  }

}

tsv_channel=Channel.fromList(inputs).view()


process filtre_ideation{

	input:
	tuple val(id), path(tsv) from tsv_channel
	
	output:
	tuple val(id), path("${id}_filtre_ideation_tx.tsv") into filtre_channel

	shell:
	"""
	mkdir -p !{params.outDir}
	head -1 !{tsv_file} > ./!{id}_filtre_ideation_tx.tsv 
	awk '{if (\$28 >= 5 && \$29 >= 0.1 && \$36 <= 500 && \$7 == 1) print }' !{tsv} >> ./!{id}_filtre_ideation_tx.tsv 
	cp ./!{id}_filtre_ideation_tx.tsv !{params.outDir} 
	"""
}

process Peptide_neo{

	input:
	tuple val(id), path(filtre_ideation_tx_tsv) from filtre_channel
	
	output:
	tuple val(id), path("${id}_uniq_CMHI_MT_Epitopes_Seq.txt") into peptide_channel

	shell:
	"""
	mkdir -p !{params.outDir_pep}/CMHI
	mkdir -p !{params.outDir_pep}/CMHII

	grep -v "NetMHCIIpan" !{filtre_ideation_tx_tsv} | cut -f 19 | sort | uniq > ./!{id}_uniq_CMHI_MT_Epitopes_Seq.txt
	awk '{if (\$21 == "NetMHCIIpan") print }' !{filtre_ideation_tx_tsv} | cut -f 19 | sort | uniq > ./!{id}_uniq_CMHII_MT_Epitopes_Seq.txt

	cp ./!{id}_uniq_CMHI_MT_Epitopes_Seq.txt !{params.outDir_pep}/CMHI
	cp ./!{id}_uniq_CMHII_MT_Epitopes_Seq.txt !{params.outDir_pep}/CMHII
	
	"""
}
