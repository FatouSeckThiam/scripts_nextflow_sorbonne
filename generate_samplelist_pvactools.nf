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
params.samples = "/home/seck-thiam/samplelist_vcf_neo_run_10"
//params.samples = "/home/seck-thiam/samplelist_vcf_neo_run_13"
//params.samples = "/home/seck-thiam/samplelist_extract_alleles_germ_run_13"
params.inputDir = "/SCVOL02/run_10/neoepitope/seq2HLA"
params.outDir = "/home/seck-thiam/HLA_alleles_run_10"

params.VEP_plugins="/home/seck-thiam/VEP_plugins"
params.cpus = 4
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
if (params.help){
	usage()
	exit(1)
} else if (params.samples == " "|| params.inputDir == " " || params.outDir == " ")

{
	println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}



samplelist=file("${params.samples}")
liste=[]

reader=samplelist.newReader()
samplelist.withReader() {
	String line
	while(line=reader.readLine()) {
	String id_sample = line.split("\t")[0]
	classe1="${params.inputDir}/"+id_sample+"-ClassI.HLAgenotype4digits"
	classe2="${params.inputDir}/"+id_sample+"-ClassII.HLAgenotype4digits"
	liste.add([id_sample,[classe1,classe2]])
	}
}

hla2_channel = Channel.fromList(liste).view()
hla1_channel = Channel.fromList(liste).view()

process allele1{

        input:
        tuple val(id), path(classe) from hla1_channel

        output:
        tuple val(id), path("${id}_hla1") into allele1_extract

        shell:
        """
        mkdir -p !{params.outDir}

        grep -v "#" !{params.inputDir}/!{classe[0]} | awk '{print \$2, \$4}' | tr -d "'" | awk '{print "HLA-"\$1, "HLA-"\$2}'| awk 'BEGIN { ORS = " " } { print }' > ./!{id}_hla1
        sed -i 's/ /,/g' !{id}_hla1
        sed -i 's/^/!{id}_/' !{id}_hla1
	cp ./!{id}_hla1 !{params.outDir}
        """
}

process allele2{

	input:
	tuple val(id), path(classe) from hla2_channel
	
	output:
	tuple val(id), path("${id}_hla2") into allele2_extract

	shell:
	"""
	mkdir -p !{params.outDir}
	grep -v "#" !{params.inputDir}/!{classe[1]} | awk '{print \$2, \$4}' | tr -d "'" | awk 'BEGIN { ORS = " " } { print }' > ./!{id}_hla2
	sed -i 's/ /,/g' !{id}_hla2 
	sed -i 's/,\$//' !{id}_hla2
	sed -i 's/^/!{id}_/' !{id}_hla2
	cp ./!{id}_hla2 !{params.outDir}
	"""
}

/*
group1=allele1_extract.groupTuple(by: 0).view()
process concat{
	
	input:
	tuple val(id), path(hla1) from group1

	output:
	file "fichier1.txt"
	
	script:
	"""
	for hla in ${params.outDir}
	do
		cat \$hla > file1.txt
		 
	done
	
	cp fichier1.txt ${params.outDir} 
	"""
	
}

*/
