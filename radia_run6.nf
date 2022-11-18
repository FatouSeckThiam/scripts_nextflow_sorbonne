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
params.samplelist = "/home/seck-thiam/samplelist_radia_run_11"
params.inputDir = "/SCVOL02/Fatou/Bam_test_run_11"
params.inputDirRna = "/SCVOL02/run_11/Bam_rnaseq_run_11"

params.outDir="/SCVOL02/run_11/vcf_radia/"
params.refGenome="/refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa"

params.cpus = 5
params.help=false


// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline est fait pour réaliser des appels de variants somatiques avec l'outil Mutect 1")
        println(   " Arguments réquis :\n")
        println("  --samplelist: liste des échantillons")
        println("  --inputDir : [PATH] chemin des fichiers bam")
        println("  --output : [PATH] chemin des fichiers de sortie")
        println("  --chr_vcf_gz : [PATH] chemin vcf")
        println("  --run_snp_pileup_facets_R = [PATH] chemin script R facets")
        println(" Arguments optionnels :\n")
        println("  --cpus [INT] = nombre de cpus (4 par défaut)")

}
if (params.help){
	usage()
	exit(1)
 } else if (params.samplelist == " " || params.inputDir == " " || params.outDir == "" || params.refGenome == ""){
	println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")

}

samples=file("${params.samplelist}")
inputs=[]

reader = samples.newReader()
samples.withReader {
	String line
	while(line = reader.readLine() ){
	String normal_bn = line.split(" ")[0].split("_")[0].replaceAll(".bam","")
        String tumor_bn = line.split(" ")[1].split("_")[0].replaceAll(".bam","")
	String rnaseq_tum = line.split(" ")[2].split("_")[0].replaceAll(".bam","")
        String normal = line.split(" ")[0]
        //String normal_dir = line.split(" ")[0].split("_")[0]
        normal_bam = "${params.inputDir}/" + normal //normal
        normal_idx = normal_bam.replaceAll("bam","bam.bai")
        //normal_pileup = "${params.inputDir}/${normal_dir}/${normal_dir}.pileups"
        String tumoral = line.split(" ")[1]
        //String tumor_dir = line.split(" ")[1].split("/")[0]
        tumor_bam = "${params.inputDir}/" + tumoral
        tumor_idx = tumor_bam.replaceAll("bam","bam.bai")
	String rnaseq = line.split(" ")[2]
	rnaseq_bam = "${params.inputDirRna}/"+rnaseq
	rnaseq_idx = rnaseq_bam.replaceAll("bam","bam.bai")
        //tumor_pileup = "${params.inputDir}/${tumor_dir}/${tumor_dir}.pileups"
        //inputs_norm.add([normal_bn, normal_bam, normal_idx])
        inputs.add([normal_bn, [normal_bam, normal_idx], tumor_bn, [tumor_bam, tumor_idx], rnaseq_tum, [rnaseq_bam, rnaseq_idx]])
        //inputs_tumor.add([tumor_bn, [tumor_bam, tumor_idx]])

   }

}

liste_chr=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]

BamChannel=Channel.fromList(inputs).view()
ChromChannel=Channel.fromList(liste_chr).view()
ChromChannelFilter=Channel.fromList(liste_chr).view()

process run_radia{
	
	module "samtools:radia/1.1.4"
	memory "30GB"

	input:
	tuple val(id_norm), path(normal_bam), val(id_tum), path(tumor_bam), val(id_rna), path(rnaseq_bam) from BamChannel 
	each chr from ChromChannel

	output:
	tuple val(id_norm), path("${id_norm}_radia_unfiltered_${chr}.vcf") into radia_channel

	script:
	"""	
	mkdir -p ${params.outDir}
	radia.py ${id_norm} ${chr} -n ${params.inputDir}/${normal_bam[0]} -t ${params.inputDir}/${tumor_bam[0]} -r ${params.inputDirRna}/${rnaseq_bam[0]} -f ${params.refGenome} -o ./${id_norm}_radia_unfiltered_${chr}.vcf
	
	cp ./${id_norm}_radia_unfiltered_${chr}.vcf ${params.outDir}
	"""

}

process radia_filter{

 	module "samtools:radia/1.1.4"
        memory "30GB"

	input:
	tuple val(id_norm), path("${id_norm}_radia_unfiltered_${chr}.vcf") from radia_channel
	each chr from ChromChannelFilter

	output: 
	tuple val(id_norm), path("${id_norm}_radia_filtered_${chr}.vcf")  into radia_filter_channel

	script:
	"""
	python /shared/apps/radia/1.1.4/bin/scripts/filterRadia.py ${id_norm}_radia_filtered ${chr} ${id_norm}_radia_unfiltered_${chr}.vcf ${params.outDir} /shared/apps/radia/1.1.4/bin/scripts/ -b /shared/apps/radia/1.1.4/bin/data/hg38/blacklists/1000Genomes/phase3/ -d /shared/apps/radia/1.1.4/bin/data/hg38/snp151/ -r /shared/apps/radia/1.1.4/bin/data/hg38/retroGenes/ -p /shared/apps/radia/1.1.4/bin/data/hg38/pseudoGenes/ -c /shared/apps/radia/1.1.4/bin/data/hg38/cosmic/ -t /shared/apps/radia/1.1.4/bin/data/hg38/gencode/basic/ --rnaGeneBlckFile /shared/apps/radia/1.1.4/bin/data/rnaGeneBlacklist.tab --rnaGeneFamilyBlckFile /shared/apps/radia/1.1.4/bin/data/rnaGeneFamilyBlacklist.tab --noRadar --noDarned --noSnpEff --noCosmic
	
	cp ./${id_norm}_radia_filtered_${chr}.vcf ${params.outDir}
	"""
}

/*process merge_radia{
	
	module "samtools:radia/1.1.4"
	
	input:
	tuple val(id_norm), path("${id_norm}_radia_filtered_${chr}.vcf") from radia_filter_channel

	output:
	tuple val(id_norm), path("${id_norm}_indels.rna.intersect.variants.raf.vcf"), path("${id_norm}_snv.rna.intersect.variants.raf.vcf") into radia_merge_variant_channel

	shell:
	"""
	python /shared/apps/radia/1.1.4/bin/scripts/mergeChroms.py !{id_norm}_radia_filtered !{params.outDir} !{params.outDir}
}	"""
*/











