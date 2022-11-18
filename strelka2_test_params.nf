#!/usr/bin/env nextflow

//params.inputDir = '/SCVOL02/run_10/exome/bam/*_Tumour_bqsr_realign_markdup_rg_sorted.bam'
//samples_ch = Channel.fromPath(params.inputDir).map{
        //tuple(it.name.split('_')[0], it)

//}
//samples_ch.view()




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
params.samplelist = '/home/seck-thiam/samplelist'
params.inputDir = '/SCVOL01/Projects/Ideation/run_06/exome/bam_files'
params.outDir = '/home/seck-thiam'
params.refGenome = '/refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa'
params.intervals = '/refs/references/GENOME/Homo_sapiens.hg38/MedExome_hg38_capture_targets.bed '
params.vcf = '/home/seck-thiam/test_variants_calling'
params.dbsnp = '/refs/references/GENOME/Homo_sapiens.hg38/dbsnp_146.Homo_sapiens.hg38.vcf.gz'
params.callRegions = 'home/seck-thiam/refs/MedExome_hg38_capture_targets.bed.gz'
params.cpus = 4
params.help=false

// 2) Consignes paramètres

def usage(){
	println("\nCe pipeline est fait pour réaliser des appels de variants somatiques avec l'outil Mutect 1")
        println(   " Arguments réquis :\n")
        println("  --samplelist: liste des échantillons")
        println("  --inputDir : [PATH] chemin d'accès pour les fichiers bam")
        println("  --outDir : [PATH] chemin d'accès pour les fichiers de sortie")
        println("  --refGenone : [PATH] génome de référence au format fasta")
        println("  --intervals : [PATH] Chemin d'accès à la liste des intervalles correspondant au kit de capture. Peut être un fichier lit, gtf, vcf ou intervalle")
        println("  --coverage-file : [PATH]...")
        println("  --vcf = ...")
        println("  --dbsnp = ...")
        println("  --indelCandidates = [PATH]")
        println(" Arguments optionnels :\n")
        println("  --cpsu [INT] = nombre de cpus (4 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist == " "|| params.inputDir == " " || params.outDir == " " || params.vcf == " " || params.refGenome == " " || params.intervals == " " || params.dbsnp == " " || params.dbsnp == " " || params.callRegions == " "){
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

samplist = file("${params.samplelist}")
// inputChannel = Channel.empty()
//inputs_norm = []

inputs_tumor = []
inputs = []

reader = samplist.newReader()
samplist.withReader {
    String line
    while( line = reader.readLine() ) {
        String normal_bn = line.split(" ")[0].split("_")[0].replaceAll(".bam","")
        String tumor_bn = line.split(" ")[1].split("_")[0].replaceAll(".bam","")
        String normal = line.split(" ")[0]
        //String normal_dir = line.split(" ")[0].split("_")[0]
        normal_bam = "${params.inputDir}/" + normal //normal
        normal_idx = normal_bam.replaceAll("bam","bai")
        //normal_pileup = "${params.inputDir}/${normal_dir}/${normal_dir}.pileups"
        String tumoral = line.split(" ")[1]
        //String tumor_dir = line.split(" ")[1].split("/")[0]
        tumor_bam = "${params.inputDir}/" + tumoral
        tumor_idx = tumor_bam.replaceAll("bam","bai")
        //tumor_pileup = "${params.inputDir}/${tumor_dir}/${tumor_dir}.pileups"
        //inputs_norm.add([normal_bn, normal_bam, normal_idx])
        inputs.add([normal_bn, [normal_bam, normal_idx], tumor_bn, [tumor_bam, tumor_idx]])
        inputs_tumor.add([tumor_bn, [tumor_bam, tumor_idx]])
    }
}
//normalChannel = Channel.fromList(inputs_norm)//.subscribe{ println it }

TumorChannel = Channel.fromList(inputs_tumor)
//MutectChannel = Channel.fromList(inputs)
MantaChannel = Channel.fromList(inputs)
Strelka2Channel = Channel.fromList(inputs)


process RunManta {
	module "manta/1.6.0"
        errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
        maxRetries 2
        memory { 10.GB * task.attempt}

	input:
	tuple val(norm), path(normal_bam), val(tum), path(tumor_bam) from MantaChannel
	output:
	tuple val(norm), path("${params.vcf}/${norm}_manta/results/variants/candidateSmallIndels.vcf.gz"), path("${params.vcf}/${norm}_manta/results/variants/candidateSmallIndels.vcf.gz.tbi"), path(normal_bam), val(tum), path(tumor_bam) into manta_variants

	shell:
	"""
	mkdir -p !{params.outDir}/VC_Mutect_Strelka2
	mkdir -p !{params.vcf}/!{norm}_manta
	//mkdir !{norm}_manta
	
	\$MANTA/bin/configManta.py --normalBam !{params.inputDir}/!{normal_bam[0]} --tumorBam !{params.inputDir}/!{tumor_bam[0]} --referenceFasta !{params.refGenome} --exome --runDir !{params.vcf}/!{norm}_manta
	python ${params.vcf}/!{norm}_manta/runWorkflow.py
	
	cp -r !{params.vcf}/!{norm}_manta/* . 
	
	"""
}


process RunStrelka2{
	module "strelka/2.9.10"
        errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
        maxRetries 2
        memory { 10.GB * task.attempt}

	input:
	//tuple val(norm), path(normal_bam), val(tum), path(tumor_bam) from Strelka2Channel

	tuple val(norm), path(candidateSmallIndels), path(candidateSmallIndels_index), path(normal_bam), val(tum), path(tumor_bam) from manta_variants
	
	output:
	val(norm) into strelka2_variants

	shell:
	"""
	//mkdir -p !{norm}_strelka2
	mkdir -p !{params.vcf}/!{norm}_strelka2

	\$STRELKA/bin/configureStrelkaSomaticWorkflow.py --normalBam !{params.inputDir}/!{normal_bam[0]} --tumorBam !{params.inputDir}/!{tumor_bam[0]} --referenceFasta !{params.refGenome} --indelCandidates !{candidateSmallIndels} --callRegions !{params.callRegions} --exome --runDir !{params.vcf}/!{norm}_strelka2

	python !{params.vcf}/!{norm}_strelka2/runWorkflow.py -m local -j 4
	
	zcat !{params.vcf}/!{norm}_strelka2/results/variants/somatic.indels.vcf.gz | grep "#" | sed 's/NORMAL/!{norm}_Germline/g' | sed 's/TUMOR/!{norm}_Tumour/g' > !{params.vcf}/!{norm}_strelka.indels.pass.vcfcandidateSmallIndels.vcf.gz'
	zcat !{params.vcf}/!{norm}_strelka2/results/variants/somatic.indels.vcf.gz | grep -v "#" | grep PASS >> !{params.vcf}/!{norm}_strelka.indels.pass.vcf
	cp -r !{params.vcf}/!{norm}_strelka2/* .
	"""
}


//process MergeMutectStrelka{
//	module "picard/2.21.7"
//	errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
  //      maxRetries 2
    //    memory { 10.GB * task.attempt}

//	input:
	//tuple val(norm), file ("*_strelka.indels.pass.vcf") from strelka2_variants

//	output:
//	tuple val(norm), file("*_snv.indels.wes.variants.vcf") into all_variants_ch

//	script:
//	"""
//	java -jar \$PICARD MergeVcfs INPUT=${params.vcf}/${norm}_mutect1.snvs.pass.vcf INPUT=${params.vcf}/${norm}_strelka.indels.pass.vcf OUTPUT=${params.vcf}/${norm}_snv.indels.wes.variants.vcf
//	"""
//}
