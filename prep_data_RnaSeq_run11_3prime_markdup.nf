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

params.Bamlist_markdup = "/home/seck-thiam/markdup_list_11"
params.inputDir_bam = "/SCVOL02/run_11/Bam_rnaseq_run_11"
params.refStar = "/refs/references/GENOME/Homo_sapiens.hg38/STAR"
params.refGenome = "/refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa"
params.known_sites = "/home/seck-thiam/dbsnp_146.Homo_sapiens.hg38.vcf.gz" 
params.known_sites_indels ="/home/seck-thiam/Mills_and_1000G_gold_standard.indels.Homo_sapiens.hg38.vcf.gz" //process Basecallibrator_Tumour
params.refGenome_targets_list="/refs/references/GENOME/Homo_sapiens.hg38/MedExome_Homo_sapiens.hg38_capture_targets.interval_list" //process DOC
 
//rennomer bam en bam_files pour ne pas prendre en compte l'extension bam et donc éviter cette erreur //SCVOL01/Projects/Ideation/run_06/exome/bai

//params.outDir_Fastp = "/SCVOL02/run_12/Fastp_rnaseq3_run_11"
params.outDir_Bam = "/SCVOL02/run_12/Bam_rnaseq3_run_11"
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
        println("  --cpsu [INT] = nombre de cpus (4 par défaut)")

}
if (params.help) {
        usage()
        exit(1)
} else if ( params.RNAlist == " "|| params.inputDir == " " || params.outDir_Fastp == " " || params.outDir_Bam == " " || params.refGenome == " " || params.known_sites == " " || params.known_sites_indels == " " || params.refGenome_targets_list == " " || params.refStar == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}


samplist_bam = file("${params.Bamlist_markdup}")

input_bam_markdup = []
reader= samplist_bam.newReader()
samplist_bam.withReader {
        String line
        while (line = reader.readLine()	){ 
        String id_Bam=line.split("_")[0].replaceAll(".bam","")
        Bam = "${params.inputDir_bam}/"+line
	//Bai = Bam.replaceAll("bam","bai")
        input_bam_markdup.add([id_Bam,[Bam]])
  }
}

BamChannel_markdup = Channel.fromList(input_bam_markdup).view()



process Markupduplicate_rna {
        module "picard/2.18.25:samtools"
        cpus params.cpus

        input:
	tuple val(id), path(Bam) from BamChannel_markdup

        output:
        tuple val(id), path("./${id}_rnaseq_markdup.bam"), path("./${id}_rnaseq_markdup.bai"), path("./${id}_marked_dup_metrics.txt") into mardup_rna_ch

        shell:
	"""
	java -jar \$PICARD MarkDuplicates I=!{Bam} O=./!{id}_rnaseq_markdup.bam M=./!{id}_marked_dup_metrics.txt CREATE_INDEX=TRUE TMP_DIR=!{params.outDir_Bam}

        cp ./!{id}_marked_dup_metrics.txt !{params.outDir_Bam}
        """
}


process SplitNCigarReads_rna{

	module "gatk:samtools"
        cpus params.cpus

        input:
	tuple val(id), path(rnaseq_markdup_bam), path(rnaseq_markdup_bai), path(marked_dup_metrics_txt) from mardup_rna_ch
	
        output:
	tuple val(id), path("./${id}_rnaseq_sorted_markdup_splitn.bam") into sng_rna_ch
	
        script:
        """
	gatk SplitNCigarReads -I ${rnaseq_markdup_bam} -O ./${id}_rnaseq_sorted_markdup_splitn.bam -R ${params.refGenome} --tmp-dir ${params.outDir_Bam}
        """
}

process AddOrReplaceReadGroups_rna{
        module "gatk/4.1:picard/2.18.25:samtools"
        cpus params.cpus

        input:
	tuple val(id), path(rnaseq_sorted_markdup_splitn_bam) from sng_rna_ch
			
        output:
        tuple val(id), path ("./${id}_rnaseq_sorted_ar_markdup_splitn.bam"), path("./${id}_rnaseq_sorted_ar_markdup_splitn.bai") into addOrRgroup_rna_ch

        shell:
	"""
	java -jar \$PICARD AddOrReplaceReadGroups I=!{rnaseq_sorted_markdup_splitn_bam} O=./!{id}_rnaseq_sorted_ar_markdup_splitn.bam RGPL=illumina RGPU=Novaseq1 RGLB=rnaseq RGSM=!{id} CREATE_INDEX=TRUE TMP_DIR=/SCVOL01/Projects/Ideation/cache_dir
        """
}


process Basecallibrator_rna{

        module "gatk/4.1:samtools"
        cpus params.cpus
	

        input:
	tuple val(id), path(rnaseq_sorted_ar_markdup_splitn_bam), path(rnaseq_sorted_ar_markdup_splitn_bai) from addOrRgroup_rna_ch
	
        output:
        tuple val(id), path("./${id}_rnaseq_bqsr_sorted_markdup_splitn.table"), path(rnaseq_sorted_ar_markdup_splitn_bam) into basecall_rna_ch 

        shell:
        """
        gatk BaseRecalibrator --input !{rnaseq_sorted_ar_markdup_splitn_bam} --output ./!{id}_rnaseq_bqsr_sorted_markdup_splitn.table --reference !{params.refGenome} --known-sites !{params.known_sites} --known-sites !{params.known_sites_indels} --tmp-dir !{params.outDir_Bam}
        cp ./!{id}_rnaseq_bqsr_sorted_markdup_splitn.table !{params.outDir_Bam}

        """
}

process Apply_BSQR_rna{

        module "gatk/4.1:samtools/1.14"
        cpus params.cpus


        input:
	tuple val(id), path(rnaseq_bqsr_sorted_markdup_splitn_table), path(rnaseq_sorted_ar_markdup_splitn_bam) from basecall_rna_ch
	
        output:
        tuple val(id), path("./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam"), path("./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bai") into applyBQSR_ch


        script:
        """
        gatk ApplyBQSR --input ${rnaseq_sorted_ar_markdup_splitn_bam} --bqsr-recal-file ${rnaseq_bqsr_sorted_markdup_splitn_table} --output ./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam --reference ${params.refGenome} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
	cp ./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam ${params.outDir_Bam}
	cp ./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bai ${params.outDir_Bam}
        """  
}

