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
params.Tumlist = "/home/seck-thiam/sampleliste_tum_wes_run13_prep_data"
params.inputDir = "/SCVOL02/Data_transfer_globus_18_02_22/Novaseq_170222_Bouzidi_EXOMES/fastq"

params.refGenome_BWA = "/refs/references/GENOME/Homo_sapiens.hg38/BWA/Homo_sapiens.hg38.fa"
params.refGenome = "/refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa"
params.known_sites = "/home/seck-thiam/dbsnp_146.Homo_sapiens.hg38.vcf.gz" 
params.known_sites_indels ="/home/seck-thiam/Mills_and_1000G_gold_standard.indels.Homo_sapiens.hg38.vcf.gz" //process Basecallibrator_Tumour
params.refGenome_targets_list="/home/seck-thiam/Twist_ComprehensiveExome_targets_hg38.bed.gz" //process DOC
 
//rennomer bam en bam_files pour ne pas prendre en compte l'extension bam et donc éviter cette erreur //SCVOL01/Projects/Ideation/run_06/exome/bai

params.outDir_Fastp = "/SCVOL02/Fatou/Fastp_test_run_11"
params.outDir_Bam = "/SCVOL02/Fatou/Bam_test_run_11"
params.outDir_QCmetrics = "/SCVOL02/Fatou/Bam_test_run_11/QCmetrics"
params.cpus = 5
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
} else if (params.Germlist == " "|| params.Tumlist == " "|| params.inputDir == " " || params.outDir_Fastp == " " || params.outDir_Bam == " " || params.refGenome == " " || params.refGenome_BWA == " "|| params.known_sites == " " || params.known_sites_indels == " " || params.refGenome_targets_list == " " || params.outDir_QCmetrics == " ")

{
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

samplist_Tumour = file("${params.Tumlist}")
// inputChannel = Channel.empty()

inputs_fastq_tumour = []
reader = samplist_Tumour.newReader()
samplist_Tumour.withReader {
    String line
    while( line = reader.readLine() ) {
        String id_fastq1_tum = line.split("\t")[3].replaceAll(".fastq.gz","")
        String id_fastq2_tum = line.split("\t")[4].replaceAll(".fastq.gz","")
        String fq1_tum = line.split("\t")[3]
        String fq2_tum = line.split("\t")[4]
        fastq1_tum = "${params.inputDir}/"+fq1_tum
        fastq2_tum = "${params.inputDir}/"+fq2_tum
        String id_ideation = line.split("\t")[1]
        //tag_ideation_tum = line.split("\t")[2] 
        inputs_fastq_tumour.add([id_fastq1_tum, [fastq1_tum], id_fastq2_tum, [fastq2_tum], id_ideation])
       // inputs_id_ideation.add([id_ideation]])

        }	
}

FastqChannel_Tumour = Channel.fromList(inputs_fastq_tumour).view()

process Fastq_Tumor{

      	module "fastp/0.22.0"
	cpus params.cpus

        input:
        tuple val(read1_tum_id), path(fastq1_tum), val(read2_tum_id), path(fastq2_tum), val(id) from FastqChannel_Tumour

	output:
        tuple val(id), path("./${id}_Tumour_R1.fastq.gz"), path("./${id}_Tumour_R2.fastq.gz"), path("./${id}_Tumour.html"), path("./${id}_Tumour.jason") into fastq_ideation_tumor_ch

        shell:
        """
        fastp -i !{params.inputDir}/!{fastq1_tum} -o ./!{id}_Tumour_R1.fastq.gz -I !{params.inputDir}/!{fastq2_tum} -O ./!{id}_Tumour_R2.fastq.gz -h ./!{id}_Tumour.html -j ./!{id}_Tumour.jason -c
	cp ./!{id}_Tumour_R1.fastq.gz !{params.outDir_Fastp}
        cp ./!{id}_Tumour_R2.fastq.gz !{params.outDir_Fastp}
	cp ./!{id}_Tumour.html !{params.outDir_Fastp}
	cp ./!{id}_Tumour.jason !{params.outDir_Fastp}

        """
}


process Mapping_Tumour {
        module "bwa/0.7.17:samtools"
        cpus params.cpus

        input:
	tuple val(id), path(Tumour_R1), path(Tumour_R2), path(Tumour_html), path(Tumour_jason) from fastq_ideation_tumor_ch

        output:
        tuple val(id), path("./${id}_Tumour.bam") into bam_tumor_ch

        script:
	"""
	mkdir -p ${params.outDir_Bam}
        bwa mem -t 4 -M ${params.refGenome_BWA} ${Tumour_R1} ${Tumour_R2} | samtools view -Sb > ./${id}_Tumour.bam
        cp ./${id}_Tumour.bam ${params.outDir_Bam}
        """
}


process ARgroup_Tumour {

        module "picard/2.18.25:samtools"
        cpus params.cpus

        input: 
        tuple val(id), path(Tumour_bam) from bam_tumor_ch

        output:
        tuple val(id), path("./${id}_Tumour_rg_sorted.bam"), path("./${id}_Tumour_rg_sorted.bai") into ARgroup_tum_ch

        shell:
	"""
        java -jar \$PICARD AddOrReplaceReadGroups I=!{Tumour_bam} O=./!{id}_Tumour_rg_sorted.bam RGPL=illumina RGPU=Novaseq1 RGLB=lib1 RGSM=!{id}_Tumour CREATE_INDEX=TRUE TMP_DIR=/SCVOL01/Projects/Ideation/cache_dir SORT_ORDER=coordinate
	
        """
}


process Markupduplicate_Tumour{
        module "picard/2.18.25:samtools"
        cpus params.cpus

        input: 
	tuple val(id), path(Tumour_rg_sorted_bam), path(Tumour_rg_sorted_bai) from ARgroup_tum_ch
        
	output:
        tuple val(id), path("./${id}_Tumour_markdup_rg_sorted.bam"), path("./${id}_Tumour_markdup_metrics.txt"), path("./${id}_Tumour_markdup_rg_sorted.bai") into Mardup_Tum_ch

        shell:
        """
	mkdir -p !{params.outDir_QCmetrics}

        java -jar \$PICARD MarkDuplicates I=!{Tumour_rg_sorted_bam} O=./!{id}_Tumour_markdup_rg_sorted.bam M=./!{id}_Tumour_markdup_metrics.txt CREATE_INDEX=TRUE TMP_DIR=/SCVOL01/Projects/Ideation/cache_dir
        
	cp ./!{id}_Tumour_markdup_metrics.txt !{params.outDir_QCmetrics}
	"""
}


process Realign_Intervals_Tum{
        module "gatk/3.8:samtools"
        cpus params.cpus

        input:
	tuple val(id), path(Tumour_markdup_rg_sorted_bam), path(Tumour_markdup_metrics_txt), path(Tumour_markdup_rg_sorted_bai) from Mardup_Tum_ch

        output:
        tuple val(id), path ("./${id}_Tumour.realign.intervals"), path(Tumour_markdup_rg_sorted_bam) into Realign_intervals_Tum_ch

        shell :
        """
	mkdir -p !{params.outDir_QCmetrics}
        java -jar \$GATK -T RealignerTargetCreator -I !{Tumour_markdup_rg_sorted_bam} -o ./!{id}_Tumour.realign.intervals -R !{params.refGenome} -L !{params.refGenome_targets_list} -known !{params.known_sites_indels} -allowPotentiallyMisencodedQuals
        cp ./!{id}_Tumour.realign.intervals !{params.outDir_QCmetrics}
        """
}

process Apply_Realing_Tum{
        module "gatk/3.8:samtools"
        cpus params.cpus

        input:
	tuple val(id), path (Tumour_realign_intervals), path(Tumour_markdup_rg_sorted_bam) from Realign_intervals_Tum_ch
	
        output:
        tuple val(id), path ("./${id}_Tumour_realign_markdup_rg_sorted.bam") into Apply_Realign_Tum_ch

        shell:
	"""
        java -jar \$GATK -T IndelRealigner -I !{Tumour_markdup_rg_sorted_bam} -o ./!{id}_Tumour_realign_markdup_rg_sorted.bam -targetIntervals !{Tumour_realign_intervals} -R !{params.refGenome} -known !{params.known_sites_indels} -allowPotentiallyMisencodedQuals
        """
}


process Basecallibrator_Tumour{

        module "gatk/4.1:samtools"
        cpus params.cpus
	

        input:
        tuple val(id), path(Tumour_realign_markdup_rg_sorted_bam) from Apply_Realign_Tum_ch 

        output:
        tuple val(id), path("${id}_Tumour_bqsr_realign_markdup_rg_sorted.table"), path(Tumour_realign_markdup_rg_sorted_bam) into Basecall_Tum_ch 

        shell:
        """
        gatk BaseRecalibrator --input !{Tumour_realign_markdup_rg_sorted_bam} --output ./!{id}_Tumour_bqsr_realign_markdup_rg_sorted.table --reference !{params.refGenome} --known-sites !{params.known_sites} --known-sites !{params.known_sites_indels} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
        cp ./!{id}_Tumour_bqsr_realign_markdup_rg_sorted.table !{params.outDir_QCmetrics}

        """
}

process Apply_Calibrator_Tumour{

        module "gatk/4.1:samtools/1.14"
        cpus params.cpus


        input:
	tuple val(id), path(Tumour_bqsr_realign_markdup_rg_sorted_table), path(Tumour_realign_markdup_rg_sorted_bam) from Basecall_Tum_ch

        output:
        tuple val(id), path ("./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bam"),  path ("./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bai") into Apply_Calibrator_Tum_ch


        script:
        """
        gatk ApplyBQSR --input ${Tumour_realign_markdup_rg_sorted_bam} --bqsr-recal-file ${Tumour_bqsr_realign_markdup_rg_sorted_table} --output ./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bam --reference ${params.refGenome} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
	cp ./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bam ${params.outDir_Bam}
	cp ./${id}_Tumour_bqsr_realign_markdup_rg_sorted.bai ${params.outDir_Bam}
        """  
}

process DOC_Tumour{

        module "gatk/3.8:samtools/1.14"
        cpus params.cpus

        input:
	tuple val(id), path (Tumour_bqsr_realign_markdup_rg_sorted_bam), path(Tumour_bqsr_realign_markdup_rg_sorted_bai) from Apply_Calibrator_Tum_ch

        output:
        tuple val(id), path("./${id}_Tumour_depthofcoverage") into DOC_Tum_ch

        shell:
	"""
	mkdir -p !{params.outDir_QCmetrics}
        java -jar \$GATK -T DepthOfCoverage -I !{Tumour_bqsr_realign_markdup_rg_sorted_bam} -R !{params.refGenome} -L !{params.refGenome_targets_list} -o ./!{id}_Tumour_depthofcoverage
        cp ./!{id}_Tumour_depthofcoverage !{params.outDir_QCmetrics}
	cp ./!{id}_Tumour_depthofcoverage.sample_cumulative_coverage_counts !{params.outDir_QCmetrics}
	cp ./!{id}_Tumour_depthofcoverage.sample_cumulative_coverage_proportions !{params.outDir_QCmetrics}
	cp ./!{id}_Tumour_depthofcoverage.sample_interval_statistics !{params.outDir_QCmetrics}
	cp ./!{id}_Tumour_depthofcoverage.sample_interval_summary !{params.outDir_QCmetrics}
	cp ./!{id}_Tumour_depthofcoverage.sample_statistics !{params.outDir_QCmetrics}
	cp ./!{id}_Tumour_depthofcoverage.sample_summary !{params.outDir_QCmetrics}
        """
}

