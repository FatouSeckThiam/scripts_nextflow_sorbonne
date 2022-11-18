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
//params.RNAlist = "/home/seck-thiam/samplelist_rna_run13"
params.RNAlist = "/home/seck-thiam/samples_bam_star_rna_run13"
//params.inputDir = "/SCVOL02/ideation_raw_data/Novaseq_mRNA_220713_karim/fastq"
params.inputDir = "/SCVOL02/run_13/rnaseq/Bam"
params.refStar = "/refs/references/GENOME/Homo_sapiens.hg38/STAR"
params.refGenome = "/refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa"
params.known_sites = "/home/seck-thiam/dbsnp_146.Homo_sapiens.hg38.vcf.gz" 
params.known_sites_indels ="/home/seck-thiam/Mills_and_1000G_gold_standard.indels.Homo_sapiens.hg38.vcf.gz" 
params.refGenome_targets_list="/home/seck-thiam/Twist_ComprehensiveExome_targets_hg38.bed" 
 

params.outDir_Fastp = "/SCVOL02/run_13/rnaseq/Fastp_return"
params.outDir_Bam = "/SCVOL02/run_13/rnaseq/Bam"
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

samplist_rna = file("${params.RNAlist}")
// inputChannel = Channel.empty()

inputs_bam_star = []
reader = samplist_rna.newReader()
samplist_rna.withReader {
    String line
    while( line = reader.readLine() ) {
        String id_bam_star = line.split("_")[0]
        file_bam = "${params.inputDir}/"+id_bam_star+"_sorted_markdup_splitn.bam"
	file_bam_metrics = "${params.inputDir}/"+id_bam_star+"_marked_dup_metrics.txt"
        inputs_bam_star.add([id_bam_star, [file_bam, file_bam_metrics]])

        }	
}

FastqChannel_bam = Channel.fromList(inputs_bam_star).view()

/*process Fastq_rna{

      	module "fastp/0.22.0"
	cpus params.cpus
	time "72h"

        input:
        tuple val(read1_id), path(fastq1), val(read2_id), path(fastq2), val(id) from FastqChannel_rna

	output:
        tuple val(id), path("./${id}_rnaseq_R1.fastq.gz"), path("./${id}_rnaseq_R2.fastq.gz"), path("./${id}.html"), path("./${id}.jason") into fastq_ideation_rna_ch

        shell:
        """
	mkdir -p !{params.outDir_Fastp}
        fastp -i !{params.inputDir}/!{fastq1} -o ./!{id}_rnaseq_R1.fastq.gz -I !{params.inputDir}/!{fastq2} -O ./!{id}_rnaseq_R2.fastq.gz -h ./!{id}.html -j ./!{id}.jason -c
	cp ./!{id}_rnaseq_R1.fastq.gz !{params.outDir_Fastp}
        cp ./!{id}_rnaseq_R2.fastq.gz !{params.outDir_Fastp}
	cp ./!{id}.html !{params.outDir_Fastp}
	cp ./!{id}.jason !{params.outDir_Fastp}

        """
}


process Star_Mapping_rna {
        module "star/2.7.2:samtools"
        cpus params.cpus
	memory "64G"
	time "72h"

        input:
	tuple val(id), path(rnaseq_R1), path(rnaseq_R2), path(html), path(jason) from fastq_ideation_rna_ch
        
	output:
        tuple val(id), path("./${id}_rnaseq.bamAligned.out.bam") into bam_rna_ch

        script:
	"""
	mkdir -p ${params.outDir_Bam}
	mkdir -p ${id}_rnaseq.bam_STARgenome
	mkdir -p ${id}_rnaseq.bam_STARpass1
	mkdir -p ${id}_rnaseq.bam_STARtmp
	
	STAR --genomeDir ${params.refStar} --readFilesIn ${rnaseq_R1} ${rnaseq_R2} --outFileNamePrefix ${id}_rnaseq.bam --readFilesCommand zcat --runThreadN 8 --genomeChrBinNbits 11 --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignSoftClipAtReferenceEnds Yes --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1 --chimOutType Junctions SeparateSAMold WithinBAM SoftClip --chimSegmentMin 15 --outFilterIntronMotifs None --outFilterMatchNminOverLread 0.33 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterType BySJout --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --outSAMunmapped Within --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts
        cp ./${id}_rnaseq.bamAligned.out.bam ${params.outDir_Bam}
	cp ./${id}_rnaseq.bamAligned.toTranscriptome.out.bam ${params.outDir_Bam}
	cp ./${id}_rnaseq.bamLog.out ${params.outDir_Bam}
	cp ./${id}_rnaseq.bamLog.progress.out ${params.outDir_Bam}
	cp ./${id}_rnaseq.bamLog.final.out ${params.outDir_Bam}
	cp ./${id}_rnaseq.bamChimeric.out.junction ${params.outDir_Bam}
	cp ./${id}_rnaseq.bamSJ.out.tab ${params.outDir_Bam}
	cp ./${id}_rnaseq.bamReadsPerGene.out.tab ${params.outDir_Bam}
	cp -R ./${id}_rnaseq.bam_STARgenome ${params.outDir_Bam}
	cp -R ./${id}_rnaseq.bam_STARpass1 ${params.outDir_Bam}
	cp -R ./${id}_rnaseq.bam_STARtmp ${params.outDir_Bam}
	
        """
}



process SortSam_rna{

        module "picard/2.21.7:samtools"
        cpus params.cpus
	memory "10G"
	time "72h"
        
	input: 
	tuple val(id), path(file_bam) from FastqChannel_bam

        output:
        tuple val(id), path("./${id}_rnaseq_sorted.bam") into sort_sam_rna_ch

        shell:
	"""
	java -jar \$PICARD SortSam I=!{file_bam} O=./!{id}_rnaseq_sorted.bam SORT_ORDER=coordinate TMP_DIR=!{params.outDir_Bam} 	
        """
}


process Markupduplicate_rna {
        module "picard/2.18.25:samtools"
        cpus params.cpus
	time "72h"

        input: 
	tuple val(id), path(rnaseq_sorted_bam) from sort_sam_rna_ch        

	output:
        tuple val(id), path("./${id}_rnaseq_markdup.bam"), path("./${id}_rnaseq_markdup.bai"), path("./${id}_marked_dup_metrics.txt") into mardup_rna_ch

        shell:
        """
        java -jar \$PICARD MarkDuplicates I=!{rnaseq_sorted_bam} O=./!{id}_rnaseq_markdup.bam M=./!{id}_marked_dup_metrics.txt CREATE_INDEX=TRUE TMP_DIR=!{params.outDir_Bam}
	        
	cp ./${id}_marked_dup_metrics.txt !{params.outDir_Bam}

	"""
}

*/
process SplitNCigarReads_rna{
        module "gatk/4.1:samtools"
        cpus params.cpus
	time "72h"

        input:
	tuple val(id), path(rnaseq_markdup_bam) from FastqChannel_bam
	
        output:
	tuple val(id), path("./${id}_rnaseq_sorted_markdup_splitn.bam"), path("./${id}_rnaseq_sorted_markdup_splitn.bai") into sng_rna_ch
	
        shell :
        """
	gatk SplitNCigarReads -I !{params.inputDir}/!{rnaseq_markdup_bam[0]} -O ./!{id}_rnaseq_sorted_markdup_splitn.bam -R !{params.refGenome} --tmp-dir !{params.outDir_Bam}
	cp ./!{id}_rnaseq_sorted_markdup_splitn.bam !{params.outDir_Bam}	
	cp ./!{id}_rnaseq_sorted_markdup_splitn.bai !{params.outDir_Bam}
        """
}

process AddOrReplaceReadGroups_rna {
        module "gatk/4.1:picard/2.18.25:samtools"
        cpus params.cpus
	time "72h"

        input:
	tuple val(id), path(rnaseq_sorted_markdup_splitn_bam), path(rnaseq_sorted_markdup_splitn_bai) from sng_rna_ch
			
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
	time "72h"

        input:
	tuple val(id), path(rnaseq_sorted_ar_markdup_splitn_bam), path(rnaseq_sorted_ar_markdup_splitn_bai) from addOrRgroup_rna_ch
	
        output:
        tuple val(id), path("${id}_rnaseq_bqsr_sorted_markdup_splitn.table"), path(rnaseq_sorted_ar_markdup_splitn_bam), path(rnaseq_sorted_ar_markdup_splitn_bai) into basecall_rna_ch 

        shell:
        """
        gatk BaseRecalibrator --input !{rnaseq_sorted_ar_markdup_splitn_bam} --output ./!{id}_rnaseq_bqsr_sorted_markdup_splitn.table --reference !{params.refGenome} --known-sites !{params.known_sites} --known-sites !{params.known_sites_indels} --tmp-dir !{params.outDir_Bam}
        cp ./!{id}_rnaseq_bqsr_sorted_markdup_splitn.table !{params.outDir_Bam}

        """
}

process Apply_BSQR_rna{

        module "gatk/4.1:samtools/1.14"
        cpus params.cpus
	time "72h"

        input:
	tuple val(id), path(rnaseq_bqsr_sorted_markdup_splitn_table), path(rnaseq_sorted_ar_markdup_splitn_bam), path(rnaseq_sorted_ar_markdup_splitn_bai) from basecall_rna_ch
	
        output:
        tuple val(id), path("./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam"),  path("./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bai") into applyBQSR_ch


        script:
        """
        gatk ApplyBQSR --input ${rnaseq_sorted_ar_markdup_splitn_bam} --bqsr-recal-file ${rnaseq_bqsr_sorted_markdup_splitn_table} --output ./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam --reference ${params.refGenome} --tmp-dir /SCVOL01/Projects/Ideation/cache_dir
	cp ./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bam ${params.outDir_Bam}
	cp ./${id}_rnaseq_bqsr_sorted_ar_markdup_splitn.bai ${params.outDir_Bam}
        """  
}




