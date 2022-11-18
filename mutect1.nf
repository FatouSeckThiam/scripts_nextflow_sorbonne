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
params.samplelist = "/home/seck-thiam/samplelist"
params.inputDir = "/SCVOL01/Projects/Ideation/run_06/exome/Bam_files"
//rennomer bam en bam_files pour ne pas prendre en compte l'extension bam et donc éviter cette erreur //SCVOL01/Projects/Ideation/run_06/exome/bai

params.outDir = "/home/seck-thiam"
params.refGenome = "/refs/references/GENOME/Homo_sapiens.hg38/Homo_sapiens.hg38.fa"
params.intervals = "/refs/references/GENOME/Homo_sapiens.hg38/MedExome_hg38_capture_targets.bed"
params.vcf = "/SCVOL02/Fatou/TEST"
params.dbsnp = "/refs/references/GENOME/Homo_sapiens.hg38/dbsnp_146.Homo_sapiens.hg38.vcf.gz"
//params.callRegions = "/home/seck-thiam/MedExome_hg38_capture_targets.bed.gz"
params.cpus = 4
params.help=false

//SCVOL01/Projects/Ideation/run_06/exome/bai

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
MutectChannel = Channel.fromList(inputs)
MantaChannel = Channel.fromList(inputs)
Strelka2Channel = Channel.fromList(inputs)

process RunMutect1{

        module "gatk/4.1:mutect/1.1.7"
        //errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
        //maxRetries 2
        //memory { 10.GB * task.attempt}
	cpus params.cpus
        //params.inputBamM = '/SCVOL01/Projects/Ideation/run_06/exome/bam/*_sorted.bam'
        //mutect1_samples_ch = Channel.fromPath(params.inputBamM).map{it.name.split('_')[0]}

        input:
        tuple val(norm), path(normal_bam), val(tum), path(tumor_bam) from MutectChannel

        output:
        tuple val(norm), path("./${norm}_call_stats.out"), path("./${norm}_filtered.vcf"), path(tumor_bam) into mutect1_variants_ch

        shell:
        """
	mkdir -p !{params.outDir}/VC_Mutect1_Strelka2
       	mkdir -p !{params.vcf}
        
        \$JAVA7 -Xmx20g -jar \$MUTECT1 --analysis_type MuTect --reference_sequence !{params.refGenome} --dbsnp !{params.dbsnp} --intervals !{params.intervals} --input_file:normal !{params.inputDir}/!{normal_bam[0]} --input_file:tumor !{params.inputDir}/!{tumor_bam[0]} --out ./!{norm}_call_stats.out --coverage_file ./!{norm}_coverage.wig.txt --vcf ./!{norm}_filtered.vcf
	cp ./!{norm}_call_stats.out !{params.vcf}/
        cp ./!{norm}_filtered.vcf !{params.vcf}/
        cp ./!{norm}_coverage.wig.txt !{params.vcf}/
        """
 }

//filterChannel = mutect1_variants_ch.concat(TumorChannel).groupTuple(by: 0,size: 2).map {it.flatten().toList()}

process Mutect1Filter{

        module "bedtools/2.26.0:opt-python"
      	//errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
        //maxRetries 2
        //memory { 10.GB * task.attempt}
	cpus params.cpus
 
        input:
 	tuple val(norm), path(calll_stats_out), path(filtered_vcf), path(tumor_bam) from mutect1_variants_ch
	//tuple val(norm), path(tumor_bam) from TumorChannel

        output:
	tuple val(norm), file("./${norm}_call_stats.pass.out"), file("./${norm}_filtered_PASS.vcf"), path("./${norm}_mutect1.snvs.pass.vcf") into mutect1_filter_ch

        script:
        """
        grep "#" ./${norm}_filtered.vcf > ./${norm}_filtered_PASS.vcf ; grep -v "#" ./${norm}_filtered.vcf | grep PASS >> ./${norm}_filtered_PASS.vcf
 	grep KEEP ./${norm}_call_stats.out  >> ./${norm}_call_stats.pass.out

        python /home/seck-thiam/1_metalfox.py -f1 ./${norm}_call_stats.pass.out -f2 ./${norm}_filtered_PASS.vcf -f3 ${params.inputDir}/${tumor_bam[0]} > ./${norm}_filter_stat.report
	cp ./${norm}_mutect1.snvs.pass.vcf ${params.vcf}
	cp ./${norm}_filtered_PASS.vcf ${params.vcf}
        cp ./${norm}_call_stats.pass.out ${params.vcf}
        cp ./${norm}_filter_stat.report ${params.vcf}

        """
 }	


//process MergeMutectStrelka{
//	module "picard/2.21.7"
//	errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
 //       maxRetries 2
   //     memory { 10.GB * task.attempt}

//	input:
//	tuple val(norm), file ("*_strelka.indels.pass.vcf") from strelka2_variants

//	output:
//	tuple val(norm), file("*_snv.indels.wes.variants.vcf") into all_variants_ch

//	script:
//	"""
//	java -jar \$PICARD MergeVcfs INPUT=${params.vcf}/${norm}_mutect1.snvs.pass.vcf INPUT=${params.vcf}/${norm}_strelka.indels.pass.vcf OUTPUT=${params.vcf}/${norm}_snv.indels.wes.variants.vcf
//	"""
//}
