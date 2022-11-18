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
//params.samples = "/home/seck-thiam/samplelist_vcf_neo_run_09"
params.samples = "/home/seck-thiam/samplelist_vcf_neo_run_13"

params.inputDir_vcf_converted = "/SCVOL02/run_13/exome/neoepitope_vcf"
params.outDir = "/SCVOL02/run_13/neoepitope_vcf/TX"
params.vcf_radia ="/SCVOL02/run_13/rnaseq/radia"
params.input_Bam_rnaseq="/SCVOL02/run_13/rnaseq/Bam"
params.input_kallisto="/SCVOL02/run_13/rnaseq/expression"
params.refGenome = "/home/seck-thiam/Homo_sapiens.hg38.fa"
params.VEP_plugins="/SCVOL02/VEP_plugins"
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
} else if (params.samples == " "|| params.inputDir_vcf_converted == " " || params.outDir == " " || params.vcf_radia == " " || params.input_Bam_rnaseq == " " || params.input_kallisto == " " || params.refGenome == " ") {
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

samplelist=file("${params.samples}")

inputs=[]

reader=samplelist.newReader()
samplelist.withReader() {
	String line
	while( line=reader.readLine() ){
	
	String id_mutect1 = line.split("\t")[1].split("_")[0].replaceAll(".vcf","")
	String id_strelka2 = line.split("\t")[2].split("_")[0].replaceAll(".vcf","")
	String mutect1 = line.split("\t")[1]
	String strelka2 = line.split("\t")[2]
	vcf_mutect1 ="${params.inputDir_vcf_converted}/"+mutect1
	vcf_strelka2 = "${params.inputDir_vcf_converted}/"+strelka2
	String id_radia= line.split("\t")[3].split("_")[0].replaceAll(".vcf","")
	radia=line.split("\t")[3]
	radia_vcf = "${params.vcf_radia}/"+radia
	String id_rnaseq = line.split("\t")[4].split("_")[0].replaceAll(".bam","")
	rnaseq=line.split("\t")[4]
	rnaseq_bam= "${params.input_Bam_rnaseq}/"+rnaseq
	rnaseq_index = rnaseq_bam.replaceAll(".bam",".bam.bai")
	kallisto_id = line.split("\t")[0]
	kallisto_res = "${params.input_kallisto}/"+kallisto_id 
	inputs.add([id_mutect1,[vcf_mutect1], id_strelka2,[vcf_strelka2], id_radia,[radia_vcf], id_rnaseq,[rnaseq_bam,rnaseq_index], kallisto_id,[kallisto_res]])
  }
 
}

vcf_converted_mutect1_strelka2_channel=Channel.fromList(inputs).view()
vcf_converted_mutect1=Channel.fromList(inputs)


process merge_vcf_converted {

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	memory "40GB"

	input:
	tuple val(id), path(vcf_mutect1), val(id_s), path(vcf_strelka2), val(id_ra), path(radia_vcf), val(id_rna), path(rnaseq_bam), val(kal), path(kallisto_res) from vcf_converted_mutect1_strelka2_channel	
	
	output:
	tuple val(id), path("${id}_mutect1_strelka.pass.converted.vcf"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) into strelk2_mutect1_variants
	
	shell:
	"""
	mkdir -p !{params.outDir}

	java -jar \$PICARD MergeVcfs INPUT=!{params.inputDir_vcf_converted}/!{vcf_mutect1} INPUT=!{params.inputDir_vcf_converted}/!{vcf_strelka2} OUTPUT=./!{id}_mutect1_strelka.pass.converted.vcf
	
	cp ./!{id}_mutect1_strelka.pass.converted.vcf !{params.outDir}	

	"""
}

process create_vcf_radia{
	
	input:
	tuple val(id), path(mutect1_strelka_pass_converted), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) from strelk2_mutect1_variants

	output:
	tuple val(id), path("${id}_radia_id.vcf"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) into radia_vcf_channel

	shell:
	"""
	grep -w "#" !{mutect1_strelka_pass_converted} > ./!{id}_radia_id.vcf

	grep "#CHROM" !{params.vcf_radia}/!{radia_vcf} | sed -e "s/DNA_NORMAL/!{id}_Germline/g" | sed -e "s/DNA_TUMOR/!{id}_Tumour/g" | sed -e "s/RNA_TUMOR/!{id}_Tumour_RNA/g"  >> ./!{id}_radia_id.vcf
	
	grep "PASS" !{params.vcf_radia}/!{radia_vcf} | grep -v "RADIAGerm" | grep -v "RADIASomDNA" >> ./!{id}_radia_id.vcf

	cp ./!{id}_radia_id.vcf !{params.outDir}

	"""
}

process run_intersectBed{

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"

	input: 
	tuple val(id), path(radia_id_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) from radia_vcf_channel

	output:
	tuple val(id), path("${id}_wes_rna_snv_intersection.vcf"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) into intersectBed_channel
		
	shell:
	
	"""
	intersectBed -header -wa -a !{vcf_mutect1} -b !{radia_id_vcf} > ./!{id}_wes_rna_snv_intersection.vcf
	cp !{id}_wes_rna_snv_intersection.vcf !{params.outDir}
	"""
}


process mergeVCF{

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	
	input:
	tuple val(id), path(wes_rna_snv_intersection_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) from intersectBed_channel
	
	output:
	tuple val(id), path("${id}_wes_rna_all_intersection.vcf"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) into mergeVCF_channel
	
	shell:
	"""
	java -jar \$PICARD MergeVcfs INPUT=!{wes_rna_snv_intersection_vcf} INPUT=!{vcf_strelka2} OUTPUT=./!{id}_wes_rna_all_intersection.vcf
	
	cp !{id}_wes_rna_all_intersection.vcf !{params.outDir}

	"""
}
process run_Vep{

	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	
	input:
	tuple val(id), path(wes_rna_all_intersection_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) from mergeVCF_channel
	
	output:
	tuple val(id), path("${id}_wes_rna_intersection.vep.vcf"), path("${id}_sites.txt"), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) into Vep_channel

	shell:
	"""
	vep --input_file !{wes_rna_all_intersection_vcf} --output_file ./!{id}_wes_rna_intersection.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta !{params.refGenome} --offline --cache --dir_cache /SCVOL01/Tools/ensembl-vep-99.0/.vep ensembl --assembly GRCh38 --plugin Frameshift --plugin Wildtype --dir_plugins !{params.VEP_plugins} --transcript_version

	grep -v "#" ./!{wes_rna_all_intersection_vcf} | awk '{print \$1,"\t", \$2-10,"\t", \$2+10}' > ./!{id}_sites.txt
	
	cp !{id}_wes_rna_intersection.vep.vcf !{params.outDir}
	cp !{id}_sites.txt !{params.outDir}
	"""
}

process Bam_readcount{
	module "bedtools:vep/release-99:picard/2.21.7:bam-readcount"
	
	input:
	tuple val(id), path(wes_rna_intersection_vep_vcf), path(sites_txt), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) from Vep_channel
	
	output:
	tuple val(id), path("${id}_sites_readcount.vcf"), path(wes_rna_intersection_vep_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) into Bam_readcount_channel
	
	shell:
	"""
	bam-readcount -f !{params.refGenome} !{params.input_Bam_rnaseq}/!{rnaseq_bam[0]} -l ./!{id}_sites.txt  > ./!{id}_sites_readcount.vcf
	cp !{id}_sites_readcount.vcf !{params.outDir}
 	"""

}
	
process Vcf_expression_annotator{

	module "bedtools:conda3/4.7.12:vep/release-99:picard/2.21.7:bam-readcount"
	//conda '/home/seck-thiam/my_env_conda_pvactools.yaml'
	
	input: 
	tuple val(id), path(sites_readcount_vcf), path(wes_rna_intersection_vep_vcf), path(vcf_mutect1), path(vcf_strelka2), path(radia_vcf), path(rnaseq_bam), path(kallisto_res) from Bam_readcount_channel

	output:
	tuple val(id), path("${id}_wes_rna_intersection.vep.tx.vcf"), path(sites_readcount_vcf) into vcf_expression_annotator_channel

	shell:
	"""
	sed -i 's/ID=FA/ID=AF/' !{wes_rna_intersection_vep_vcf}
	conda run -p /SCVOL01/Tools/pVACtools-1.5.4 vcf-expression-annotator -i target_id -o ./!{id}_wes_rna_intersection.vep.tx.vcf -s !{id}_Tumour -e tmp --ignore-ensembl-id-version !{wes_rna_intersection_vep_vcf} !{params.input_kallisto}/!{id}/abundance.tsv kallisto transcript
	
	cp !{id}_wes_rna_intersection.vep.tx.vcf !{params.outDir}
	"""
}

process Vcf_readcount_annotator{

	module "bedtools4:conda3/4.7.12:vep/release-99:picard/2.21.7:bam-readcount"
	//conda '/home/seck-thiam/my_env_conda_pvactools.yaml'
	
	input:
	tuple val(id), path(wes_rna_intersection_vep_tx_vcf), path(sites_readcount_vcf) from vcf_expression_annotator_channel

	output:
	tuple val(id), path("${id}_wes_rna_intersection.vep.tx.raf.vcf") into vcf_readcount_annotator_channel

	shell:

	"""
	conda run -p /SCVOL01/Tools/pVACtools-1.5.4 vcf-readcount-annotator -s !{id}_Tumour -o !{id}_wes_rna_intersection.vep.tx.raf.vcf -t all !{wes_rna_intersection_vep_tx_vcf} !{sites_readcount_vcf} RNA
	cp !{id}_wes_rna_intersection.vep.tx.raf.vcf !{params.outDir}
	
	"""
}







