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
params.samplelist = "/home/seck-thiam/samplelist_cnv_run_13"
params.inputDir = "/SCVOL02/run_13/exome/Bam"
//rennomer bam en bam_files pour ne pas prendre en compte l'extension bam et donc éviter cette erreur //SCVOL01/Projects/Ideation/run_06/exome/bai

params.output="/SCVOL02/run_13/exome/cnv/"
params.errout="/SCVOL02/run_13/exome/cnv/errout"
//params.run_name="run_11"
params.chr_vcf_gz="/SCVOL01/Workspace/Labreche/Ideation/cnv/facets/00-common_all.chr.vcf.gz"
params.run_snp_pileup_facets_R="/SCVOL01/Workspace/Labreche/Ideation/cnv/facets/run_snp-pileup_facets2.R"

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
if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist == " "|| params.inputDir == " " || params.output == " " || params.errout  == " " || params.chr_vcf_gz  == " " || params.run_snp_pileup_facets_R == " "){
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

samples=file("${params.samplelist}")
// inputChannel = Channel.empty()

inputs = []

reader = samples.newReader()
samples.withReader {
    String line
    while( line = reader.readLine() ) {
        String normal_bn = line.split("\t")[0].split("_")[0].replaceAll(".bam","")
        String tumor_bn = line.split("\t")[1].split("_")[0].replaceAll(".bam","")
        String normal = line.split("\t")[0]
        normal_bam = "${params.inputDir}/" + normal //normal
        normal_idx = normal_bam.replaceAll("bam","bai")
        String tumoral = line.split("\t")[1]
        tumor_bam = "${params.inputDir}/" + tumoral
        tumor_idx = tumor_bam.replaceAll("bam","bai")
        inputs.add([normal_bn, [normal_bam, normal_idx], tumor_bn, [tumor_bam, tumor_idx]])
    }
}


//normalChannel = Channel.fromList(inputs_norm)//.subscribe{ println it }

InputChannel = Channel.fromList(inputs).view()


process run_cnv_plot{
	
	module "R/3.5.2:samtools:conda3/4.7.12"
	
	conda '/home/seck-thiam/my_env_conda.yaml'	

	cpus params.cpus	
	
	memory "20GB"

	input:
	tuple val(norm), path(normal_bam), val(tum), path(tumor_bam) from InputChannel
	
	output:
	tuple val(norm), path("./${norm}_snp-pileup.csv.gz") into cnv_ch

	shell:
	"""
	mkdir -p !{params.output}
	mkdir -p !{params.errout}

	snp-pileup -g -q20 -Q20 -r5 --min-read-counts 20 -x !{params.chr_vcf_gz} ./!{norm}_snp-pileup.csv !{params.inputDir}/!{normal_bam[0]} !{params.inputDir}/!{tumor_bam[0]}	

	cp ./!{norm}_snp-pileup.csv.gz !{params.output}
	"""
}

process run_facetS_R{
	
	module "R/3.5.2:samtools:conda3/4.7.12"
	
	input:
        tuple val(norm), path(snp_pileup_csv_gz) from cnv_ch
	
	output:
	tuple val(norm), path("./${norm}_pileup_facets.Rout") into facets_ch

	shell:
	"""
	R CMD BATCH '--args id="!{norm}" output="!{params.output}"' !{params.run_snp_pileup_facets_R} ./!{norm}_pileup_facets.Rout

	cp ./!{norm}_pileup_facets.Rout !{params.errout}

	"""
}	

//cp ./!{norm}_snpileup_facets_cnv.pdf !{params.output}	=> A NE PAS FAIRE CAR SEUL {norm}_pileup_facets.Rout ET	{norm}_snp-pileup.csv.gz SONT PROSUITS DANS LE WORK DE NEXTFLOW



