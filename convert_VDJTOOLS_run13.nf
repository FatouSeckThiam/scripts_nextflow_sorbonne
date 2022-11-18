#!/usr/bin/env nextflow

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
params.samplelist = "/home/seck-thiam/samplelist_vdjtools_run_13"
params.inputDir = "/SCVOL02/run_13/rnaseq/mixcr"
params.outDir = "/SCVOL02/run_13/convert"
params.cpus = 4
params.hepl = false

def usage() {
	println("\nEtude du microenvironnement mixcr and vdjtools")
        println("\nINFO Pipeline parameters :")
        println("inputDir $params.inputDir = path fastq")
        println("outDir $params.outDir = path outDir")
        println("cpus $params.cpus = number of process")
}


if (params.help) {
        usage()
        exit(1)
} else if (params.samplelist == " "|| params.inputDir == " " || params.outDir == " " ) {
        println("\nUn ou plusieurs arguments réquis non renseignés, merci de suivre les instructions ci-dessus")
}

samplist = file("${params.samplelist}")

inputs = []

reader = samplist.newReader()
samplist.withReader {
    String line
    while( line = reader.readLine() ) {
        String id_mixcr = line.split("\t")[0]
        IGH = "${params.inputDir}/" + id_mixcr + "/"+ id_mixcr + "_rnaseq_.clonotypes.IGH.txt"
        IGK = "${params.inputDir}/" + id_mixcr + "/"+ id_mixcr + "_rnaseq_.clonotypes.IGK.txt"
        IGL = "${params.inputDir}/" + id_mixcr + "/"+ id_mixcr + "_rnaseq_.clonotypes.IGL.txt"
        TRA = "${params.inputDir}/" + id_mixcr + "/"+ id_mixcr + "_rnaseq_.clonotypes.TRA.txt"
        TRB = "${params.inputDir}/" + id_mixcr + "/"+ id_mixcr + "_rnaseq_.clonotypes.TRB.txt"
        TRD = "${params.inputDir}/" + id_mixcr + "/"+ id_mixcr + "_rnaseq_.clonotypes.TRD.txt"
        TRG = "${params.inputDir}/" + id_mixcr + "/"+ id_mixcr + "_rnaseq_.clonotypes.TRG.txt"
        inputs.add([IGH, IGK, IGL, TRA, TRB, TRD, TRG])
    }
liste_clonotype=["ALL","IGH","IGK", "IGL", "TRA", "TRB", "TRD","TRG"]
clones=Channel.fromList(inputs).flatten().view()

/*
process Mixcr {

        module "mixcr/3.0"
        time "96h"
        memory "40GB"

        input:
	tuple val(id), file(reads) from inputChannel

        output:
        path("${id}_.clonotypes.*.txt") into clones

        shell:
	"""
	mkdir -p !{params.outDir}/!{id}
        mkdir !{id}
        mixcr analyze shotgun --species hs  --starting-material rna --only-productive !{reads[0]} !{reads[1]}  ./!{id}_
        cp ./!{id}_.clonotypes.*.txt !{params.outDir}/!{id}
        """
 }

*/


Mettre Channel.flatern() pour eviter que tous les clonotypes sont émis de manière simultanée par le cannal résultat
#!/bin/bash -ue
java -Xmx16G -jar $VDJTOOLS Convert -S mixcr FR-01-173-RP_rnaseq_.clonotypes.ALL.txt FR-01-173-RP_rnaseq_.clonotypes.IGH.txt FR-01-173-RP_rnaseq_.clonotypes.IGK.txt FR-01-173-RP_rnaseq_.clonotypes.IGL.txt FR-01>
cp ./Convert_FR-01-173-RP_rnaseq_* /SCVOL02/Fatou/Analyse_repertoire/Test_mixr_VDJtools_nextflow/FR-01-173-RP_rnaseq
*/

/*
process convert_mixcr{
        module "vdjtools/1.2.1"
        time "96h"

        input:
	path(clonotypes) from clones

        output:
        path("Convert_${id}_*") into convert_mixcr

        shell:
	"""
	java -jar \$VDJTOOLS Convert -S mixcr !{clonotypes} ./Convert_!{id}_
        cp ./Convert_!{id}_* !{params.outDir}/!{id}
        """

}

convert_mixcr.into{ ch_PlotFancyVjusage; ch_PlotFancySpectratype; ch_Plotpectratype }

process PlotFancyVJUsage{

        module "vdjtools/1.2.1"
        time "96h"

        input:
	path(Convert_id) from ch_PlotFancyVjusage

        output:
        path("PlotFancyVJUsage_${id}_*")

        shell:
	"""
        java -jar \$VDJTOOLS PlotFancyVJUsage !{Convert_id} ./PlotFancyVJUsage_!{id}_
        cp PlotFancyVJUsage_!{id}_* !{params.outDir}/!{id}
        """
}

process PlotFancySpectratype  {

        module "vdjtools/1.2.1"
        time "96h"

        input:
	path(Convert_id) from ch_PlotFancySpectratype

        output:
        path("PlotFancySpectratype_${id}_*")

        shell:
	"""
        java -jar \$VDJTOOLS PlotFancySpectratype -top 10 !{Convert_id} ./PlotFancySpectratype_!{id}_
        cp PlotFancySpectratype_!{id}_* !{params.outDir}/!{id}
        """
}

process Plotpectratype  {
        module "vdjtools/1.2.1"

        input:
	path(Convert_id) from ch_Plotpectratype

        output:
        path("PlotSpectratypeV_${id}_*")

        shell:
	"""
	java -jar \$VDJTOOLS PlotSpectratypeV !{file} !{params.outDir}/PlotSpectratypeV_!{id}_
        cp PlotSpectratypeV_!{id}_* !{params.outDir}/!{id}
        """
}

