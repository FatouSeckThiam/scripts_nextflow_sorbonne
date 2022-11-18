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



params.samplelist = "/SCVOL02/Fatou/samplelist_vdjtools_run_13"
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
}
liste_clonotype=["ALL","IGH","IGK", "IGL", "TRA", "TRB", "TRD","TRG"]
clones=Channel.fromList(inputs).flatten().view()

/*
Mettre Channel.flatern() pour eviter que tous les clonotypes sont émis de manière simultanée par le cannal résultat
java -Xmx16G -jar $VDJTOOLS Convert -S mixcr FR-01-173-RP_rnaseq_.clonotypes.ALL.txt FR-01-173-RP_rnaseq_.clonotypes.IGH.txt FR-01-173-RP_rnaseq_.clonotypes.IGK.txt FR-01-173-RP_rnaseq_.clonotypes.IGL.t>
cp ./Convert_FR-01-173-RP_rnaseq_* /SCVOL02/Fatou/Analyse_repertoire/Test_mixr_VDJtools_nextflow/FR-01-173-RP_rnaseq
*/


process convert_mixcr{
        module "vdjtools/1.2.1"
        time "96h"

        input:
	path(clonotypes) from clones

        output:
        path("Convert_${clonotypes}_*") into convert_mixcr

        shell:
	"""
        java -jar \$VDJTOOLS Convert -S mixcr !{clonotypes} ./Convert_!{Convert_clonotypes}_
        cp ./Convert_!{clonotypes}_* !{params.outDir}/
        """

}


convert_mixcr.into{ ch_PlotFancyVjusage; ch_PlotFancySpectratype; ch_Plotpectratype }

process PlotFancyVJUsage{

        module "vdjtools/1.2.1"
        time "96h"

        input:
	path(Convert_clonotypes) from ch_PlotFancyVjusage

        output:
        path("PlotFancyVJUsage_${clonotypes}_*")

        shell:
	"""
	java -jar \$VDJTOOLS PlotFancyVJUsage !{Convert_clonotypes} ./PlotFancyVJUsage_!{Convert_clonotypes}_
        cp PlotFancyVJUsage_!{Convert_clonotypes}_* !{params.outDir}
        """
}

process PlotFancySpectratype{

        module "vdjtools/1.2.1"
        time "96h"

        input:
	path(Convert_clonotypes) from ch_PlotFancySpectratype

        output:
        path("PlotFancySpectratype_${clonotypes}_*")

        shell:
	"""
	java -jar \$VDJTOOLS PlotFancySpectratype -top 10 !{Convert_clonotypes} ./PlotFancySpectratype_!{Convert_clonotypes}_
        cp PlotFancySpectratype_!{Convert_clonotypes}_* !{params.outDir}
        """
}

process Plotpectratype{ 

        module "vdjtools/1.2.1"

        input:
	path(Convert_clonotypes) from ch_Plotpectratype

        output:
        path("PlotSpectratypeV_${clonotypes}_*")

        shell:
	"""
	java -jar \$VDJTOOLS PlotSpectratypeV !{Convert_clonotypes} ./PlotSpectratypeV_!{Convert_clonotypes}_
        cp PlotSpectratypeV_!{Convert_clonotypes}_* !{params.outDir}
        """
}



			
