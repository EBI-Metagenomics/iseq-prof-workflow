#!/usr/bin/env nextflow

params.accfile = "$baseDir/accessions.txt"
// params.datadir = "$baseDir/data"
params.scriptdir = "$baseDir/script"

accessions = Channel.from(file(params.accfile).readLines())
// datadir = file(params.datadir)
scriptdir = file(params.scriptdir)

//datadir.mkdirs()
//datadir.isDirectory

//accessions = Channel
//  .fromPath(params.accfile)
//  .readLines()
//  .each { println it.trim() }

// acc_ch = Channel.fromPath(params.accfile)
// .splitFasta(by: params.chunkSize, file:true)
//    .set { fasta_ch }
//
//    myFile.eachLine {  str ->
//        println "line ${count++}: $str"
//    }

//process myProc {
//    input:
//    file acc from acc_ch
//
//    output:
//    stdout result
//
//    """
//    echo $acc
//    """
//}

// process myTask {
//
//     input:
//     val str from 'Hello', 'Hola', 'Bonjour'
//
//     shell:
//     '''
//     echo User $USER says !{str}
//     '''
//
// }

//Channel
//    .fromList( ['a', 'b', 'c', 'd'] )
//    .view { "value: $it" }

process download_genbank {
    input:
    val acc from accessions

    output:
    tuple path("${acc}.gb"), path("${acc}.fasta") into genbank_ch

    """
    $scriptdir/download_genbank.py $acc gb
    $scriptdir/download_genbank.py $acc fasta
    """
}

process checksum {
    input:
    tuple path(gb), path(fasta) from genbank_ch

    output:
    stdout result

    """
    sha256sum $gb
    sha256sum $fasta
    """
}

result.view { "$it" }
