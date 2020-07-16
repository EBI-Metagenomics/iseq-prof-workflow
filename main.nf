#!/usr/bin/env nextflow

params.outdir = "/Users/horta/code/iseq-profmark/result"
params.accfile = "$baseDir/accessions.txt"
params.hmmfile = "/Users/horta/code/iseq-profmark/Pfam-A_small.hmm"
params.scriptdir = "$baseDir/script"
params.chunkSize = 4

Channel
    .fromList(file(params.accfile).readLines())
    .into { acc_ch1; acc_ch2 }

scriptdir = file(params.scriptdir)

process download_hmmfile {
    input:
    path hmmfile from params.hmmfile

    output:
    path "db.hmm" into hmmfile_ch1, hmmfile_ch2

    script:
    """
    cp $hmmfile db.hmm
    """
}

process press_hmmfile {
    input:
    path hmmfile from hmmfile_ch1

    output:
    path "db.*", includeInputs: true into hmmdb_ch1, hmmdb_ch2

    script:
    """
    hmmpress $hmmfile
    """
}

process download_genbank_gb {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    val acc from acc_ch1

    output:
    tuple val(acc), path("genbank.gb") into genbank_gb_ch

    script:
    """
    $scriptdir/download_genbank.py $acc gb genbank.gb
    """
}

process download_genbank_fasta {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    val acc from acc_ch2

    output:
    tuple val(acc), path("genbank.fasta") into genbank_fasta_ch

    script:
    """
    $scriptdir/download_genbank.py $acc fasta genbank.fasta
    """
}

process extract_cds {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    tuple val(acc), path(gb) from genbank_gb_ch

    output:
    tuple val(acc), path("cds_amino.fasta") into cds_amino_ch
    tuple val(acc), path("cds_nucl.fasta") into cds_nucl_ch1, cds_nucl_ch2

    script:
    """
    $scriptdir/extract_cds.py $gb cds_amino.fasta cds_nucl.fasta
    """
}

process hmmscan {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    path hmmdb from hmmdb_ch1
    tuple val(acc), path(amino) from cds_amino_ch

    output:
    tuple val(acc), path("domtblout.txt") into hmmscan_output_ch

    script:
    """
    hmmscan --cut_ga --domtblout domtblout.txt db.hmm $amino
    """
}

cds_nucl_ch1
    .splitFasta(by: params.chunkSize, file:true, elem:1)
    .set { cds_nucl_split_ch }

process iseq_scan {
    input:
    path hmmdb from hmmdb_ch2
    tuple val(acc), path(nucl) from cds_nucl_split_ch

    output:
    path "${acc}_output.gff" into iseq_output_split_ch
    path "${acc}_oamino.fasta" into iseq_oamino_split_ch
    path "${acc}_ocodon.fasta" into iseq_ocodon_split_ch

    script:
    chunk = nucl.name.toString().tokenize('.')[-2]
    """
    out=""
    iseq pscan2 db.hmm $nucl --hit-prefix chunk_${chunk}_item\
        --output ${acc}_output.gff --oamino ${acc}_oamino.fasta\
        --ocodon ${acc}_ocodon.fasta --quiet
    """
}

iseq_output_split_ch
   .collectFile(keepHeader:true, skip:1)
   .map { [it.name.toString().tokenize('_')[0], it] }
   .set { iseq_output_ch1 }

iseq_oamino_split_ch
   .collectFile()
   .map { [it.name.toString().tokenize('_')[0], it] }
   .set { iseq_oamino_ch }

iseq_ocodon_split_ch
   .collectFile()
   .map { [it.name.toString().tokenize('_')[0], it] }
   .set { iseq_ocodon_ch }

iseq_output_ch1
   .join(iseq_oamino_ch)
   .join(iseq_ocodon_ch)
   .set { iseq_results_ch }

process save_iseq_results {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    tuple val(acc), path(output), path(oamino), path(ocodon) from iseq_results_ch

    output:
    tuple val(acc), path("output.gff") into iseq_output_ch2
    path("oamino.fasta")
    path("ocodon.fasta")

    script:
    """
    mv $output output.gff
    mv $oamino oamino.fasta
    mv $ocodon ocodon.fasta
    """
}

cds_nucl_ch2
   .join(hmmscan_output_ch)
   .join(iseq_output_ch2)
   .set { results_ch }

process profmark {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    path hmmfile from hmmfile_ch2
    tuple val(acc), path(nuclfile), path("domtblout.txt"), path("output.gff") from results_ch

    output:
    path("profmark.pkl")

    script:
    """
    $scriptdir/evaluate.py $hmmfile $nuclfile domtblout.txt output.gff profmark.pkl
    """
}
