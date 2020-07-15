#!/usr/bin/env nextflow

params.accfile = "$baseDir/accessions.txt"
params.hmmfile = "/Users/horta/code/iseq-profmark/Pfam-A_small.hmm"
params.scriptdir = "$baseDir/script"

Channel
    .fromList(file(params.accfile).readLines())
    .into { acc_ch1; acc_ch2 }

scriptdir = file(params.scriptdir)

process download_hmmfile {
    input:
    path hmmfile from params.hmmfile

    output:
    path "db.hmm" into hmmfile_ch

    script:
    """
    cp $hmmfile db.hmm
    """
}

process press_hmmfile {
    input:
    path hmmfile from hmmfile_ch

    output:
    path "db.*", includeInputs: true into hmmdb_ch1, hmmdb_ch2

    script:
    """
    hmmpress $hmmfile
    """
}

process download_genbank_gb {
    input:
    val acc from acc_ch1

    output:
    tuple val("$acc"), path("${acc}.gb") into genbank_gb_ch

    script:
    """
    $scriptdir/download_genbank.py $acc gb ${acc}.gb
    """
}

process download_genbank_fasta {
    input:
    val acc from acc_ch2

    output:
    tuple val("$acc"), path("${acc}.fasta") into genbank_fasta_ch

    script:
    """
    $scriptdir/download_genbank.py $acc fasta ${acc}.fasta
    """
}

process extract_cds {
    input:
    tuple acc, path(gb) from genbank_gb_ch

    output:
    path "${acc}_cds_amino.fasta" into cds_amino_ch
    path "${acc}_cds_nucl.fasta" into cds_nucl_ch

    script:
    """
    $scriptdir/extract_cds.py $gb ${acc}_cds_amino.fasta ${acc}_cds_nucl.fasta
    """
}

process hmmscan {
    input:
    path hmmdb from hmmdb_ch1
    path amino from cds_amino_ch

    output:
    path "${acc}_domtblout.txt" into hmmscan_output_ch

    script:
    acc = amino.name.toString().tokenize('_').get(0)
    """
    hmmscan --cut_ga --domtblout ${acc}_domtblout.txt db.hmm $amino
    """
}

process iseq_scan {
    input:
    path hmmdb from hmmdb_ch2
    path nucl from cds_nucl_ch

    output:
    path "${acc}_output.gff" into iseq_output_ch

    script:
    acc = nucl.name.toString().tokenize('_').get(0)
    """
    iseq pscan2 db.hmm $nucl --output ${acc}_output.gff
    """
}
