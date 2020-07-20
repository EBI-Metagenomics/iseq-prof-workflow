#!/usr/bin/env nextflow

params.outdir = "/Users/horta/code/iseq-profmark/result"
params.accfile = "$baseDir/accessions.txt"
params.hmmfile = "https://iseq-py.s3.eu-west-2.amazonaws.com/Pfam-A_24.hmm"
params.scriptdir = "$baseDir/script"
params.chunkSize = 100

process save_params {
    publishDir params.outdir, mode:"copy"

    output:
    path "params.txt"

    script:
    params_cmd = params.all()
    """
    echo $params_cmd > params.txt
    """
}

process download_organism_names
{
   output:
   tuple path("archaea.txt"), path("bacteria.txt") into organisms_ch

   script:
   """
   $basedir/download_organism_names.py
   file archaea.txt
   file bacteria.txt
   """
}

process download_genbank_catalog
{
    memory "12 GB"

    output:
    path "gb238.catalog.all.tsv" into gb_catalog_ch1

    script:
    """
    new_header="Version\tMolType\tBasePairs\tOrganism\tTaxID"
    old_header="Accession\tVersion\tID\tMolType\tBasePairs\tOrganism\tTaxID\tDB\tBioProject\tBioSample"

    echo -e $new_header > gb238.catalog.all.tsv

    for db in est gss other;
    do
       echo $db
       curl -s ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/gb238.catalog.${db}.txt.gz \
          | gunzip -c  \
          | cut -d$'\t' -f2,4,5,6,7 \
          | grep --invert-match $'\tmRNA\t' \
          | grep --invert-match $'\trRNA\t' \
          | grep --invert-match $'\ttRNA\t' \
          | grep --invert-match $'\tncRNA\t' \
          | grep --invert-match $'\tNoTaxID' \
          >> gb238.catalog.all.tsv
    done
    """
}

process unique_genbank_organisms
{
    memory "60 GB"

    input:
    path "gb238.catalog.all.tsv" from gb_catalog_ch1

    output:
    path "gb238.catalog.unique.feather" into gb_catalog_ch2

    script:
    """
    $basedir/unique_genbank_catalog.py gb238.catalog.unique.feather gb238.catalog.unique.feather
    """
}

Channel
    .fromList(file(params.accfile).readLines())
    .into { acc_ch1; acc_ch2 }

scriptdir = file(params.scriptdir)
hmmfile_ch = Channel.value(file(params.hmmfile))

process copy_hmmfile {
    publishDir params.outdir, mode:"copy"

    input:
    path hmmfile from hmmfile_ch

    output:
    path "*.hmm", includeInputs: true

    script:
    """
    """
}

process press_hmmfile {
    input:
    path hmmfile from hmmfile_ch

    output:
    path "*.hmm*", includeInputs: true into hmmdb_ch1, hmmdb_ch2

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
    tuple val(acc), path("cds_nucl.fasta") into cds_nucl_ch

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
    hmmfile=\$(ls *.hmm)
    hmmscan --cut_ga --domtblout domtblout.txt \$hmmfile $amino
    """
}

cds_nucl_ch
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
    hmmfile=\$(ls *.hmm)
    iseq pscan2 \$hmmfile $nucl --hit-prefix chunk_${chunk}_item\
        --output ${acc}_output.gff --oamino ${acc}_oamino.fasta\
        --ocodon ${acc}_ocodon.fasta --quiet
    """
}

iseq_output_split_ch
   .collectFile(keepHeader:true, skip:1)
   .map { [it.name.toString().tokenize('_')[0], it] }
   .set { iseq_output_ch }

iseq_oamino_split_ch
   .collectFile()
   .map { [it.name.toString().tokenize('_')[0], it] }
   .set { iseq_oamino_ch }

iseq_ocodon_split_ch
   .collectFile()
   .map { [it.name.toString().tokenize('_')[0], it] }
   .set { iseq_ocodon_ch }

iseq_output_ch
   .join(iseq_oamino_ch)
   .join(iseq_ocodon_ch)
   .set { iseq_results_ch }

process save_iseq_results {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    tuple val(acc), path(output), path(oamino), path(ocodon) from iseq_results_ch

    output:
    path("output.gff")
    path("oamino.fasta")
    path("ocodon.fasta")

    script:
    """
    mv $output output.gff
    mv $oamino oamino.fasta
    mv $ocodon ocodon.fasta
    """
}
