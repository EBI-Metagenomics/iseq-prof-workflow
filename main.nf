#!/usr/bin/env nextflow

params.store = "/hps/research/finn/horta/db"
params.outdir = "/Users/horta/code/iseq-profmark/result"
params.hmmfile = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz"
params.scriptdir = "$baseDir/script"
params.chunkSize = 50
params.seed = 98748
params.downsample_cds = 0
params.domains = "archaea:4,bacteria:10"

scriptdir = file(params.scriptdir)

Channel
    .fromList(params.domains.tokenize(","))
    .set { domain_specs_ch }

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

process download_pfam_hmm {
    publishDir params.outdir, mode:"copy"
    storeDir "$params.store/pfam"

    output:
    path "db.hmm" into hmmfile_ch

    script:
    """
    curl -s $params.hmmfile | gunzip -c > db.hmm
    """
}

process download_organism_names {
    input:
    val domain_spec from domain_specs_ch

    output:
    tuple path("${domain}.txt"), val(nsamples) into domain_files_spec_ch

    script:
    domain = domain_spec.tokenize(":")[0]
    nsamples = domain_spec.tokenize(":")[1]
    """
    $scriptdir/download_organism_names.py $domain
    """
}

process download_genbank_catalog {
    storeDir "$params.store/genbank"

    input:
    val db from Channel.fromList(["gss", "other"])

    output:
    path "gb238.catalog.${db}.tsv" into gb_catalog_ch1

    script:
    """
    $scriptdir/download_genbank_catalog.sh $db gb238.catalog.${db}.tsv
    """
}

process merge_genbank_catalogs {
    storeDir "$params.store/genbank"

    input:
    path "*" from gb_catalog_ch1.collect()

    output:
    path "gb238.catalog.all.tsv" into gb_catalog_ch2

    script:
    """
    $scriptdir/merge_genbank_catalog.sh *.tsv gb238.catalog.all.tsv
    """
}

process unique_genbank_organisms {
    memory "30 GB"
    storeDir "$params.store/genbank"

    input:
    path "gb238.catalog.all.tsv" from gb_catalog_ch2

    output:
    path "gb238.catalog.tsv" into gb_catalog_ch3

    script:
    """
    $scriptdir/unique_genbank_catalog.py gb238.catalog.all.tsv gb238.catalog.tsv
    """
}

process sample_accessions {
    input:
    path "gb238.catalog.tsv" from gb_catalog_ch3
    tuple path(domaintxt), val(nsamples) from domain_files_spec_ch

    output:
    path "accessions" into acc_file_ch

    script:
    domain = domaintxt.name.toString().tokenize(".")[0]
    """
    $scriptdir/sample_accessions.py gb238.catalog.tsv $domaintxt accessions $nsamples $params.seed
    """
}

acc_file_ch
    .splitText() { it.trim() }
    .into { acc_ch1; acc_ch2 }

process press_hmmfile {
    input:
    path hmmfile from hmmfile_ch

    output:
    path "*.hmm*", includeInputs: true into hmmdb_ch

    script:
    """
    hmmpress $hmmfile
    """
}

process download_genbank_gb {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }
    storeDir "$params.store/genbank"
    maxForks 1

    input:
    val acc from acc_ch1

    output:
    tuple val(acc), path("${acc}.gb") into genbank_gb_ch

    script:
    """
    $scriptdir/download_genbank.py $acc gb ${acc}.gb
    """
}

process download_genbank_fasta {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }
    storeDir "$params.store/genbank"
    maxForks 1

    input:
    val acc from acc_ch2

    output:
    tuple val(acc), path("${acc}.fasta") into genbank_fasta_ch

    script:
    """
    $scriptdir/download_genbank.py $acc fasta ${acc}.fasta
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
    if [[ "$params.downsample_cds" != "0" ]];
    then
       $scriptdir/downsample_fasta.py cds_amino.fasta .cds_amino.fasta $params.downsample_cds
       mv .cds_amino.fasta cds_amino.fasta
       $scriptdir/downsample_fasta.py cds_nucl.fasta .cds_nucl.fasta $params.downsample_cds
       mv .cds_nucl.fasta cds_nucl.fasta
    fi
    """
}

process hmmscan {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }
    cpus 4

    input:
    path hmmdb from hmmdb_ch.collect()
    tuple val(acc), path(amino) from cds_amino_ch

    output:
    tuple val(acc), path("domtblout.txt") into hmmscan_output_ch

    script:
    """
    hmmfile=\$(ls *.hmm)
    hmmscan --noali --cut_ga --domtblout domtblout.txt --cpu ${task.cpus} \$hmmfile $amino
    """
}

process create_true_false_profiles {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    path hmmdb from hmmdb_ch.collect()
    tuple val(acc), path("domtblout.txt") from hmmscan_output_ch

    output:
    tuple val(acc), path("*.hmm*") into dbspace_ch

    script:
    """
    hmmfile=\$(ls *.hmm)
    $scriptdir/create_true_false_profiles.py domtblout.txt \$hmmfile accspace.txt dbspace.hmm $params.seed
    hmmfetch --index dbspace.hmm
    """
}

cds_nucl_ch
    .join(dbspace_ch)
    .splitFasta(by:params.chunkSize, file:true, elem:1)
    .set { cds_nucl_db_split_ch }

process iseq_scan {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/chunks/$name" }

    errorStrategy "retry"
    maxRetries 2
    memory "8 GB"

    input:
    tuple val(acc), path(nucl), path(dbspace) from cds_nucl_db_split_ch

    output:
    tuple val(acc), path("output.*.gff") into iseq_output_split_ch
    tuple val(acc), path("oamino.*.fasta") into iseq_oamino_split_ch
    tuple val(acc), path("ocodon.*.fasta") into iseq_ocodon_split_ch

    script:
    chunk = nucl.name.toString().tokenize('.')[-2]
    """
    hmmfile=\$(ls *.hmm)
    iseq pscan3 \$hmmfile $nucl --hit-prefix chunk_${chunk}_item\
        --output output.${chunk}.gff --oamino oamino.${chunk}.fasta\
        --ocodon ocodon.${chunk}.fasta\
        --no-cut-ga --no-heuristic
    """
}

iseq_output_split_ch
    .collectFile(keepHeader:true, skip:1)
    .map { ["output.gff", it] }
    .set { iseq_output_ch }

iseq_oamino_split_ch
    .collectFile()
    .map { ["oamino.fasta", it] }
    .set { iseq_oamino_ch }

iseq_ocodon_split_ch
    .collectFile()
    .map { ["ocodon.fasta", it] }
    .set { iseq_ocodon_ch }

iseq_output_ch
    .mix(iseq_oamino_ch, iseq_ocodon_ch)
    .set { iseq_results_ch }

process save_output {
    publishDir params.outdir, mode:"copy", saveAs: { name -> "${acc}/${name}" }

    input:
    tuple val(name), path(acc) from iseq_results_ch

    output:
    path(name)

    script:
    """
    mv $acc $name
    """
}
