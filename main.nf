#!/usr/bin/env nextflow

println " "
println "Nextflow pipeline"
println " "
println "using profile: $workflow.profile"
println "using params: $params"
println " "

scriptDir = file("$projectDir/script")
groupRoot = "/horta/$workflow.runName"

Channel
    .fromList(params.domains.tokenize(","))
    .set { domain_specs_ch }

process save_params {
    publishDir params.outputDir, mode:"copy"
    clusterOptions "-g $groupRoot/save_params"

    output:
    path "params.txt"

    script:
    params_cmd = params.all()
    """
    echo "$params_cmd" > params.txt
    """
}

process download_pfam_hmm {
    clusterOptions "-g $groupRoot/download_pfam_hmm"
    storeDir "$params.storageDir/pfam"

    output:
    path "db.hmm" into hmmfile_pre_filter_ch

    script:
    """
    curl -s $params.hmmfile | gunzip -c > db.hmm
    """
}

process download_pfam_clans {
    clusterOptions "-g $groupRoot/download_pfam_clans"
    publishDir params.outputDir, mode:"copy"
    storeDir "$params.storageDir/pfam"

    output:
    path "clans.sto.gz" into clans_sto_ch

    script:
    """
    curl -s $params.clanfile > clans.sto.gz
    """
}

process create_clan_profile_assoc {
    clusterOptions "-g $groupRoot/create_clan_profile_assoc"
    publishDir params.outputDir, mode:"copy"
    storeDir "$params.storageDir/pfam"

    input:
    path "clans.sto.gz" from clans_sto_ch

    output:
    path "clans.csv" into clans_csv_ch

    script:
    """
    iseq-prof compute-clans clans.sto.gz clans.csv
    """
}

process filter_pfam_hmm {
    clusterOptions "-g $groupRoot/filter_pfam_hmm"
    memory "16 GB"
    publishDir params.outputDir, mode:"copy"

    input:
    path "db.hmm.in" from hmmfile_pre_filter_ch
    path "clans.csv" from clans_csv_ch

    output:
    path "db.hmm" into hmmfile_ch

    script:
    """
    iseq-prof hmm-filter db.hmm.in clans.csv db.hmm --by-clan "$params.filterClan" --no-case-sensitive
    """
}

process download_organism_names {
    clusterOptions "-g $groupRoot/download_organism_names"
    storeDir "$params.storageDir/genbank"

    input:
    val domain_spec from domain_specs_ch

    output:
    tuple val(domain), val(nsamples), path("${domain}.txt") into domain_files_spec_ch

    script:
    domain = domain_spec.tokenize(":")[0]
    nsamples = domain_spec.tokenize(":")[1]
    """
    $scriptDir/download_organism_names.py $domain
    """
}

process download_genbank_catalog {
    clusterOptions "-g $groupRoot/download_genbank_catalog"
    storeDir "$params.storageDir/genbank"

    input:
    val db from Channel.fromList(["gss", "other"])

    output:
    path "gb238.catalog.${db}.tsv" into gb_catalog_ch1

    script:
    """
    $scriptDir/download_genbank_catalog.sh $db gb238.catalog.${db}.tsv
    """
}

process merge_genbank_catalogs {
    clusterOptions "-g $groupRoot/merge_genbank_catalogs"
    storeDir "$params.storageDir/genbank"

    input:
    path "*" from gb_catalog_ch1.collect()

    output:
    path "gb238.catalog.all.tsv" into gb_catalog_ch2

    script:
    """
    $scriptDir/merge_genbank_catalog.sh *.tsv gb238.catalog.all.tsv
    """
}

process unique_genbank_organisms {
    clusterOptions "-g $groupRoot/unique_genbank_organisms"
    memory "30 GB"
    storeDir "$params.storageDir/genbank"

    input:
    path "gb238.catalog.all.tsv" from gb_catalog_ch2

    output:
    path "gb238.catalog.tsv" into gb_catalog_ch3

    script:
    """
    $scriptDir/unique_genbank_catalog.py gb238.catalog.all.tsv gb238.catalog.tsv
    """
}

process sample_accessions {
    clusterOptions "-g $groupRoot/sample_accessions"
    errorStrategy "retry"
    maxForks 1
    maxRetries 2
    storeDir "$params.storageDir/genbank"

    input:
    path "gb238.catalog.tsv" from gb_catalog_ch3
    tuple val(domain), val(nsamples), path(domaintxt) from domain_files_spec_ch

    output:
    path "${domain}_${nsamples}_accessions" into acc_file_ch

    script:
    """
    $scriptDir/sample_accessions.py gb238.catalog.tsv $domaintxt ${domain}_${nsamples}_accessions $nsamples $params.seed
    """
}

acc_file_ch
    .splitText() { it.trim() }
    .filter ( ~"$params.filterAcc" )
    .into { acc_ch1; acc_ch2 }

filterClanHash = params.filterClan.digest('SHA-256')
filterClanHash = filterClanHash[0..3] + filterClanHash[4..7]

process press_hmmfile {
    clusterOptions "-g $groupRoot/press_hmmfile"
    storeDir "$params.storageDir/pfam/$filterClanHash/press"

    input:
    path hmmfile from hmmfile_ch

    output:
    path "*.hmm*", includeInputs: true into hmmdb_ch

    script:
    """
    if [ -s $hmmfile ]
    then
        hmmpress $hmmfile
    fi
    """
}

process download_genbank_gb {
    clusterOptions "-g $groupRoot/download_genbank_gb"
    errorStrategy "retry"
    maxForks 1
    maxRetries 2
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/$name" }
    storeDir "$params.storageDir/genbank"

    input:
    val acc from acc_ch1

    output:
    tuple val(acc), path("${acc}.gb") into genbank_gb_ch

    script:
    """
    $scriptDir/download_genbank.py $acc gb ${acc}.gb
    """
}

process download_genbank_fasta {
    clusterOptions "-g $groupRoot/download_genbank_fasta"
    errorStrategy "retry"
    maxForks 1
    maxRetries 2
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/$name" }
    storeDir "$params.storageDir/genbank"

    input:
    val acc from acc_ch2

    output:
    tuple val(acc), path("${acc}.fasta") into genbank_fasta_ch

    script:
    """
    $scriptDir/download_genbank.py $acc fasta ${acc}.fasta
    """
}

process extract_cds {
    clusterOptions "-g $groupRoot/extract_cds"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    tuple val(acc), path(gb) from genbank_gb_ch

    output:
    tuple val(acc), path("cds_amino.fasta") into cds_amino_ch
    tuple val(acc), path("cds_nucl.fasta") into cds_nucl_ch

    script:
    """
    $scriptDir/extract_cds.py $gb cds_amino.fasta cds_nucl.fasta
    if [[ "$params.downsampleCDS" != "0" ]];
    then
       $scriptDir/downsample_fasta.py cds_amino.fasta .cds_amino.fasta $params.downsampleCDS
       mv .cds_amino.fasta cds_amino.fasta
       $scriptDir/downsample_fasta.py cds_nucl.fasta .cds_nucl.fasta $params.downsampleCDS
       mv .cds_nucl.fasta cds_nucl.fasta
    fi
    """
}

process hmmscan {
    clusterOptions "-g $groupRoot/hmmscan -R 'rusage[scratch=5120]'"
    cpus 4
    memory "8 GB"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/$name" }
    scratch true
    stageInMode "copy"

    input:
    path hmmdb from hmmdb_ch.collect()
    tuple val(acc), path(amino) from cds_amino_ch

    output:
    tuple val(acc), path("domtblout.txt") into hmmscan_output_ch

    script:
    """
    hmmfile=\$(echo *.hmm)
    if [ -s \$hmmfile ]
    then
        hmmscan -o /dev/null --noali --cut_ga --domtblout domtblout.txt --cpu ${task.cpus} \$hmmfile $amino
    else
        touch domtblout.txt
    fi
    """
}

process create_solution_space {
    clusterOptions "-g $groupRoot/create_solution_space -R 'rusage[scratch=5120]'"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/$name" }

    input:
    path hmmdb from hmmdb_ch.collect()
    tuple val(acc), path("domtblout.txt") from hmmscan_output_ch

    output:
    tuple val(acc), path("*.hmm*") into dbspace_ch

    script:
    """
    hmmfile=\$(echo *.hmm)
    if [ -s \$hmmfile ]
    then
        $scriptDir/create_solution_space.py domtblout.txt \$hmmfile accspace.txt dbspace.hmm $params.seed
        if [ -s dbspace.hmm ]
        then
            hmmfetch --index dbspace.hmm
        fi
    else
        touch dbspace.hmm
    fi
    """
}

cds_nucl_ch
    .join(dbspace_ch)
    .splitFasta(by:params.chunkSize, file:true, elem:1)
    .set { cds_nucl_db_split_ch }

process iseq_scan {
    clusterOptions "-g $groupRoot/iseq_scan -R 'rusage[scratch=${task.attempt * 5120}]'"
    errorStrategy "retry"
    maxRetries 4
    memory { 6.GB * task.attempt }
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/chunks/$name" }
    scratch true
    stageInMode "copy"

    input:
    tuple val(acc), path(nucl), path(dbspace) from cds_nucl_db_split_ch

    output:
    tuple val(acc), path("output.*.gff") into iseq_output_split_ch
    tuple val(acc), path("oamino.*.fasta") into iseq_oamino_split_ch
    tuple val(acc), path("ocodon.*.fasta") into iseq_ocodon_split_ch

    script:
    chunk = nucl.name.toString().tokenize('.')[-2]
    """
    hmmfile=\$(echo *.hmm)
    if [ -s \$hmmfile ]
    then
        iseq pscan3 \$hmmfile $nucl --hit-prefix chunk_${chunk}_item\
            --output output.${chunk}.gff --oamino oamino.${chunk}.fasta\
            --ocodon ocodon.${chunk}.fasta\
            --no-cut-ga --quiet
    else
        echo "##gff-version 3" > output.${chunk}.gff
        touch oamino.${chunk}.fasta
        touch ocodon.${chunk}.fasta
    fi
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
    clusterOptions "-g $groupRoot/save_output"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/${name}" }, overwrite: true

    input:
    tuple val(name), path(acc) from iseq_results_ch

    output:
    path(name)

    script:
    """
    mv $acc $name
    """
}
