#!/usr/bin/env nextflow

println " "
println "Nextflow pipeline"
println " "
println "using profile: $workflow.profile"
println "using params: $params"
println " "

scriptDir = file("$projectDir/script")
groupRoot = "/horta/$workflow.runName"
profAccsFile = file("$projectDir/profile_accessions.txt")

assembly_ch = Channel.fromPath(params.assemblyFile).collect()

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
    path "db.hmm" into original_hmmfile_ch

    script:
    """
    curl -s $params.hmmfile | gunzip -c > db.hmm
    """
}

process profile_accessions {
    clusterOptions "-g $groupRoot/profile_accessions"
    memory "1 GB"
    publishDir params.outputDir, mode:"copy"

    input:
    path "db.hmm" from original_hmmfile_ch
    path "accs.txt" from profAccsFile

    output:
    path "profile_accessions.txt" into profile_accessions_ch

    script:
    """
    for acc in \$(cat accs.txt)
    do
        $scriptDir/show_profile_accession.py db.hmm \$acc >> profile_accessions.txt
    done
    """
}

process filter_pfam_hmm {
    clusterOptions "-g $groupRoot/filter_pfam_hmm"
    memory "16 GB"
    publishDir params.outputDir, mode:"copy"

    input:
    path "db.hmm.in" from original_hmmfile_ch
    path "profile_accessions.txt" from profile_accessions_ch

    output:
    path "db.hmm" into hmmfile_ch

    script:
    """
    $scriptDir/hmmfetch.py db.hmm.in profile_accessions.txt > db.hmm
    """
}

process press_hmmfile {
    clusterOptions "-g $groupRoot/press_hmmfile"
    publishDir params.outputDir, mode:"copy"

    input:
    path hmmfile from hmmfile_ch

    output:
    path "db.hmm*", includeInputs: true into hmmdb_ch

    script:
    """
    hmmpress $hmmfile
    hmmfetch --index db.hmm
    """
}

process prokka_assembly {
    clusterOptions "-g $groupRoot/prokka_assembly -R 'rusage[scratch=5120]'"
    cpus "${ Math.min(2, params.maxCPUs as int) }"
    memory "8 GB"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "$name" }
    /* scratch true */
    /* stageInMode "copy" */

    input:
    path hmmdb from hmmdb_ch.collect()
    path assembly from assembly_ch

    output:
    path "prokka/*" into assembly_prokka_ch

    script:
    """
    prokka --outdir prokka --prefix assembly --cpus ${task.cpus} $assembly
    """
}

process alignment {
    clusterOptions "-g $groupRoot/alignment"
    memory "6 GB"
    cpus "${ Math.min(12, params.maxCPUs as int) }"
    publishDir params.outputDir, mode:"copy"

    input:
    path assembly from assembly_ch

    output:
    path "alignment.sam" into alignment_sam_ch
    path "alignment.bam" into alignment_bam_ch
    path "alignment.fasta" into alignment_fasta_ch

    script:
    targets = params.targetsFile
    """
    minimap2 -ax map-ont -t ${task.cpus} $assembly $targets | grep --invert-match "^@PG" > unsorted.sam
    samtools sort -O BAM --no-PG -n --threads ${task.cpus} unsorted.sam > all.bam
    samtools view -O SAM --no-PG -h -q 60 all.bam > alignment.sam
    samtools view -O BAM --no-PG alignment.sam > alignment.bam
    samtools fasta alignment.bam > alignment.fasta
    """
}

alignment_fasta_ch
    .splitFasta(by:params.chunkSize, file:true)
    .set { targets_split_ch }

process iseq_scan {
    clusterOptions "-g $groupRoot/iseq_scan -R 'rusage[scratch=${task.attempt * 5120}]'"
    /* errorStrategy "retry" */
    /* maxRetries 4 */
    memory { 6.GB * task.attempt }
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "chunks/$name" }
    /* scratch true */
    /* stageInMode "copy" */

    input:
    path targets from targets_split_ch
    path hmmdb from hmmdb_ch.collect()

    output:
    path("output.*.gff") into iseq_output_split_ch
    path("oamino.*.fasta") into iseq_oamino_split_ch
    path("ocodon.*.fasta") into iseq_ocodon_split_ch

    script:
    chunk = targets.name.toString().tokenize('.')[-2]
    """
    iseq pscan3 db.hmm $targets --hit-prefix chunk_${chunk}_item\
        --output output.${chunk}.gff --oamino oamino.${chunk}.fasta\
        --ocodon ocodon.${chunk}.fasta\
        --no-cut-ga --quiet
    """
}

iseq_output_split_ch
    .collectFile(name: "output.gff", skip:1)
    .set { iseq_output_ch }


iseq_oamino_split_ch
    .collectFile(name: "oamino.fasta")
    .set { iseq_oamino_ch }

iseq_ocodon_split_ch
    .collectFile(name: "ocodon.fasta")
    .set { iseq_ocodon_ch }

iseq_output_ch
    .mix(iseq_oamino_ch, iseq_ocodon_ch)
    .set { iseq_results_ch }

process save_output {
    clusterOptions "-g $groupRoot/save_output"
    publishDir params.outputDir, mode:"copy", overwrite: true

    input:
    path result from iseq_results_ch

    output:
    path result, includeInputs: true

    script:
    """
    """
}

/* process hmmscan { */
/*     clusterOptions "-g $groupRoot/hmmscan -R 'rusage[scratch=5120]'" */
/*     cpus 4 */
/*     memory "8 GB" */
/*     publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/$name" } */
/*     scratch true */
/*     stageInMode "copy" */

/*     input: */
/*     path hmmdb from hmmdb_ch.collect() */
/*     tuple val(acc), path(amino) from cds_amino_ch */

/*     output: */
/*     tuple val(acc), path("domtblout.txt") into hmmscan_output_ch */

/*     script: */
/*     """ */
/*     hmmfile=\$(echo *.hmm) */
/*     if [ -s \$hmmfile ] */
/*     then */
/*         hmmscan -o /dev/null --noali --cut_ga --domtblout domtblout.txt --cpu ${task.cpus} \$hmmfile $amino */
/*     else */
/*         touch domtblout.txt */
/*     fi */
/*     """ */
/* } */

/* process create_solution_space { */
/*     clusterOptions "-g $groupRoot/create_solution_space -R 'rusage[scratch=5120]'" */
/*     publishDir params.outputDir, mode:"copy", saveAs: { name -> "${acc}/$name" } */

/*     input: */
/*     path hmmdb from hmmdb_ch.collect() */
/*     tuple val(acc), path("domtblout.txt") from hmmscan_output_ch */

/*     output: */
/*     tuple val(acc), path("*.hmm*") into dbspace_ch */

/*     script: */
/*     """ */
/*     hmmfile=\$(echo *.hmm) */
/*     if [ -s \$hmmfile ] */
/*     then */
/*         $scriptDir/create_solution_space.py domtblout.txt \$hmmfile accspace.txt dbspace.hmm $params.seed */
/*         if [ -s dbspace.hmm ] */
/*         then */
/*             hmmfetch --index dbspace.hmm */
/*         fi */
/*     else */
/*         touch dbspace.hmm */
/*     fi */
/*     """ */
/* } */

/* cds_nucl_ch */
/*     .join(dbspace_ch) */
/*     .splitFasta(by:params.chunkSize, file:true, elem:1) */
/*     .set { cds_nucl_db_split_ch } */


/* iseq_output_split_ch */
/*     .collectFile(keepHeader:true, skip:1) */
/*     .map { ["output.gff", it] } */
/*     .set { iseq_output_ch } */

/* iseq_oamino_split_ch */
/*     .collectFile() */
/*     .map { ["oamino.fasta", it] } */
/*     .set { iseq_oamino_ch } */

/* iseq_ocodon_split_ch */
/*     .collectFile() */
/*     .map { ["ocodon.fasta", it] } */
/*     .set { iseq_ocodon_ch } */

/* iseq_output_ch */
/*     .mix(iseq_oamino_ch, iseq_ocodon_ch) */
/*     .set { iseq_results_ch } */
