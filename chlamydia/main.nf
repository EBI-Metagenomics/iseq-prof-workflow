#!/usr/bin/env nextflow

println " "
println "Nextflow pipeline"
println " "
println "using profile: $workflow.profile"
println "using params: $params"
println " "

scriptDir = file("$projectDir/script")
groupRoot = "/horta/$workflow.runName"
/* profAccsFile = file("$projectDir/profile_accessions.txt") */

assembly_ch = Channel.fromPath(params.assemblyFile).collect()

process save_params {
    clusterOptions "-g $groupRoot/save_params"
    publishDir params.outputDir, mode:"copy"

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
    path "Pfam-A.hmm" into pfam_hmmfile_ch

    script:
    """
    curl -s $params.hmmfile | gunzip -c > Pfam-A.hmm
    """
}

/* Pfam-A.hmm  Pfam-A.hmm.h3f  Pfam-A.hmm.h3i  Pfam-A.hmm.h3m  Pfam-A.hmm.h3m.ssi  Pfam-A.hmm.h3p */
process press_pfam_hmmfile {
    clusterOptions "-g $groupRoot/press_pfam_hmmfile"
    storeDir "$params.storageDir/pfam"

    input:
    path hmmfile from pfam_hmmfile_ch

    output:
    path "${hmmfile}.h3f" into pfam_hmm3f_ch
    path "${hmmfile}.h3i" into pfam_hmm3i_ch
    path "${hmmfile}.h3m" into pfam_hmm3m_ch
    path "${hmmfile}.h3p" into pfam_hmm3p_ch
    path "${hmmfile}.h3m.ssi" into pfam_hmm_ssi_ch

    script:
    """
    hmmpress $hmmfile
    hmmfetch --index $hmmfile
    """
}

pfam_hmmfile_ch.combine(pfam_hmm3f_ch)
               .combine(pfam_hmm3i_ch)
               .combine(pfam_hmm3m_ch)
               .combine(pfam_hmm3p_ch)
               .combine(pfam_hmm_ssi_ch)
               .set { pfam_hmmdb_ch }

process pfam_metafile {
    clusterOptions "-g $groupRoot/pfam_metafile"
    storeDir "$params.storageDir/pfam"

    input:
    path hmmfile from pfam_hmmfile_ch

    output:
    path "${hmmfile}.meta.pkl.gz" into pfam_meta_ch

    script:
    """
    #!/usr/bin/env python

    from hmmer_reader import fetch_metadata

    df = fetch_metadata("$hmmfile")
    df.to_pickle("${hmmfile}.meta.pkl.gz")
    """
}

process prokka_assembly {
    clusterOptions "-g $groupRoot/prokka_assembly"
    cpus "${ Math.min(2, params.maxCPUs as int) }"
    memory "8 GB"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "prokka/$name" }
    storeDir "$params.storageDir/prokka"

    input:
    path assembly from assembly_ch

    output:
    path "assembly.err" into assembly_err_ch
    path "assembly.ffn" into assembly_ffn_ch
    path "assembly.fsa" into assembly_fsa_ch
    path "assembly.gff" into assembly_gff_ch
    path "assembly.sqn" into assembly_sqn_ch
    path "assembly.tsv" into assembly_tsv_ch
    path "assembly.faa" into assembly_faa_ch
    path "assembly.fna" into assembly_fna_ch
    path "assembly.gbk" into assembly_gbk_ch
    path "assembly.log" into assembly_log_ch
    path "assembly.tbl" into assembly_tbl_ch
    path "assembly.txt" into assembly_txt_ch

    script:
    """
    prokka --outdir . --force --prefix assembly --cpus ${task.cpus} $assembly
    """
}

process uniprotkb_accessions {
    clusterOptions "-g $groupRoot/uniprotkb_accessions"
    memory "8 GB"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "$name" }

    input:
    path assembly_gff from assembly_gff_ch

    output:
    path "uniprotkb_accessions.txt" into uniprotkb_accessions_ch

    script:
    """
    #!/usr/bin/env python
    import re

    from gff_io import read_gff

    uniprotkb_accs = set()
    pattern = re.compile(".*UniProtKB:([^:,]+)")

    for gff in read_gff("${assembly_gff}"):
        atts = gff.attributes_asdict()

        inference = atts.get("inference", None)
        if inference is None:
            continue

        m = re.match(pattern, inference)
        if m is None:
            continue

        assert len(m.groups()) == 1
        uniprotkb_accs.add(m.groups()[0])

    with open("uniprotkb_accessions.txt", "w") as file:
        for acc in uniprotkb_accs:
            file.write(acc + "\\n")
    """
}

/* process alignment { */
/*     clusterOptions "-g $groupRoot/alignment" */
/*     memory "6 GB" */
/*     cpus "${ Math.min(12, params.maxCPUs as int) }" */
/*     publishDir params.outputDir, mode:"copy" */

/*     input: */
/*     path assembly from assembly_ch */

/*     output: */
/*     path "alignment.sam" into alignment_sam_ch */
/*     path "alignment.bam" into alignment_bam_ch */
/*     path "alignment.fasta" into alignment_fasta_ch */

/*     script: */
/*     targets = params.targetsFile */
/*     """ */
/*     minimap2 -ax map-ont -t ${task.cpus} $assembly $targets | grep --invert-match "^@PG" > unsorted.sam */
/*     samtools sort -O BAM --no-PG -n --threads ${task.cpus} unsorted.sam > all.bam */
/*     samtools view -O SAM --no-PG -h -q 60 all.bam > alignment.sam */
/*     samtools view -O BAM --no-PG alignment.sam > alignment.bam */
/*     samtools fasta alignment.bam > alignment.fasta */
/*     """ */
/* } */

/* alignment_fasta_ch */
/*     .splitFasta(by:params.chunkSize, file:true) */
/*     .set { targets_split_ch } */

/* process iseq_scan { */
/*     clusterOptions "-g $groupRoot/iseq_scan -R 'rusage[scratch=${task.attempt * 5120}]'" */
/*     /1* errorStrategy "retry" *1/ */
/*     /1* maxRetries 4 *1/ */
/*     memory { 6.GB * task.attempt } */
/*     publishDir params.outputDir, mode:"copy", saveAs: { name -> "chunks/$name" } */
/*     /1* scratch true *1/ */
/*     /1* stageInMode "copy" *1/ */

/*     input: */
/*     path targets from targets_split_ch */
/*     path hmmdb from hmmdb_ch.collect() */

/*     output: */
/*     path("output.*.gff") into iseq_output_split_ch */
/*     path("oamino.*.fasta") into iseq_oamino_split_ch */
/*     path("ocodon.*.fasta") into iseq_ocodon_split_ch */

/*     script: */
/*     chunk = targets.name.toString().tokenize('.')[-2] */
/*     """ */
/*     iseq pscan3 db.hmm $targets --hit-prefix chunk_${chunk}_item\ */
/*         --output output.${chunk}.gff --oamino oamino.${chunk}.fasta\ */
/*         --ocodon ocodon.${chunk}.fasta\ */
/*         --no-cut-ga --quiet */
/*     """ */
/* } */

/* iseq_output_split_ch */
/*     .collectFile(name: "output.gff", skip:1) */
/*     .set { iseq_output_ch } */


/* iseq_oamino_split_ch */
/*     .collectFile(name: "oamino.fasta") */
/*     .set { iseq_oamino_ch } */

/* iseq_ocodon_split_ch */
/*     .collectFile(name: "ocodon.fasta") */
/*     .set { iseq_ocodon_ch } */

/* iseq_output_ch */
/*     .mix(iseq_oamino_ch, iseq_ocodon_ch) */
/*     .set { iseq_results_ch } */

/* process save_output { */
/*     clusterOptions "-g $groupRoot/save_output" */
/*     publishDir params.outputDir, mode:"copy", overwrite: true */

/*     input: */
/*     path result from iseq_results_ch */

/*     output: */
/*     path result, includeInputs: true */

/*     script: */
/*     """ */
/*     """ */
/* } */

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
