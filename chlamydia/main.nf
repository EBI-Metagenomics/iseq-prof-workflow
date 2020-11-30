#!/usr/bin/env nextflow

scriptDir = file("$projectDir/script")
groupRoot = "/horta/$workflow.runName"

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
    path "${hmmfile}.h3m.ssi" into pfam_hmmssi_ch

    script:
    """
    hmmpress $hmmfile
    hmmfetch --index $hmmfile
    """
}

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
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "prokka_assembly/$name" }

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

process hmmscan_assembly {
    clusterOptions "-g $groupRoot/hmmscan_assembly -R 'rusage[scratch=5120]'"
    cpus "${ Math.min(4, params.maxCPUs as int) }"
    memory "8 GB"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "hmmscan_assembly/$name" }
    /* scratch true */
    /* stageInMode "copy" */

    input:
    path pfam_hmmfile from pfam_hmmfile_ch
    path pfam_hmm3f from pfam_hmm3f_ch
    path pfam_hmm3i from pfam_hmm3i_ch
    path pfam_hmm3m from pfam_hmm3m_ch
    path pfam_hmm3p from pfam_hmm3p_ch
    path pfam_hmmssi from pfam_hmmssi_ch
    path assembly_faa from assembly_faa_ch

    output:
    path "domtblout.txt" into assembly_domtblout_ch

    script:
    """
    hmmscan -o /dev/null --noali --cut_ga --domtblout domtblout.txt \
        --cpu ${task.cpus} $pfam_hmmfile $assembly_faa
    """
}

process iseq_scan_assembly {
    clusterOptions "-g $groupRoot/iseq_scan_assembly -R 'rusage[scratch=5120]'"
    memory "8 GB"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "iseq_scan_assembly/$name" }
    /* scratch true */
    /* stageInMode "copy" */

    when:
    params.scanAssembly == "yes"

    input:
    path pfam_hmmfile from pfam_hmmfile_ch
    path pfam_hmm3f from pfam_hmm3f_ch
    path pfam_hmm3i from pfam_hmm3i_ch
    path pfam_hmm3m from pfam_hmm3m_ch
    path pfam_hmm3p from pfam_hmm3p_ch
    path pfam_hmmssi from pfam_hmmssi_ch
    path assembly from assembly_ch

    output:
    path "output.gff" into assembly_iseq_output_ch
    path "oamino.fasta" into assembly_iseq_oamino_ch
    path "ocodon.fasta" into assembly_iseq_ocodon_ch

    script:
    """
    iseq pscan3 $pfam_hmmfile $assembly --hit-prefix item\
        --output output.gff --oamino oamino.fasta\
        --ocodon ocodon.fasta\
        --no-cut-ga --quiet --epsilon $params.iseqEpsilon
    """
}

process create_hmmdb_solution_space {
    clusterOptions "-g $groupRoot/create_hmmdb_solution_space -R 'rusage[scratch=5120]'"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "$name" }

    input:
    path pfam_hmmfile from pfam_hmmfile_ch
    path pfam_hmm3f from pfam_hmm3f_ch
    path pfam_hmm3i from pfam_hmm3i_ch
    path pfam_hmm3m from pfam_hmm3m_ch
    path pfam_hmm3p from pfam_hmm3p_ch
    path pfam_hmmssi from pfam_hmmssi_ch
    path assembly_domtblout from assembly_domtblout_ch
    path pfam_meta from pfam_meta_ch

    output:
    path "db.hmm" into db_hmmfile_ch

    script:
    """
    #!/usr/bin/env python

    from pathlib import Path
    from hmmer import HMMER, read_domtbl
    import hmmer_reader
    import pandas as pd
    from numpy.random import RandomState

    random = RandomState($params.seed)
    meta_filepath = Path("$pfam_meta")
    dombtbl_filepath = Path("$assembly_domtblout")

    meta = pd.read_pickle(meta_filepath)
    rows = read_domtbl(dombtbl_filepath)

    true_profiles = [row.target.accession for row in rows]
    all_false_profiles = set(meta["ACC"].tolist()) - set(true_profiles)
    false_profiles = list(random.choice(list(all_false_profiles), size=100, replace=False))

    hmmer = HMMER("$pfam_hmmfile")
    with open("db.hmm", "w") as file:
        file.write(hmmer.fetch(true_profiles + false_profiles))
    """
}

process press_db_hmmfile {
    clusterOptions "-g $groupRoot/press_db_hmmfile"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "$name" }

    input:
    path hmmfile from db_hmmfile_ch

    output:
    path "${hmmfile}.h3f" into db_hmm3f_ch
    path "${hmmfile}.h3i" into db_hmm3i_ch
    path "${hmmfile}.h3m" into db_hmm3m_ch
    path "${hmmfile}.h3p" into db_hmm3p_ch
    path "${hmmfile}.h3m.ssi" into db_hmmssi_ch

    script:
    """
    hmmpress $hmmfile
    hmmfetch --index $hmmfile
    """
}

process alignment {
    clusterOptions "-g $groupRoot/alignment"
    memory "6 GB"
    cpus "${ Math.min(12, params.maxCPUs as int) }"
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "alignment/$name" }

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

alignment_fasta_ch.into{ targets_ch1; targets_ch2 }

targets_ch1
    .splitFasta(by:params.chunkSize, file:true)
    .set { targets_chunk_ch }

process iseq_scan_targets {
    clusterOptions "-g $groupRoot/iseq_scan_targets -R 'rusage[scratch=${task.attempt * 5120}]'"
    /* errorStrategy "retry" */
    /* maxRetries 4 */
    memory { 6.GB * task.attempt }
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "iseq_scan_targets/chunks/$name" }
    /* scratch true */
    /* stageInMode "copy" */

    input:
    path targets from targets_chunk_ch
    path hmmfile from db_hmmfile_ch.collect()
    path "${hmmfile}.h3f" from db_hmm3f_ch.collect()
    path "${hmmfile}.h3i" from db_hmm3i_ch.collect()
    path "${hmmfile}.h3m" from db_hmm3m_ch.collect()
    path "${hmmfile}.h3p" from db_hmm3p_ch.collect()
    path "${hmmfile}.h3m.ssi" from db_hmmssi_ch.collect()

    output:
    path("output.*.gff") into iseq_output_split_ch
    path("oamino.*.fasta") into iseq_oamino_split_ch
    path("ocodon.*.fasta") into iseq_ocodon_split_ch

    script:
    chunk = targets.name.toString().tokenize('.')[-2]
    """
    iseq pscan3 $hmmfile $targets --hit-prefix chunk_${chunk}_item\
        --output output.${chunk}.gff --oamino oamino.${chunk}.fasta\
        --ocodon ocodon.${chunk}.fasta\
        --no-cut-ga --quiet
    """
}

process prokka_targets {
    clusterOptions "-g $groupRoot/prokka_targets"
    memory { 4.GB }
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "prokka_targets/$name" }

    input:
    path targets from targets_ch2

    output:
    path "targets.err" into targets_err_ch
    path "targets.ffn" into targets_ffn_ch
    path "targets.fsa" into targets_fsa_ch
    path "targets.gff" into targets_gff_ch
    path "targets.sqn" into targets_sqn_ch
    path "targets.tsv" into targets_tsv_ch
    path "targets.faa" into targets_faa_ch
    path "targets.fna" into targets_fna_ch
    path "targets.gbk" into targets_gbk_ch
    path "targets.log" into targets_log_ch
    path "targets.tbl" into targets_tbl_ch
    path "targets.txt" into targets_txt_ch

    script:
    """
    prokka --outdir . --force --prefix targets --cpus ${task.cpus} $targets
    """
}

targets_faa_ch.into{ amino_targets_ch }

process hmmscan_targets {
    clusterOptions "-g $groupRoot/hmmscan_targets -R 'rusage[scratch=${task.attempt * 5120}]'"
    /* errorStrategy "retry" */
    /* maxRetries 4 */
    memory { 6.GB * task.attempt }
    publishDir params.outputDir, mode:"copy", saveAs: { name -> "hmmscan_targets/$name" }
    /* scratch true */
    /* stageInMode "copy" */

    input:
    path targets from amino_targets_ch
    path hmmfile from db_hmmfile_ch.collect()
    path "${hmmfile}.h3f" from db_hmm3f_ch.collect()
    path "${hmmfile}.h3i" from db_hmm3i_ch.collect()
    path "${hmmfile}.h3m" from db_hmm3m_ch.collect()
    path "${hmmfile}.h3p" from db_hmm3p_ch.collect()
    path "${hmmfile}.h3m.ssi" from db_hmmssi_ch.collect()

    output:
    path("domtblout.txt") into hmmscan_domtblout_ch

    script:
    """
    hmmscan -o /dev/null --noali --cut_ga --domtblout domtblout.txt --cpu ${task.cpus} \
        $hmmfile $targets
    """
}

iseq_output_split_ch
    .collectFile(name: "output.gff", keepHeader:true, skip:1)
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

process merge_iseq_chunks {
    clusterOptions "-g $groupRoot/merge_iseq_chunks"
    publishDir params.outputDir, mode:"copy", overwrite: true, saveAs: { name -> "iseq_scan_targets/$name" }

    input:
    path result from iseq_results_ch

    output:
    path result, includeInputs: true

    script:
    """
    """
}
