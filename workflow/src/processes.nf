// Bulk RNA-seq alignment, using IMGT reference
process mixcr_align_bulk{
    cpus "$params.cpus"
    memory { 16.GB * task.attempt }
    time { 12.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "$SampleID"
    label "mixcr" 

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(R1), file(R2), file(imgt)

    output:
    tuple val(SampleID),file(imgt), path("${SampleID}.vdjca"), emit: vdjca

    script:
    """
    mixcr align \
    --library imgt.202214-2.sv8 \
    --preset rna-seq \
    --species $params.species \
    -OvParameters.geneFeatureToAlign="VTranscriptWithP" \
    -OvParameters.parameters.floatingLeftBound=false \
    -OjParameters.parameters.floatingRightBound=false \
    -OcParameters.parameters.floatingRightBound=false \
    -OallowPartialAlignments=true \
    -OallowNoCDR3PartAlignments=true \
    --threads $params.threads \
    --force-overwrite \
    --report ${SampleID}.report.txt \
    --json-report ${SampleID}.report.json \
    ${R1} ${R2} ${SampleID}.vdjca
    """
}

// 10X 5' GEX scRNA alignment, using IMGT reference
process mixcr_align_10X_5p{
    cpus "$params.cpus"
    memory { 16.GB * task.attempt }
    time { 12.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "$SampleID"
    label "mixcr" 

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(R1), file(R2), file(imgt)

    output:
    tuple val(SampleID),file(imgt), path("${SampleID}.vdjca"), emit: vdjca

    script:
    """
    mixcr align \
    --library imgt.202214-2.sv8 \
    --preset 10x-sc-5gex \
    --species $params.species \
    -OvParameters.geneFeatureToAlign="VTranscriptWithout5UTRWithP" \
    -OvParameters.parameters.floatingLeftBound=false \
    -OjParameters.parameters.floatingRightBound=false \
    -OcParameters.parameters.floatingRightBound=false \
    -OallowPartialAlignments=true \
    -OallowNoCDR3PartAlignments=true \
    --threads $params.threads \
    --force-overwrite \
    --report ${SampleID}.report.txt \
    --json-report ${SampleID}.report.json \
    ${R1} ${R2} ${SampleID}.vdjca
    """
}

// TODO: implement 10X 3' GEX, not necessary to implement at this stage:
process mixcr_10X_3p_GEX{
    cpus "$params.cpus"
    memory { 16.GB * task.attempt }
    time { 12.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true
    input:
    tuple val(SampleID), file(R1), file(R2), file(imgt)

    output:
    path("*")

    script:
    """
    mixcr analyze 10x-sc-5gex \
    --species $params.species \
    --rna \
    --tag-pattern "${'^(CELL:N{16})(UMI:N{12})\\^(R2:*)'}" \
    --library imgt.202214-2.sv8 \
    --impute-germline-on-export \
    -MrefineTagsAndSort.whitelists.CELL=builtin:3M-febrary-2018 \
    --assemble-clonotypes-by 'CDR3' \
    --assemble-contigs-by VDJRegion \
    --append-export-clones-field -cloneId \
    --append-export-clones-field -readCount \
    --append-export-clones-field -aaFeatureImputed VDJRegion \
    --append-export-clones-field -nFeatureImputed VDJRegion \
    --append-export-clones-field -nFeatureImputed FR1 \
    --append-export-clones-field -aaFeatureImputed FR1 \
    --append-export-clones-field -allNFeaturesImputed CDR1Begin FR4End \
    --append-export-clones-field -allAAFeaturesImputed CDR1Begin FR4End \
    --append-export-clones-field -nMutations VRegion \
    --append-export-clones-field -nMutations DRegion \
    --append-export-clones-field -nMutations JRegion \
    --export-productive-clones-only \
    ${R1} \
    ${R2} \
    ${SampleID}
    """
}

// TODO: Single-cell barcode refining only works when barcodes match the whitelist structure.
process mixcr_refineTags{
    cpus "$params.cpus"
    memory "16.GB"
    time "12.hour"

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(imgt), path(vdjca)

    output:
    tuple val(SampleID), file(imgt), path("${SampleID}.corrected.vdjca"), emit: vdjca_corrected

    script:
    """
    #Barcode correction for single-cell 5' GEX data
    #Presets taken from https://github.com/milaboratory/mixcr/blob/develop/src/main/resources/presets/protocols/10x.yaml

    mixcr refineTagsAndSort \
    -p 0.001 \
    -s 0.001 \
    -i 1.0e-05 \
    -q 12 \
    --max-substitutions 2 \
    --max-indels 2 \
    --max-errors 3 \
    --force-overwrite \
    ${vdjca} \
    ${SampleID}.corrected.vdjca
    """
}

// Partial assembly for short reads, round 1:
process mixcr_assemblePartial_1{
    cpus "$params.cpus"
    memory { 16.GB * task.attempt }
    time { 12.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(imgt), path(vdjca)

    output:
    tuple val(SampleID), file(imgt), path("${SampleID}.passembled.1.vdjca"), emit: vdjca_1

    script:
    """
    #assemble overlapping fragmented sequencing reads:
    ## First round:
    mixcr assemblePartial \
      --report ${SampleID}.report.txt \
      --json-report ${SampleID}.report.json \
      --force-overwrite \
      ${vdjca} \
      ${SampleID}.passembled.1.vdjca
    """
}

// Partial assembly for short reads, round 2:
process mixcr_assemblePartial_2{
    cpus "$params.cpus"
    memory { 16.GB * task.attempt }
    time { 12.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "$SampleID"
    label "mixcr" 

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(imgt), path(vdjca_1)

    output:
    tuple val(SampleID), file(imgt), path("${SampleID}.passembled.2.vdjca"), emit: vdjca_2

    script:
    """
    #assemble overlapping fragmented sequencing reads:
    ## Second round
    mixcr assemblePartial \
      --report ${SampleID}.report.txt \
      --json-report ${SampleID}.report.json \
      --force-overwrite \
      ${vdjca_1} \
      ${SampleID}.passembled.2.vdjca
    """
}

// Assemble partial assemblies into final assembly:
process mixcr_assemble{
    cpus "$params.cpus"
    memory { 16.GB * task.attempt }
    time { 12.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(imgt), path(vdjca_2)

    output:
    tuple val(SampleID), file(imgt), path("${SampleID}.clna"), emit: clna

    script:
    """
    #Full assemble:
    mixcr assemble \
        -OassemblingFeatures="CDR3" \
        -OseparateByV=true \
        -OseparateByJ=true \
        --write-alignments \
        --force-overwrite \
        --report ${SampleID}.assembled.report.txt \
        --json-report ${SampleID}.assembled.json \
        ${vdjca_2} \
        ${SampleID}.clna
    """
}

// Assemble contigs:
process mixcr_assembleContigs {
    cpus "$params.cpus"
    memory { 16.GB * task.attempt }
    time { 12.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "$SampleID"
    label "mixcr"
    errorStrategy 'ignore' // common bug with assembleContigs for single-cell data, for now ignore samples that fail

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID),file(imgt), path(clna)

    output:
    tuple val(SampleID),file(imgt), path("${SampleID}.clns"), emit: clns

    script:
    """
    #Assemble contigs based on VDJ region:
    mixcr assembleContigs \
        -OsubCloningRegions=VDJRegion \
        --report ${SampleID}.report.txt \
        --json-report ${SampleID}.report.json \
        --force-overwrite \
        ${clna} \
        ${SampleID}.clns
    """
}

// Export clones to standardized format:
process mixcr_exportClones {
    cpus "$params.cpus"
    memory { 4.GB * task.attempt }
    time { 3.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "$SampleID"
    label "mixcr" 

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true
    input:
    tuple val(SampleID),file(imgt), path(clns)

    output:
    path("*"), emit: all_clones
    val true, emit: ready

    script:
    """
    #Export clones:
    mixcr exportClones \
            --chains IGH \
            --impute-germline-on-export \
            -allAAFeaturesImputed CDR1Begin FR4End \
            -allNFeaturesImputed CDR1Begin FR4End \
            -nFeatureImputed FR1 \
            -aaFeatureImputed FR1 \
            --force-overwrite \
            -nMutations VRegion \
            -nMutations JRegion \
            -nMutations DRegion \
        ${clns} \
        ${SampleID}.tsv

    #Export MiXCR Reports:
    #json-formatted:
    mixcr exportReports \
    --json \
    ${clns} \
    ${SampleID}_MiXCR.report.json

    #txt-formatted:
    mixcr exportReports \
    ${clns} \
    ${SampleID}_MiXCR.report.txt
    """
}

// Downstream processing and figures using downstream:
process downstream {
    cpus "$params.cpus"
    memory { 16.GB * task.attempt }
    time { 12.hour }
    errorStrategy { task.exitStatus in [1, 137, 138, 139, 140, 143] ? 'retry' : 'finish' } // broader error-handling due to nonspecific R error codes
    maxRetries 3

    label 'downstream'
    tag "downstream"

    publishDir "$params.outdir",
        pattern: '*',
        mode: 'copy',
        overwrite: true
    input:
    val ready
    tuple path(inputDir), path(sampleInfo), file(downstream), file(utils)

    output:
    path("*")

    script:
    """
    Rscript ${downstream} -s ${params.species} -m ${sampleInfo} -t ${params.threshold} -d ${params.downsample}
    """
}
