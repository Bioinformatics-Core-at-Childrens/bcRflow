#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import processes from src:
include {
    mixcr_align_bulk; mixcr_assemblePartial_1; mixcr_assemblePartial_2;
    mixcr_assemble; mixcr_assembleContigs; mixcr_exportClones; 
    mixcr_align_10X_5p; mixcr_refineTags; downstream
} from './src/processes.nf'

// Print pipeline parameters to stdout:
log.info """\
bcRflow: RNASeq - BCR Pipeline
                                          
===================================

Welcome to bcRflow! This pipeline extracts, aligns, and analyzes BCR repertoires from non-targeted RNA-sequencing data.
The bcRflow pipeline was developed by Brent T. Schlegel and the team at the Bioinformatics Core @ the UPMC Children's Hospital of Pittsburgh.

The use of MiXCR and IMGT in this pipeline are for strictly non-commerical, academic purposes.
In order to obtain a non-commercial license to use MiXCR, please see: https://licensing.milaboratories.com/

For inquiries or technical assistance, please open a new Issue on the bcRflow GitHub page (https://github.com/Bioinformatics-Core-at-Childrens/bcRflow/).
Any reproductions or publications utilizing this pipeline should attribute all proper credit for development to the Developers (see above).

===================================

Project parameters:

- Project Name          : ${params.project_name}
- Sample List           : ${params.input}
- IMGT Library          : ${params.imgt}
- Species               : ${params.species}
- Chain                 : ${params.chain}
- Input Directory       : ${params.indir}
- Output Directory      : ${params.outdir}
- Single Cell?          : ${params.is_sc}

"""

// Bulk RNAseq Workflow:
workflow mixcr_bulk {

    take: 
        samples_channel
        downstream_channel
    main:
        mixcr_align_bulk(samples_channel) | 
        mixcr_assemblePartial_1 | 
        mixcr_assemblePartial_2 | 
        mixcr_assemble | 
        mixcr_assembleContigs | 
        mixcr_exportClones 
        
        downstream(mixcr_exportClones.out.ready.collect().unique(), downstream_channel)
}
// 10X 5' GEX Workflow:
workflow mixcr_10X_GEX {

    take: 
        samples_channel
        downstream_channel
    main:
        mixcr_align_10X_5p(samples_channel) | 
        mixcr_assemblePartial_1 | 
        mixcr_assemblePartial_2 |
        mixcr_assemble |
        mixcr_assembleContigs |
        mixcr_exportClones
        downstream(mixcr_exportClones.out.ready.collect().unique(), downstream_channel)
}

workflow {

    // Input validation
    def valid_species = ['hsa','mmu'] // pipeline only supports Human and Mouse samples
    is_valid_specie = params.species in valid_species
    if (!is_valid_specie) {
        log.error "`params.species` must be one of ${valid_species}"
        exit 1,   "`params.species` must be one of ${valid_species}"
    }

    Channel.fromPath("${params.input}") // check that the master CSV file exists, pass to Input channel
        .ifEmpty { exit 1, "File not foud: ${params.input}" }
        .set { sampleInfoChannel }

   /*
   * Create a channel that emits tuples containing three elements:
   * the SampleID, the first read-pair file and the second read-pair file
   */
    samples_channel = sampleInfoChannel
        .splitCsv(header:true)                  // Read in the CSV
        .map { row -> tuple (row.SampleID,      // Map to tuple, Sample ID
        file(row.R1.trim()),                    // R1 FASTQ
        file(row.R2.trim()),                    // R2 FASTQ
        file(params.imgt))                      // IGMT database file (MiXCR allele reference)
    }

    /*
    * Create a channel that emits a tuple containing four elements:
    * the location of the MiXCR output, the sample info file,
    * the downstream downstream script, and a script containing utility functions 
    * for processing the MiXCR results
    */
    Channel.fromPath("${params.outdir}/MiXCR")                      // Input directory
      .combine(Channel.fromPath("${params.input}"))                 // Reads file (.csv)
      .combine(Channel.fromPath("${baseDir}/src/downstream.R"))     // Path to downstream.R
      .combine(Channel.fromPath("${baseDir}/src/utils.R"))          // Path to utils.R
      .map { outdir, input, downstreamScript, utilsScript ->        // Map to tuple
          tuple(outdir, input, downstreamScript, utilsScript)
      }                                                         
    .set { downstream_channel }
    
    // If the sample is 10X 5' GEX, use single-cell MiXCR preset, else use the bulk preset:
    params.is_sc ? mixcr_10X_GEX(samples_channel, downstream_channel) : mixcr_bulk(samples_channel, downstream_channel)
}

workflow.onComplete {
  log.info "bcRflow completed at: ${workflow.complete}"
  log.info "Time taken: ${workflow.duration}"
  log.info "Execution status: ${workflow.success ? 'success' : 'failed'}"
}

