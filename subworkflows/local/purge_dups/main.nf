include { PURGEDUPS_SPLITFA } from '../../../modules/nf-core/purgedups/splitfa/main'
include { PURGEDUPS_PBCSTAT } from '../../../modules/nf-core/purgedups/pbcstat/main'
include { PURGEDUPS_CALCUTS } from '../../../modules/nf-core/purgedups/calcuts/main'
include { PURGEDUPS_PURGEDUPS } from '../../../modules/nf-core/purgedups/purgedups/main'
include { PURGEDUPS_GETSEQS } from '../../../modules/nf-core/purgedups/getseqs/main'
include { MINIMAP2_ALIGN as MINIMAP2_READS } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_SELF } from '../../../modules/nf-core/minimap2/align/main'

workflow PURGE_DUPS {
    take:
    assembly      // meta, fasta
    reads         // meta, reads
    
    main:
    Channel.empty().set { ch_versions }
    
    // Step 1: Split assembly
    PURGEDUPS_SPLITFA(assembly)
    ch_versions = ch_versions.mix(PURGEDUPS_SPLITFA.out.versions)
    
    // Step 2: Align reads to assembly for coverage
    MINIMAP2_READS(
        reads,
        assembly.map { meta, fasta -> [fasta] },
        true,  // bam_format
        false, // cigar_paf
        false  // cigar_bam
    )
    ch_versions = ch_versions.mix(MINIMAP2_READS.out.versions)
    
    // Step 3: Self-align split assembly
    MINIMAP2_SELF(
        PURGEDUPS_SPLITFA.out.split_fasta,
        PURGEDUPS_SPLITFA.out.split_fasta.map { meta, fasta -> [fasta] },
        false, // bam_format
        false, // cigar_paf
        false  // cigar_bam
    )
    
    // Step 4: Calculate coverage stats
    PURGEDUPS_PBCSTAT(MINIMAP2_READS.out.paf)
    ch_versions = ch_versions.mix(PURGEDUPS_PBCSTAT.out.versions)
    
    // Step 5: Calculate cutoffs
    PURGEDUPS_CALCUTS(PURGEDUPS_PBCSTAT.out.stat)
    ch_versions = ch_versions.mix(PURGEDUPS_CALCUTS.out.versions)
    
    // Step 6: Identify duplicates
    PURGEDUPS_PBCSTAT.out.basecov
        .join(PURGEDUPS_CALCUTS.out.cutoff)
        .join(MINIMAP2_SELF.out.paf)
        .set { purge_input }
    
    PURGEDUPS_PURGEDUPS(purge_input)
    ch_versions = ch_versions.mix(PURGEDUPS_PURGEDUPS.out.versions)
    
    // Step 7: Get purged sequences
    assembly
        .join(PURGEDUPS_PURGEDUPS.out.bed)
        .set { getseqs_input }
    
    PURGEDUPS_GETSEQS(getseqs_input)
    ch_versions = ch_versions.mix(PURGEDUPS_GETSEQS.out.versions)
    
    emit:
    purged_assembly = PURGEDUPS_GETSEQS.out.purged
    haplotigs       = PURGEDUPS_GETSEQS.out.haplotigs
    versions        = ch_versions
}
