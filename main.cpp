#include "pgenlib_read.h"
#include "pgenlib_write.h"
#include <vector>
#include <cstdio>
#include <algorithm>
#include <memory>
#include <iostream>
using namespace std;

#ifdef __cplusplus
using namespace plink2;
#endif

int32_t main(int32_t argc, char** argv) {
    PglErr reterr = kPglRetSuccess;
    PgenHeaderCtrl header_ctrl;
    PgenFileInfo pgfi;
    PgenReader pgr;
    PreinitPgfi(&pgfi);
    PreinitPgr(&pgr);

    uint32_t sample_ct = 0xffffffffU;
    char errstr_buf[kPglErrstrBufBlen];
    uintptr_t cur_alloc_cacheline_ct;

    // Initialize Phase 1
    reterr = PgfiInitPhase1("plink2.pgen", nullptr, 0xffffffffU, sample_ct, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, errstr_buf);
    if (reterr) {
        fputs(errstr_buf, stderr);
        return 1;
    }

    sample_ct = pgfi.raw_sample_ct;
    if (!sample_ct) {
        fprintf(stderr, "error: sample_ct == 0\n");
        return 1;
    }

    const uint32_t variant_ct = pgfi.raw_variant_ct;
    if (!variant_ct) {
        fprintf(stderr, "error: variant_ct == 0\n");
        return 1;
    }

    printf("%u variant%s and %u sample%s detected.\n", variant_ct, (variant_ct == 1) ? "" : "s", sample_ct, (sample_ct == 1) ? "" : "s");

    // Allocate memory for pgfi using vector
    size_t pgfi_alloc_size = cur_alloc_cacheline_ct * kCacheline;
    std::vector<unsigned char> pgfi_alloc(pgfi_alloc_size);


    // Initialize Phase 2 - corrected parameters
    uint32_t max_vrec_width = 0;// pgfi.const_vrec_width;  // This will be filled by the function

    reterr = PgfiInitPhase2(header_ctrl,
        0,              
        0,              
        0,              
        0,              
        variant_ct,    
        &max_vrec_width,
        &pgfi,
        pgfi_alloc.data(),
        &cur_alloc_cacheline_ct,
        errstr_buf);

    if (reterr) {
        fputs(errstr_buf, stderr);
        return 1;
    }







    // Calculate memory requirements for reader
    const uint32_t raw_sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t genovec_bytes_req = raw_sample_ctl * sizeof(intptr_t);
    cur_alloc_cacheline_ct = DivUp(genovec_bytes_req, kCacheline);

    // Allocate memory for pgr using vector
    size_t pgr_alloc_size = cur_alloc_cacheline_ct * kCacheline;
    std::vector<unsigned char> pgr_alloc(pgr_alloc_size);



    // Initialize reader - corrected parameters
    reterr = PgrInit("plink2.pgen", max_vrec_width, &pgfi, &pgr, pgr_alloc.data());
    if (reterr) {
        fprintf(stderr, "PgrInit failed.\n");
        CleanupPgfi(&pgfi, &reterr);
        return 1;
    }

    // Allocate memory for genotype buffer using vector
    std::vector<uintptr_t> genovec(raw_sample_ctl);

    // Define subset index for all samples (no subsetting)
    PgrSampleSubsetIndex null_subset_index;
    PgrSetSampleSubsetIndex(nullptr, &pgr, &null_subset_index);

    // Define the chunk size (number of variants to read at once)
    const uint32_t chunk_size = 100;  // Adjust based on your memory constraints
    bool error_occurred = false;

    // Process the file in chunks
    for (uint32_t variant_idx = 0; variant_idx < variant_ct && !error_occurred; variant_idx += chunk_size) {
        uint32_t cur_chunk_size = std::min(chunk_size, variant_ct - variant_idx);
        printf("Processing variants %u to %u...\n", variant_idx, variant_idx + cur_chunk_size - 1);

        // Process each variant in the current chunk
        for (uint32_t i = 0; i < cur_chunk_size && !error_occurred; ++i) {
            uint32_t cur_variant_idx = variant_idx + i;

            // Read hardcalls for the current variant - corrected parameters
            reterr = PgrGet(nullptr, null_subset_index, sample_ct, cur_variant_idx, &pgr, genovec.data());
            if (reterr) {
                fprintf(stderr, "Error reading variant %u\n", cur_variant_idx);
                error_occurred = true;
                break;
            }

            // Process genovec data as needed
            // Example: Count allele frequencies
            uint64_t allele_counts[4] = { 0 };
            for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
                uintptr_t geno_word = genovec[widx];

                // Each genotype uses 2 bits: 
                // 00 = homozygous ref, 01 = het, 10 = homozygous alt, 11 = missing
                for (uint32_t bit_pos = 0; bit_pos < kBitsPerWord; bit_pos += 2) {
                    uint32_t genotype = (geno_word >> bit_pos) & 3;
                    allele_counts[genotype]++;
                }
            }

            // Adjust the counts for the last word if sample_ct is not a multiple of 32
            if (sample_ct % 32) {
                uint32_t extra_samples = sample_ct % 32;
                uint32_t last_word_idx = raw_sample_ctl - 1;

                // Mask out bits beyond the last sample
                for (uint32_t bit_pos = extra_samples * 2; bit_pos < kBitsPerWord; bit_pos += 2) {
                    uint32_t genotype = (genovec[last_word_idx] >> bit_pos) & 3;
                    allele_counts[genotype]--;
                }
            }

            printf("Variant %u: Hom Ref: %llu, Het: %llu, Hom Alt: %llu, Missing: %llu\n",
                cur_variant_idx,
                (unsigned long long)allele_counts[0],
                (unsigned long long)allele_counts[1],
                (unsigned long long)allele_counts[2],
                (unsigned long long)allele_counts[3]);
        }
    }

    // Clean up resources
    CleanupPgr(&pgr, &reterr);
    CleanupPgfi(&pgfi, &reterr);

    return (reterr != kPglRetSuccess || error_occurred) ? 1 : 0;
}