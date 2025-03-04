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

int32_t main(int32_t argc, char** argv)
{
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

    // Allocate memory for pgfi
    size_t pgfi_alloc_size = cur_alloc_cacheline_ct * kCacheline;
    std::vector<unsigned char> pgfi_alloc(pgfi_alloc_size);

    // Initialize Phase 2
    uint32_t max_vrec_width = 0; // Will be filled by the function

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

    // Increase memory allocation for reader to prevent buffer overflows
    uintptr_t pgr_alloc_cacheline_ct = DivUp(raw_sample_ctl * sizeof(uintptr_t), kCacheline);
    // CHANGE 1: More conservative allocation for reader buffers (10x instead of 6x)
    pgr_alloc_cacheline_ct += DivUp(10 * max_vrec_width, kCacheline);
    // CHANGE 2: Add more space for working buffers (12x instead of 8x)
    pgr_alloc_cacheline_ct += DivUp(raw_sample_ctl * sizeof(uintptr_t) * 12, kCacheline);

    // Allocate memory for pgr using vector
    size_t pgr_alloc_size = pgr_alloc_cacheline_ct * kCacheline;
    std::vector<unsigned char> pgr_alloc(pgr_alloc_size);

    // Initialize reader - make sure we're using the actual file path
    const char* pgen_filename = "plink2.pgen";

    // Print debug info for memory allocation
    printf("Debug: max_vrec_width=%u, pgr_alloc_size=%zu bytes\n",
        max_vrec_width, pgr_alloc_size);

    // Make sure our allocation is big enough
    if (pgr_alloc.size() < pgr_alloc_size) {
        fprintf(stderr, "Internal error: pgr_alloc buffer too small\n");
        CleanupPgfi(&pgfi, &reterr);
        return 1;
    }

    reterr = PgrInit(pgen_filename, max_vrec_width, &pgfi, &pgr, pgr_alloc.data());
    if (reterr) {
        fprintf(stderr, "PgrInit failed.\n");
        CleanupPgfi(&pgfi, &reterr);
        return 1;
    }

    std::vector<uintptr_t> genovec(sample_ct / 8);

    // Define subset index for all samples (no subsetting)
    PgrSampleSubsetIndex null_subset_index;
    PgrClearSampleSubsetIndex(&pgr, &null_subset_index);  // Properly initialize
    PgrSetSampleSubsetIndex(nullptr, &pgr, &null_subset_index);

    // Define the chunk size (number of variants to read at once)
    const uint32_t chunk_size = 100;  // Adjust based on your memory constraints
    bool error_occurred = false;

    // Process the file in chunks
    for (uint32_t variant_idx = 0; variant_idx < variant_ct && !error_occurred; variant_idx += chunk_size) {
        uint32_t cur_chunk_size = std::min(chunk_size, variant_ct - variant_idx);
        printf("Processing variants %u to %u...\n", variant_idx, variant_idx + cur_chunk_size - 1);

        // Process each variant in the current chunk
        for (uint32_t i = 0; i < cur_chunk_size && !error_occurred; ++i)
        {
            uint32_t cur_variant_idx = variant_idx + i;

            // CHANGE 4: Add safeguard against out-of-bounds variant access
            if (cur_variant_idx >= variant_ct) {
                fprintf(stderr, "Error: attempting to read variant %u but variant_ct is %u\n",
                    cur_variant_idx, variant_ct);
                error_occurred = true;
                break;
            }

            // clear genovec before reading - prevent potential memory issues
            std::fill(genovec.begin(), genovec.end(), 0);

            // Read hardcalls for the current variant
            reterr = PgrGet(nullptr, null_subset_index, sample_ct, cur_variant_idx, &pgr, genovec.data());
            if (reterr)
            {
                fprintf(stderr, "Error reading variant %u (reterr=%d)\n", cur_variant_idx, (int)reterr);
                error_occurred = true;
                break;
            }

            // Sanity check the data after reading
            // In debug mode, check that there are no illegal values in the genovec
            bool has_illegal_value = false;
            for (uint32_t widx = 0; widx < raw_sample_ctl && !has_illegal_value; ++widx) {
                uintptr_t geno_word = genovec[widx];
                // Check if any 2-bit value is > 3
                // This is a bit-hack to detect invalid values (each genotype is 2 bits)
                uintptr_t detect = (geno_word & (geno_word >> 1)) & UINTPTR_MAX / 3;
                if (detect) {
                    has_illegal_value = true;
                }
            }
            if (has_illegal_value) {
                fprintf(stderr, "Warning: Detected illegal genotype value in variant %u\n", cur_variant_idx);
                // Continue anyway, as this is just a warning
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

            // CHANGE 5: Correctly adjust the counts for the last word
            if (sample_ct % (kBitsPerWord / 2)) {
                uint32_t excess_entries = (raw_sample_ctl * (kBitsPerWord / 2)) - sample_ct;
                uint32_t last_word_idx = raw_sample_ctl - 1;

                // Calculate the starting bit position for padding bits
                uint32_t padding_start_bit = (sample_ct % (kBitsPerWord / 2)) * 2;

                // Safely decrement counts for padding bits
                for (uint32_t bit_pos = padding_start_bit; bit_pos < kBitsPerWord; bit_pos += 2) {
                    uint32_t genotype = (genovec[last_word_idx] >> bit_pos) & 3;
                    if (allele_counts[genotype] > 0) {
                        allele_counts[genotype]--;
                    }
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