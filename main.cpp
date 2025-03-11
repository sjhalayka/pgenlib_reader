#include "pgenlib_read.h"
#include "pgenlib_write.h"

uint32_t min_uint32_t(uint32_t a, uint32_t b) {
    return (a < b) ? a : b;
}

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

    reterr = PgfiInitPhase1("plink2.pgen", nullptr, 0xffffffffU, sample_ct, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, errstr_buf);
    if (reterr) {
        fputs(errstr_buf, stderr);
        return 1;
    }

    sample_ct = pgfi.raw_sample_ct;
    uint32_t variant_ct = pgfi.raw_variant_ct;
    if (!sample_ct || !variant_ct) {
        fprintf(stderr, "error: sample_ct == 0 or variant_ct == 0\n");
        return 1;
    }

    printf("%u variants and %u samples detected.\n", variant_ct, sample_ct);

    size_t pgfi_alloc_size = cur_alloc_cacheline_ct * kCacheline;
    unsigned char* pgfi_alloc = (unsigned char*)malloc(pgfi_alloc_size);

    uint32_t max_vrec_width = 0;
    reterr = PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, &cur_alloc_cacheline_ct, errstr_buf);
    if (reterr) {
        fputs(errstr_buf, stderr);
        return 1;
    }

    const uint32_t raw_sample_ctl = BitCtToWordCt(sample_ct);
    uintptr_t pgr_alloc_cacheline_ct = DivUp(raw_sample_ctl * sizeof(uintptr_t), kCacheline) + DivUp(10 * max_vrec_width, kCacheline) + DivUp(raw_sample_ctl * sizeof(uintptr_t) * 12, kCacheline);
    size_t pgr_alloc_size = pgr_alloc_cacheline_ct * kCacheline;
    unsigned char* pgr_alloc = (unsigned char*)malloc(pgr_alloc_size);

    const char* pgen_filename = "plink2.pgen";
    printf("Debug: max_vrec_width=%u, pgr_alloc_size=%zu bytes\n", max_vrec_width, pgr_alloc_size);

    reterr = PgrInit(pgen_filename, max_vrec_width, &pgfi, &pgr, pgr_alloc);
    if (reterr) {
        fprintf(stderr, "PgrInit failed.\n");
        CleanupPgfi(&pgfi, &reterr);
        return 1;
    }

    const uint32_t variant_chunk_size = 100;
    const uint32_t sample_chunk_size = 500;  // Define sample chunk size

    int error_occurred = 0;
    uintptr_t* genovec = (uintptr_t*)malloc(sizeof(uintptr_t) * sample_chunk_size / 8);

    for (uint32_t variant_idx = 0; variant_idx < variant_ct && !error_occurred; variant_idx += variant_chunk_size) {
        uint32_t cur_variant_chunk = min_uint32_t(variant_chunk_size, variant_ct - variant_idx);
        printf("Processing variants %u to %u...\n", variant_idx, variant_idx + cur_variant_chunk - 1);


        PgrSampleSubsetIndex sample_subset_index;
        PgrClearSampleSubsetIndex(&pgr, &sample_subset_index);
        PgrSetSampleSubsetIndex(nullptr, &pgr, &sample_subset_index);


        for (uint32_t sample_idx = 0; sample_idx < sample_ct && !error_occurred; sample_idx += sample_chunk_size) 
        {
            uint32_t cur_sample_chunk = min_uint32_t(sample_chunk_size, sample_ct - sample_idx);

            for (uint32_t i = 0; i < cur_variant_chunk && !error_occurred; ++i) 
            {
                uint32_t cur_variant_idx = variant_idx + i;
                reterr = PgrGet(nullptr, sample_subset_index, cur_sample_chunk, cur_variant_idx, &pgr, genovec);
                if (reterr) {
                    fprintf(stderr, "Error reading variant %u for samples %u-%u\n", cur_variant_idx, sample_idx, sample_idx + cur_sample_chunk - 1);
                    error_occurred = 1;
                    break;
                }

                uint64_t allele_counts[4] = { 0 };
                for (uint32_t widx = 0; widx < BitCtToWordCt(cur_sample_chunk); ++widx) {
                    uintptr_t geno_word = genovec[widx];
                    for (uint32_t bit_pos = 0; bit_pos < kBitsPerWord; bit_pos += 2) {
                        uint32_t genotype = (geno_word >> bit_pos) & 3;
                        allele_counts[genotype]++;
                    }
                }

                printf("Variant %u (Samples %u-%u): Hom Ref: %llu, Het: %llu, Hom Alt: %llu, Missing: %llu\n",
                    cur_variant_idx, sample_idx, sample_idx + cur_sample_chunk - 1,
                    (unsigned long long)allele_counts[0],
                    (unsigned long long)allele_counts[1],
                    (unsigned long long)allele_counts[2],
                    (unsigned long long)allele_counts[3]);
            }
        }
    }

    free(pgfi_alloc);
    free(pgr_alloc);
    free(genovec);
    CleanupPgr(&pgr, &reterr);
    CleanupPgfi(&pgfi, &reterr);
    return (reterr != kPglRetSuccess || error_occurred) ? 1 : 0;
}