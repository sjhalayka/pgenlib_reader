#include "pgenlib_read.h"
#include "pgenlib_write.h"

#ifdef __cplusplus
using namespace plink2;
#endif

int32_t main(int32_t argc, char** argv) 
{
	PglErr reterr = kPglRetSuccess;
	unsigned char* pgfi_alloc = nullptr;
	unsigned char* pgr_alloc = nullptr;
	PgenHeaderCtrl header_ctrl;
	PgenFileInfo pgfi;
	PgenReader pgr;
	PreinitPgfi(&pgfi);
	PreinitPgr(&pgr);


	uint32_t sample_ct = 0xffffffffU;

	char errstr_buf[kPglErrstrBufBlen];
	uintptr_t cur_alloc_cacheline_ct;
	reterr = PgfiInitPhase1("data2.pgen", nullptr, 0xffffffffU, sample_ct, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, errstr_buf);
	if (reterr) {
		fputs(errstr_buf, stderr);
		return 0;
	}
	sample_ct = pgfi.raw_sample_ct;
	if (!sample_ct) {
		fprintf(stderr, "error: sample_ct == 0\n");
		return 0;
	}
	const uint32_t variant_ct = pgfi.raw_variant_ct;
	if (!variant_ct) {
		fprintf(stderr, "error: variant_ct == 0\n");
		return 0;
	}

	printf("%u variant%s and %u sample%s detected.\n", variant_ct, (variant_ct == 1) ? "" : "s", sample_ct, (sample_ct == 1) ? "" : "s");

	return 0;

	CleanupPgr(&pgr, &reterr);
	CleanupPgfi(&pgfi, &reterr);

	if (pgfi_alloc)
		aligned_free(pgfi_alloc);

	if (pgr_alloc)
		aligned_free(pgr_alloc);

	return 1;
}
