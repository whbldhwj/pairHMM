#include <stdio.h>
#include <ap_int.h>
#include <math.h>
#include <hls_math.h>
#include <string.h>

#define DEBUG

// typedef ap_uint<8> uint8;
typedef uint8_t uint8;
typedef double prob_t;

#define MAX_HAPLOTYPE_LENGTH 100 // bases
#define MAX_READ_LENGTH 100 // bases
#define PADDED_MAX_HAPLOTYPE_LENGTH MAX_HAPLOTYPE_LENGTH+1 // bases
#define PADDED_MAX_READ_LENGTH MAX_READ_LENGTH+1 // bases
#define HAPLOTYPE_NUM 10
#define READ_NUM 10

// -----------------------
// PairHMMmodel.java -> PairHMMModel
// -----------------------
#define TRANS_PROB_ARRAY_LENGTH 6 // length of the standard transition probability array
#define MATCH_TO_MATCH 0 // position in the transition probability array for the Match-to-Match transition
#define INDEL_TO_MATCH 1 // position in the transition probability array for the Indel-to-Match transition
#define MATCH_TO_INSERTION 2
#define INSERTION_TO_INSERTION 3
#define MATCH_TO_DELETION 4
#define DELETION_TO_DELETION 5

// -----------------------
// N2MemoryPairHMM.java
// -----------------------
#define doNotUseTristeCorrection false

// -----------------------
// Log10PairHMM.java
// -----------------------
const prob_t NEGATIVE_INFINITY = -1;

// -----------------------
// PairHMMModel.java
// -----------------------
#define MAX_QUAL 254

void pairHMM_kernel(uint8* haplotype, // byte
        			uint8* read,
					int* haplotypeLength, // bases
					int* readLength, // bases
					int haplotypeNum,
					int readNum,
					uint8* readBaseQualities,
					uint8* readBaseInsertionQualities,
					uint8* readBaseDeletionQualities,
					uint8* readGCP,
					prob_t* mLikelihoodArray);
