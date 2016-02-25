#include <stdio.h>
#include <ap_int.h>
#include <math.h>

typedef ap_uint<8> uint8_t;
typedef double prob_t;

#define MAX_HAPLOTYPE_LENGTH 250 // bases
#define MAX_READ_LENGTH 250 // bases
#define PADDED_MAX_HAPLOTYPE_LENGTH MAX_HAPLOTYPE_LENGTH+1 // bases
#define PADDED_MAX_READ_LENGTH MAX_READ_LENGTH+1 // bases
#define HAPLOTYPE_NUM 1000
#define READ_NUM 1000

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
const prob_t NEGATIVE_INFINITY = -1

// -----------------------
// PairHMMModel.java
// -----------------------
#define MAX_QUAL = 254

void pairHMM_kernel(uint8_t haplotype[HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH], // byte
                    uint8_t read[READ_NUM][MAX_READ_LENGTH],
                    int haplotypeLength[HAPLOTYPE_NUM], // bases
                    int readLength[READ_NUM], // bases
                    uint8_t readBaseQualities[READ_NUM][MAX_READ_LENGTH],
                    uint8_t readBaseInsertionQualities[READ_NUM][MAX_READ_LENGTH],
                    uint8_t readDeletionQualities[READ_NUM][MAX_READ_LENGTH],
                    uint8_t readGCP[READ_NUM][MAX_READ_LENGTH]);