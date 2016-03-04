#include "pairHMM.h"

void computeLikelihoods(// uint8 haplotype[HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH],
                        // uint8 read[READ_NUM][MAX_READ_LENGTH],
                        // int haplotypeLength[HAPLOTYPE_NUM], // bases
                        // int readLength[READ_NUM], // bases
                        // uint8 readBaseQualities[READ_NUM][MAX_READ_LENGTH],
                        // uint8 readBaseInsertionQualities[READ_NUM][MAX_READ_LENGTH],
                        // uint8 readBaseDeletionQualities[READ_NUM][MAX_READ_LENGTH],
                        // uint8 readGCP[READ_NUM][MAX_READ_LENGTH],
                        // double mLikelihoodArray[READ_NUM][HAPLOTYPE_NUM]);
						uint8* haplotype,
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
prob_t computeReadLikelihoodGivenHaplotypeLog10(uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH],
                                                uint8 readBases[MAX_READ_LENGTH],
                                                int haplotypeLength, // bases
                                                int readLength, // bases
                                                uint8 readQuals[MAX_READ_LENGTH],
                                                uint8 insertionGOP[MAX_READ_LENGTH],
                                                uint8 deletionGOP[MAX_READ_LENGTH],
                                                uint8 overallGCP[MAX_READ_LENGTH],
                                                bool recacheReadValues,
                                                uint8 nextHaplotypeBases[MAX_HAPLOTYPE_LENGTH],
												bool isNextHaplotypeExist,
                                                prob_t transition[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH],
                                                prob_t prior[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH]);
double subComputeReadLikelihoodGivenHaplotypeLog10(uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH],
                                                   uint8 readBases[MAX_READ_LENGTH],
                                                   int haplotypeLength,
                                                   int readLength,
                                                   uint8 readQuals[MAX_READ_LENGTH],
                                                   uint8 insertionGOP[MAX_READ_LENGTH],
                                                   uint8 deletionGOP[MAX_READ_LENGTH],
                                                   uint8 overallGCP[MAX_READ_LENGTH],
                                                   int hapStartIndex,
                                                   bool recacheReadValues,
                                                   int nextHapStartIndex,
                                                   prob_t transition[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH],
                                                   prob_t prior[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                   prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                   prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                   prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH]);
int findFirstPositionWhereHaplotypesDiffer(uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH], uint8 nextHaplotypeBases[MAX_HAPLOTYPE_LENGTH]);
void initializeMatrixValues(// uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH],
                            prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                            int haplotypeLength);
void initializePriors(uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH],
                      uint8 readBases[MAX_READ_LENGTH],
                      uint8 readQuals[MAX_READ_LENGTH],
                      int startIndex,
                      prob_t prior[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                      int haplotypeLength,
                      int readLength
                      );
prob_t qualToErrorProbLog10(uint8 qual);
prob_t qualToProbLog10(uint8 qual);
prob_t qualToErrorProb(prob_t qual);
void initializeProbabilities(uint8 insertionGOP[MAX_READ_LENGTH],
                             uint8 deletionGOP[MAX_READ_LENGTH],
                             uint8 overallGCP[MAX_READ_LENGTH],
                             prob_t transition[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH]);
void qualToTransProbsLog10(prob_t dest[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH],
                           uint8 insQual[MAX_READ_LENGTH],
                           uint8 delQual[MAX_READ_LENGTH],
                           uint8 gcp[MAX_READ_LENGTH]
                           );
void updateCell(int indI, int indJ, prob_t* prior, prob_t transition[TRANS_PROB_ARRAY_LENGTH],
                prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH]);
prob_t myLog10SumLog10(prob_t var1, prob_t var2);
prob_t finalLikelihoodCalculation(prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                  prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                  int haplotypeLength,
                                  int readLength
                                  );
prob_t exp_cal(prob_t value);
prob_t log_cal(prob_t value);
prob_t log10_cal(prob_t value);
prob_t pow_cal(prob_t value1, prob_t value2);
prob_t min_cal(prob_t value1, prob_t value2);
prob_t matchToMatchProbLog10(uint8 insQual, uint8 delQual);

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
					prob_t* mLikelihoodArray)
{
    // prob_t mLikelihoodArray[READ_NUM][HAPLOTYPE_NUM];
    
    // -------------------------------
    // PairHMM.java --> initialize()
    // -------------------------------
    // In the original GATK tools, we will travese the haplotype list and read list to find out
    // the maxHaplotypeLength and maxReadLength
    // here, we will assign two variables with fixed values
    // maxHaplotypeLength = MAX_HAPLOTYPE_LENGTH;
    // maxReadLength = MAX_READ_LENGTH;
    // paddedHaplotypeLength = maxHaplotypeLength + 1;
    // paddedReadLength = maxReadLength + 1;
    
    // constantsAreInitialized = true;
    // initialized = true;
    
    // -------------------------------
    // PairHMM.java --> computeLikelihoods()
    // @param likelihoods where to store the likelihoods where position [a][r] is reversed for the likelihoods of {@code reads[r]} conditional to {@code alleles[a]}
    // @param processedReads reads to analyze instead of the ones present in the destination read-likelihoods
    // @param gcp penalty for gap continuation base array map for processed reads
    // -------------------------------
    computeLikelihoods(haplotype, read, haplotypeLength, readLength, haplotypeNum, readNum, readBaseQualities, readBaseInsertionQualities, readBaseDeletionQualities, readGCP, mLikelihoodArray);
    return;
}

void computeLikelihoods(uint8* haplotype,
        				uint8* read,
						int* haplotypeLength, // bases
						int* readLength, // bases
						int haplotypeNum,
						int readNum,
						uint8* readBaseQualities,
						uint8* readBaseInsertionQualities,
						uint8* readBaseDeletionQualities,
						uint8* readGCP,
						prob_t* mLikelihoodArray)
{
    int readIndex = 0;
    for (readIndex = 0; readIndex < readNum; readIndex++)
    {
        bool isFirstHaplotype = true;
        bool isNextHaplotypeExist = false;
        // -------------------------------
        // N2MemoryPairHMM.java --> initialize()
        // -------------------------------
        // a matrix of {@code maxReadLength + 1} by {@link #TRANS_PROB_ARRAY_LENGTH} positions.
        prob_t transition[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH];
        prob_t prior[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH];
        prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH];
        prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH];
        prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH];
        // -------------------------------
        // Log10PairHMM.java --> initialize()
        // -------------------------------
        for (int iii = 0; iii < PADDED_MAX_READ_LENGTH; iii++)
        {
            for (int jjj = 0; jjj < PADDED_MAX_HAPLOTYPE_LENGTH; jjj++)
            {
                matchMatrix[iii][jjj] = NEGATIVE_INFINITY;
                insertionMatrix[iii][jjj] = NEGATIVE_INFINITY;
                deletionMatrix[iii][jjj] = NEGATIVE_INFINITY;
            }
        }
        uint8 readCur[MAX_READ_LENGTH];
        uint8 readBaseQualitiesCur[MAX_READ_LENGTH];
        uint8 readBaseInsertionQualitiesCur[MAX_READ_LENGTH];
        uint8 readBaseDeletionQualitiesCur[MAX_READ_LENGTH];
        uint8 readGCPCur[MAX_READ_LENGTH];
        memcpy(readCur, read + readIndex * MAX_READ_LENGTH, MAX_READ_LENGTH);
        memcpy(readBaseQualitiesCur, readBaseQualities + readIndex * MAX_READ_LENGTH, MAX_READ_LENGTH);
        memcpy(readBaseInsertionQualitiesCur, readBaseInsertionQualities + readIndex * MAX_READ_LENGTH, MAX_READ_LENGTH);
        memcpy(readBaseDeletionQualitiesCur, readBaseDeletionQualities + readIndex * MAX_READ_LENGTH, MAX_READ_LENGTH);
        memcpy(readGCPCur, readGCP + readIndex * MAX_READ_LENGTH, MAX_READ_LENGTH);
        int readLengthCur = *(readLength + readIndex);
        for (int haplotypeIndex = 0; haplotypeIndex < haplotypeNum; haplotypeIndex++)
        {
        	uint8 haplotypeCur[MAX_HAPLOTYPE_LENGTH];
        	memcpy(haplotypeCur, haplotype + haplotypeIndex * MAX_HAPLOTYPE_LENGTH, MAX_HAPLOTYPE_LENGTH);
        	uint8 nextHaplotypeBases[MAX_HAPLOTYPE_LENGTH];
        	if (haplotypeIndex < HAPLOTYPE_NUM - 1){
        		memcpy(nextHaplotypeBases, haplotype + (haplotypeIndex + 1) * MAX_HAPLOTYPE_LENGTH, MAX_HAPLOTYPE_LENGTH);
        		isNextHaplotypeExist = true;
        	}
        	else {
        		isNextHaplotypeExist = false;
        	}
        	int haplotypeLengthCur = *(haplotypeLength + haplotypeIndex);
            prob_t likelihood = computeReadLikelihoodGivenHaplotypeLog10(haplotypeCur,
                                                                         readCur,
                                                                         haplotypeLengthCur,
                                                                         readLengthCur,
                                                                         readBaseQualitiesCur,
                                                                         readBaseInsertionQualitiesCur,
                                                                         readBaseDeletionQualitiesCur,
                                                                         readGCPCur,
                                                                         isFirstHaplotype,
																		 nextHaplotypeBases,
																		 isNextHaplotypeExist,
                                                                         transition,
                                                                         prior,
                                                                         matchMatrix,
                                                                         insertionMatrix,
                                                                         deletionMatrix);
            mLikelihoodArray[readIndex * HAPLOTYPE_NUM + haplotypeIndex] = likelihood;
            isFirstHaplotype = false;
        }
    }
}

// -------------------------------
// PairHMM.java --> computeReadLikelihoodGivenHaplotypeLog10()
// Compute the total probability of read arising from haplotypeBases given base substitution, insertion, and deletion probabilities
// @param haplotypeBases the full sequence of the haplotype
// @param readBases the bases of the read
// @param readQuals the phred-scaled per base substitution quality socres of read
// @param insertionGOP the phred-scaled per base insertion quality socres of read
// @param deletionGOP the phred-scaled per base deletion quality scores of read
// @param overallGCP the phred-scaled gap continuation penalties scores of read
// @param recacheReadValues if false, we don't recalculate any cached results, assuming that readBases and its associated parameters are the same, and only the haplotype bases are changing underneath us
// @return the log10 probability of read coming from the haplotype under the provided error model
// -------------------------------
prob_t computeReadLikelihoodGivenHaplotypeLog10(uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH],
                                                uint8 readBases[MAX_READ_LENGTH],
                                                int haplotypeLength, // bases
                                                int readLength, // bases
                                                uint8 readQuals[MAX_READ_LENGTH],
                                                uint8 insertionGOP[MAX_READ_LENGTH],
                                                uint8 deletionGOP[MAX_READ_LENGTH],
                                                uint8 overallGCP[MAX_READ_LENGTH],
                                                bool recacheReadValues,
                                                uint8 nextHaplotypeBases[MAX_HAPLOTYPE_LENGTH],
												bool isNextHaplotypeExist,
                                                prob_t transition[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH],
                                                prob_t prior[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH])
{
    static int hapStartIndex;
    hapStartIndex = (recacheReadValues) ? 0 : hapStartIndex;
    // Pre-compute the difference between the current haplotype and the next one to be run
    int nextHapStartIndex = (isNextHaplotypeExist == false)? 0 : findFirstPositionWhereHaplotypesDiffer(haplotypeBases, nextHaplotypeBases);
    		// (nextHaplotypeBases == NULL) ? 0 : findFirstPositionWhereHaplotypesDiffer(haplotypeBases, nextHaplotypeBases);
    prob_t result = subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases,
    															readBases,
																haplotypeLength,
																readLength,
																readQuals,
																insertionGOP,
																deletionGOP,
																overallGCP,
																hapStartIndex,
																recacheReadValues,
																nextHapStartIndex,
																transition,
																prior,
																matchMatrix,
																insertionMatrix,
																deletionMatrix);
    //uint8 previousHaplotypeBases = haplotypeBases;
    hapStartIndex = (nextHapStartIndex < hapStartIndex) ? 0 : nextHapStartIndex;
    return result;
}

int findFirstPositionWhereHaplotypesDiffer(uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH], uint8 nextHaplotypeBases[MAX_HAPLOTYPE_LENGTH])
{
    for (int iii = 0; iii < MAX_HAPLOTYPE_LENGTH; iii++)
    {
        if (haplotypeBases[iii] != nextHaplotypeBases[iii])
        {
            return iii;
        }
    }
    return MAX_HAPLOTYPE_LENGTH; // This is wrong, it should be max(len(hap1), len(hap2))
}

// -------------------------------
// Log10PairHMM.java --> subComputeReadLikelihoodGivenHaplotypeLog10()
// -------------------------------
double subComputeReadLikelihoodGivenHaplotypeLog10(uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH],
                                                   uint8 readBases[MAX_READ_LENGTH],
                                                   int haplotypeLength,
                                                   int readLength,
                                                   uint8 readQuals[MAX_READ_LENGTH],
                                                   uint8 insertionGOP[MAX_READ_LENGTH],
                                                   uint8 deletionGOP[MAX_READ_LENGTH],
                                                   uint8 overallGCP[MAX_READ_LENGTH],
                                                   int hapStartIndex,
                                                   bool recacheReadValues,
                                                   int nextHapStartIndex,
                                                   prob_t transition[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH],
                                                   prob_t prior[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                   prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                   prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                                   prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH])
{
    // --------------------------------
    // Log10PairHMM.java -> initializeProbabilities
    // --------------------------------
    if (recacheReadValues)
    {
        initializeProbabilities(insertionGOP, deletionGOP, overallGCP, transition);
    }
    static int previousHaplotypeLength = 0;
    if (recacheReadValues)
        initializeMatrixValues(deletionMatrix, haplotypeLength);
    else if (previousHaplotypeLength != haplotypeLength)
        initializeMatrixValues(deletionMatrix, haplotypeLength);
    
    initializePriors(haplotypeBases, readBases, readQuals, hapStartIndex, prior, haplotypeLength, readLength);
    int paddedHaplotypeLength, paddedReadLength;
    paddedHaplotypeLength = haplotypeLength + 1;
    paddedReadLength = readLength + 1;
    for (int i = 1; i < PADDED_MAX_READ_LENGTH; i++)
    {
        if (i < paddedReadLength)
        {
            for (int j = hapStartIndex + 1; j < PADDED_MAX_READ_LENGTH; j++)
            {
                if (j < paddedHaplotypeLength)
                {
                    updateCell(i, j, &prior[i][j], transition[i], matchMatrix, insertionMatrix, deletionMatrix);
                }
            }
        }

    }
    previousHaplotypeLength = haplotypeLength;
    // final probability is the log10 sum of the last element in the Match and Insertion state arrays
    // this way we ignore all paths that ended in deletions
    // but we have to sum all the paths ending in the M and I matrices, because they're no longer extended.
    return finalLikelihoodCalculation(matchMatrix, insertionMatrix, haplotypeLength, readLength);
 }

// -------------------------------
// Log10PairHMM.java --> initializeMatrixValues()
// -------------------------------
void initializeMatrixValues(// uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH],
                            prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                            int haplotypeLength)
{
    prob_t initialValue = log10_cal(1.0 / haplotypeLength);
    for (int j = 0; j < PADDED_MAX_HAPLOTYPE_LENGTH; j++)
    {
        deletionMatrix[0][j] = initialValue;
    }
}

// -------------------------------
// Log10PairHMM.java --> initializePriors()
// -------------------------------
void initializePriors(uint8 haplotypeBases[MAX_HAPLOTYPE_LENGTH],
                      uint8 readBases[MAX_READ_LENGTH],
                      uint8 readQuals[MAX_READ_LENGTH],
                      int startIndex,
                      prob_t prior[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                      int haplotypeLength,
                      int readLength
                      )
{
    prob_t log10_3 = log10_cal(3.0);
    for (int i = 0; i < MAX_READ_LENGTH; i++)
    {
        if (i < readLength)
        {
            uint8 x = readBases[i];
            uint8 qual = readQuals[i];
            for (int j = startIndex; j< MAX_HAPLOTYPE_LENGTH; j++)
            {
                if (j < haplotypeLength)
                {
                    uint8 y = haplotypeBases[j];
                    prior[i + 1][j + 1] = (x == y || x == (uint8)'N' || y == (uint8)'N' ?
                                           qualToProbLog10(qual) : (qualToErrorProbLog10(qual) - log10_3));

                }
            }
        }
    }
}

// -------------------------------
// QualityUtils.java --> qualToErrorProbLog10()
// -------------------------------
prob_t qualToErrorProbLog10(uint8 qual)
{
    prob_t qual_p = qual & 0xFF;
    return qual_p * -0.1;
}

prob_t qualToProbLog10(uint8 qual)
{
    static prob_t qualToErrorProbCache[MAX_QUAL+1];
    static prob_t qualToProbLog10Cache[MAX_QUAL+1];
    for (int i = 0; i < MAX_QUAL; i++)
    {
        qualToErrorProbCache[i] = qualToErrorProb((prob_t)i);
        qualToProbLog10Cache[i] = log10_cal(1.0 - qualToErrorProbCache[i]);
    }
    return qualToProbLog10Cache[(int)(qual & 0xFF)];
}

prob_t qualToErrorProb(prob_t qual)
{
    return pow_cal(10.0, (qual / -10.0));
}

// -------------------------------
// Log10PairHMM.java --> initializeProbabilities()
// -------------------------------
void initializeProbabilities(uint8 insertionGOP[MAX_READ_LENGTH],
                             uint8 deletionGOP[MAX_READ_LENGTH],
                             uint8 overallGCP[MAX_READ_LENGTH],
                             prob_t transition[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH])
{
    qualToTransProbsLog10(transition, insertionGOP, deletionGOP, overallGCP);
    // constantsAreInitialized = true;
}

// -------------------------------
// PairHMMModel.java --> qualToTransProbsLog10()
// @param insQual the insertion quality score as as byte.
// @param delQual the deletion quality score as a byte.
// @param gcp the gap-continuation-penalty score as a byte.
// -------------------------------
void qualToTransProbsLog10(prob_t dest[PADDED_MAX_READ_LENGTH][TRANS_PROB_ARRAY_LENGTH],
                           uint8 insQual[MAX_READ_LENGTH],
                           uint8 delQual[MAX_READ_LENGTH],
                           uint8 gcp[MAX_READ_LENGTH]
                           )
{
    for (int i = 0; i < MAX_READ_LENGTH; i++)
    {
        dest[i+1][MATCH_TO_MATCH] = matchToMatchProbLog10(insQual[i], delQual[i]);
        dest[i+1][MATCH_TO_INSERTION] = qualToErrorProbLog10(insQual[i]);
        dest[i+1][MATCH_TO_DELETION] = qualToErrorProbLog10(delQual[i]);
        dest[i+1][INDEL_TO_MATCH] = qualToProbLog10(gcp[i]);
        dest[i+1][INSERTION_TO_INSERTION] = qualToErrorProbLog10(gcp[i]);
        dest[i+1][DELETION_TO_DELETION] = dest[i][INSERTION_TO_INSERTION];
    }
}

// -------------------------------
// PairHMMModel.java --> matchToMatchProb()
// returns the probability that neither of two event takes place.
// -------------------------------
prob_t matchToMatchProbLog10(uint8 insQual, uint8 delQual)
{
    prob_t LN10 = log10_cal(10);
    prob_t INV_LN10 = 1.0 / LN10;
    prob_t matchToMatchLog10[(MAX_QUAL + 1) * (MAX_QUAL + 2) >> 1];
    int offset = 0;
    for (int i = 0; i <= MAX_QUAL; i++)
    {
    	offset += i;
        for (int j = 0; j <= i; j++)
        {
            // prob_t log10Sum = approximateLog10SumLog10(-.1 * i, -.1 * j);
            prob_t log10Sum = myLog10SumLog10(-.1 * i, -.1 * j);
            matchToMatchLog10[offset + j] = log10_cal((prob_t)1 + (prob_t)(-min_cal((prob_t)1, (prob_t)(pow_cal(10, log10Sum))))) * INV_LN10;
            // matchToMatchProb[offset + j] = pow(10, matchToMatchLog10[offset + j]);
        }
    }
    int insQual_int = (insQual & 0xFF);
    int delQual_int = (delQual & 0xFF);
    int minQual;
    int maxQual;
    if (insQual <= delQual)
    {
        minQual = insQual;
        maxQual = delQual;
    } else
    {
        minQual = delQual;
        maxQual = insQual;
    }
    /* return (MAX_QUAL < maxQual) ? log10_cal(1 + (-min(1, 10 ** approximateLog10SumLog10(-.1 * minQual, -.1 * maxQual)))) * INV_LN10 :
                                matchToMatchLog10[(maxQual * (maxQual + 1)) >> 1 + minQual]; */
    prob_t value = (MAX_QUAL < maxQual) ? log10_cal((prob_t)1 + (prob_t)(-min_cal((prob_t)1, (prob_t)(pow_cal(10, myLog10SumLog10(-.1 * minQual, -.1 * maxQual)))))) * INV_LN10 :
            matchToMatchLog10[(maxQual * (maxQual + 1)) >> 1 + minQual];
    return value;
}

prob_t pow_cal(prob_t value1, prob_t value2)
{
	prob_t res = exp_cal(value2 * log_cal(value1));
	return res;
}

prob_t min_cal(prob_t value1, prob_t value2)
{
	return fmin(value1, value2);
}

prob_t log_cal(prob_t value)
{
	return log(value);
	// return logf(value);
}

prob_t log10_cal(prob_t value)
{
    return log10(value); // double
    // return log10f(value); // float
}

prob_t exp_cal(prob_t value)
{
    return exp(value); // double
    // return expf(value); // float
}

// -------------------------------
// Log10PairHMM.java --> finalLikelihoodCalculation()
// -------------------------------
prob_t finalLikelihoodCalculation(prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                  prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                                  int haplotypeLength,
                                  int readLength
                                  )
{
    int endI = (readLength + 1) - 1;
    int paddedHaplotypeLength = haplotypeLength + 1;
    prob_t finalSumProbabilities = myLog10SumLog10(matchMatrix[endI][1], insertionMatrix[endI][1]);
    for (int j = 2; j < PADDED_MAX_HAPLOTYPE_LENGTH; j++)
    {
        if (j < paddedHaplotypeLength)
        {
            finalSumProbabilities += myLog10SumLog10(matchMatrix[endI][j], insertionMatrix[endI][j]);
        }
    }
    return finalSumProbabilities;
}

// -------------------------------
// Compute the log10SumLog10 of the values
// @param vlaues an array of log10 probabilities that need to be summed
// @return the lgo10 of the sum of the probabilities
// -------------------------------
prob_t myLog10SumLog10(prob_t var1, prob_t var2)
{
    prob_t res1 = pow_cal(10, var1);
    prob_t res2 = pow_cal(10, var2);
    prob_t value = log10_cal(res1 + res2);
    return value;
}

void updateCell(int indI, int indJ, prob_t* prior, prob_t transition[TRANS_PROB_ARRAY_LENGTH],
                prob_t matchMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                prob_t insertionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH],
                prob_t deletionMatrix[PADDED_MAX_READ_LENGTH][PADDED_MAX_HAPLOTYPE_LENGTH])
{
    matchMatrix[indI][indJ] = *prior +
        myLog10SumLog10(matchMatrix[indI - 1][indJ - 1] + transition[MATCH_TO_MATCH],
                        myLog10SumLog10(insertionMatrix[indI - 1][indJ - 1] + transition[INDEL_TO_MATCH], deletionMatrix[indI - 1][indJ - 1] + transition[INDEL_TO_MATCH]));
    insertionMatrix[indI][indJ] = myLog10SumLog10(matchMatrix[indI - 1][indJ] + transition[MATCH_TO_INSERTION], insertionMatrix[indI - 1][indJ] + transition[INSERTION_TO_INSERTION]);
    deletionMatrix[indI][indJ] = myLog10SumLog10(matchMatrix[indI][indJ - 1] + transition[MATCH_TO_DELETION], deletionMatrix[indI][indJ - 1] + transition[DELETION_TO_DELETION]);
}