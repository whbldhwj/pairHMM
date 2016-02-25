#include "pairHMM.h"

int main()
{
    uint8_t haplotype[HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];
    uint8_t read[READ_NUM][MAX_READ_LENGTH];
    int haplotypeLength[HAPLOTYPE_NUM]; // bases
    int readLength[READ_NUM]; // bases
    uint8_t readBaseQualities[READ_NUM][MAX_READ_LENGTH];
    uint8_t readBaseInsertionQualities[READ_NUM][MAX_READ_LENGTH];
    uint8_t readDeletionQualities[READ_NUM][MAX_READ_LENGTH];
    uint8_t readGCP[READ_NUM][MAX_READ_LENGTH];
    prob_t mLikelihoodArray_sw[READ_NUM][HAPLOTYPE_NUM];
    prob_t mLikelihoodArray_hw[READ_NUM][HAPLOTYPE_NUM];
    
#pragma cmost task name = "pairHMM"
    pairHMM_kernel(haplotype, // byte
                   read,
                   haplotypeLength, // bases
                   readLength, // bases
                   readBaseQualities,
                   readBaseInsertionQualities,
                   readDeletionQualities,
                   readGCP,
                   mLikelihoodArray_hw
                   );
    
    int err = 0;
    pairHMM_kernel_sw(haplotype, // byte
                      read,
                      haplotypeLength, // bases
                      readLength, // bases
                      readBaseQualities,
                      readBaseInsertionQualities,
                      readDeletionQualities,
                      readGCP,
                      mLikelihoodArray_sw
                      );
    int i, j;
    for (i = 0; i < READ_NUM; i++)
    {
        for (j = 0; j < HAPLOTYPE_NUM; j++)
        {
            if (fabs(mLikelihoodArray_sw[i][j] - mLikelihoodArray_hw[i][j]) > 1e-15)
                err ++;
        }
    }
    
    if (err > 0)
    {
        printf("Error numbers = %d\n", err);
    }
    else
    {
        printf("The test passed successfully.\n");
    }
    return err;
}