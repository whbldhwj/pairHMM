#include "pairHMM.h"
#include "stdio.h"
#include "stdlib.h"

int main()
{
    // uint8 haplotype[HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];
    // uint8 read[READ_NUM][MAX_READ_LENGTH];
    // int haplotypeLength[HAPLOTYPE_NUM]; // bases
    // int readLength[READ_NUM]; // bases
    // uint8 readBaseQualities[READ_NUM][MAX_READ_LENGTH];
    // uint8 readInsertionQualities[READ_NUM][MAX_READ_LENGTH];
    // uint8 readDeletionQualities[READ_NUM][MAX_READ_LENGTH];
    // uint8 readGCP[READ_NUM][MAX_READ_LENGTH];
    // prob_t mLikelihoodArray_sw[READ_NUM][HAPLOTYPE_NUM];
    // prob_t mLikelihoodArray_hw[READ_NUM][HAPLOTYPE_NUM];

	uint8* haplotype;
	uint8* read;
	int* haplotypeLength; // bases
	int* readLength; // bases
	uint8* readBaseQualities;
	uint8* readInsertionQualities;
	uint8* readDeletionQualities;
	uint8* readGCP;
	prob_t* mLikelihoodArray_sw;
	prob_t* mLikelihoodArray_hw;

    haplotype = (uint8*) malloc(HAPLOTYPE_NUM * MAX_HAPLOTYPE_LENGTH);
    read = (uint8*) malloc(READ_NUM * MAX_READ_LENGTH);
    haplotypeLength = (int*) malloc(HAPLOTYPE_NUM * sizeof(int));
    readLength = (int*) malloc(READ_NUM * sizeof(int));
    readBaseQualities = (uint8*) malloc(READ_NUM * MAX_READ_LENGTH);
    readInsertionQualities = (uint8*) malloc(READ_NUM * MAX_READ_LENGTH);
    readDeletionQualities = (uint8*) malloc(READ_NUM * MAX_READ_LENGTH);
    readGCP = (uint8*) malloc(READ_NUM * MAX_READ_LENGTH);
    mLikelihoodArray_sw = (prob_t*) malloc(READ_NUM * HAPLOTYPE_NUM * sizeof(prob_t));
    mLikelihoodArray_hw = (prob_t*) malloc(READ_NUM * HAPLOTYPE_NUM * sizeof(prob_t));

    // Initialization
    for (int i = 0; i < HAPLOTYPE_NUM; i++){
    	for (int j = 0; j < MAX_HAPLOTYPE_LENGTH; j++){
    		haplotype[i * MAX_HAPLOTYPE_LENGTH + j] = (uint8)'I';
    	}
    }
    for (int i = 0; i < READ_NUM; i++){
    	for (int j = 0; j< MAX_READ_LENGTH; j++){
    		read[i * MAX_READ_LENGTH + j] = (uint8)'I';
    	}
    }
    for (int i = 0; i < HAPLOTYPE_NUM; i++){
    	haplotypeLength[i] = 0;
    }
    for (int i = 0; i < READ_NUM; i++){
    	readLength[i] = 0;
    }
    for (int i = 0; i < READ_NUM; i++){
    	for (int j = 0; j < MAX_READ_LENGTH; j++){
    		readBaseQualities[i * MAX_READ_LENGTH + j] = (uint8)'I';
    		readInsertionQualities[i * MAX_READ_LENGTH + j] = (uint8)'I';
    		readDeletionQualities[i * MAX_READ_LENGTH + j] = (uint8)'I';
    		readGCP[i * MAX_READ_LENGTH + j] = (uint8)'I';
    	}
    }
    for(int i = 0; i <READ_NUM; i++) {
    	for (int j = 0; j < HAPLOTYPE_NUM; j++) {
    		mLikelihoodArray_sw[i * HAPLOTYPE_NUM + j] = (prob_t)0;
    		mLikelihoodArray_hw[i * HAPLOTYPE_NUM + j] = (prob_t)0;
    	}
    }
    // Read data
    FILE* fp;
    fp = fopen("./haplotype.txt", "r");
    if (fp == NULL){
    	printf("Error, cannot open file!\n");
    	return -1;
    }
    char base = fgetc(fp);
    int length = 0;
    int haplotypeNum = 0;
    int readNum = 0;
    int num = 0;
    while(!feof(fp)){
    	if (base == '\n'){
    		haplotypeLength[num] = length;
    		length = 0;
    		num++;
    		haplotypeNum++;
    	}
    	else {
    		haplotype[num * MAX_HAPLOTYPE_LENGTH + length] = (uint8)base;
    		length++;
    	}
    }
    fclose(fp);

    fp = fopen("./read.txt", "r");
    if (fp == NULL){
    	printf("Error, cannot open file!\n");
    	return -1;
    }
    base = fgetc(fp);
    length = 0;
    num = 0;
    while(!feof(fp)){
    	if (base == '\n'){
    		readLength[num] = length;
    		length = 0;
    		num++;
    		readNum++;
    	}
    	else {
    		read[num * MAX_READ_LENGTH + length] = (uint8)base;
    		length++;
    	}
    }
    fclose(fp);

    fp = fopen("./readBaseQualities.txt", "r");
    if (fp == NULL){
    	printf("Error, cannot open file!\n");
    	return -1;
    }
    base = fgetc(fp);
    length = 0;
    num = 0;
    while(!feof(fp)){
    	if (base == '\n'){
    		length = 0;
    		num++;
    	}
    	else {
    		readBaseQualities[num * MAX_READ_LENGTH + length] = (uint8)base;
    		length++;
    	}
    }
    fclose(fp);

    fp = fopen("./readInsertionQualities.txt", "r");
    if (fp == NULL){
    	printf("Error, cannot open file!\n");
    	return -1;
    }
    base = fgetc(fp);
    length = 0;
    num = 0;
    while(!feof(fp)){
    	if (base == '\n'){
    		length = 0;
    		num++;
    	}
    	else {
    		readInsertionQualities[num * MAX_READ_LENGTH + length] = (uint8)base;
    		length++;
    	}
    }
    fclose(fp);

    fp = fopen("./readDeletionQualities.txt", "r");
    if (fp == NULL){
    	printf("Error, cannot open file!\n");
    	return -1;
    }
    base = fgetc(fp);
    length = 0;
    num = 0;
    while(!feof(fp)){
    	if (base == '\n'){
    		length = 0;
    		num++;
    	}
    	else {
    		readDeletionQualities[num * MAX_READ_LENGTH + length] = (uint8)base;
    		length++;
    	}
    }
    fclose(fp);

    fp = fopen("./readGCP.txt", "r");
    if (fp == NULL){
    	printf("Error, cannot open file!\n");
    	return -1;
    }
    base = fgetc(fp);
    length = 0;
    num = 0;
    while(!feof(fp)){
    	if (base == '\n'){
    		length = 0;
    		num = num+1;
    	}
    	else {
    		readGCP[num * MAX_READ_LENGTH + length] = (uint8)base;
    		length++;
    	}
    }
    fclose(fp);

    fp = fopen("./mLikelihoodArray.txt", "r");
    if (fp == NULL){
    	printf("Error, cannot open file!\n");
    	return -1;
    }
    double data;
    fscanf(fp, "%lf", &data);
    length = 0;
    num = 0;
    while(!feof(fp)){
    	if (base == '\n'){
    		length = 0;
    		num++;
    	}
    	else {
    		mLikelihoodArray_sw[num * HAPLOTYPE_NUM + length] = data;
    		length++;
    	}
    }
    fclose(fp);

#ifdef DEBUG
    // uint8 haplotype[HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];
    // uint8 read[READ_NUM][MAX_READ_LENGTH];
    // int haplotypeLength[HAPLOTYPE_NUM]; // bases
    // int readLength[READ_NUM]; // bases
    // uint8 readBaseQualities[READ_NUM][MAX_READ_LENGTH];
    // uint8 readInsertionQualities[READ_NUM][MAX_READ_LENGTH];
    // uint8 readDeletionQualities[READ_NUM][MAX_READ_LENGTH];
    // uint8 readGCP[READ_NUM][MAX_READ_LENGTH];
    // prob_t mLikelihoodArray_sw[READ_NUM][HAPLOTYPE_NUM];
    // prob_t mLikelihoodArray_hw[READ_NUM][HAPLOTYPE_NUM];
    for (int i = 0; i < haplotypeNum; i++){
    	for (int j = 0; j < haplotypeLength[i]; j++){
    		printf("%c ", haplotype[i * MAX_HAPLOTYPE_LENGTH + j]);
    	}
    	printf("\n");
    }
#endif

    pairHMM_kernel(haplotype, // byte
                   read,
                   haplotypeLength, // bases
                   readLength, // bases
				   haplotypeNum,
				   readNum,
                   readBaseQualities,
                   readInsertionQualities,
                   readDeletionQualities,
                   readGCP,
                   mLikelihoodArray_hw
                   );
    
    int err = 0;

    int i, j;
    for (i = 0; i < READ_NUM; i++)
    {
        for (j = 0; j < HAPLOTYPE_NUM; j++)
        {
            if (fabs(mLikelihoodArray_sw[i * HAPLOTYPE_NUM + j] - mLikelihoodArray_hw[i * HAPLOTYPE_NUM + j]) > 1e-15)
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
