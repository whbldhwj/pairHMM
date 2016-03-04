//; $PE_NUM = param_define("PE_NUM", 41);
//; $D_W = param_define("D_W", 32);
#include <updateCell.h>

void updateCell()
{
	
#pragma HLS DATAFLOW

	PE_1();
	//; for ($i = 1 + 1; $i < $PE_NUM - 1 + 1; $i++) {
	PE_`$i`();
	//; }
	PE_`$PE_NUM`();
}

void PE_1()
{
	
}

//; for ($i = 1 + 1; $i < $PE_NUM - 1 + 1; $i++) {
void PE_`$i`()
{
	
}
//; }

void PE_`$PE_NUM`()
{
	
}