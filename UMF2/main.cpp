#include "SLAU.h"

using namespace std;



int main()
{
	SLAU<double> slv = SLAU<double>("grid", "FirstCondition", "time");
	slv.do_smth();
	system("pause");
	return 0;
}
