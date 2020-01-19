#include <stdlib.h>
#include <stdio.h>
#include "ising.h"
#include "ising3D.h"


int main()
{
    srand((unsigned) time(NULL));
    testN6();
    testN10();
    testN20();
    testN63D();
    testN103D();
    testN203D();
    return 0;
}





