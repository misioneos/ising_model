#pragma once
#define N 6
#define O 10
#define P 20
#define R 40
#define S 100


//spis tresci dla 3D
int *** spinGen3D(int s);
double magnetization3D(int (***spin), int s);
double site_energy3D(int i, int j, int k, int(***spin), double J, int s);
double total_energy3D(int(***spin), double J, int s);
void step3D(int i, int j, int k, int(***spin), double J, double kT, int s);
void phase3D(int(***spin), double J, double kT, int s);

//testy koncowe
void testN63D();
void testN103D();
void testN203D();

void dependence3D(int (***spin), double J, int s);
