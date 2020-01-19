#pragma once
#define N 6
#define O 10
#define P 20
#define R 40
#define S 100

//spis tresci
int ** spinGen(int s);
double magnetization(int (**spin), int s);
double site_energy(int i, int j, int(**spin), double J, int s);
double total_energy(int(**spin), double J, int s);
void step(int i, int j, int(**spin), double J, double kT, int s);
void phase(int(**spin), double J, double kT, int s);

//testy smieciowe
void dependence(int (**spin), double J, int s);
void test1();
void test2();

//testy koncowe
void testN6();
void testN10();
void testN20();
