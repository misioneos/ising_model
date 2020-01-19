#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ising.h"

int ** spinGen(int s){
    int **out;
    out = (int**)malloc(sizeof(int*) * s);
    for (int i = 0; i < s; ++i) {
        out[i] = malloc(sizeof(int) * s);
    }
    for(int i=0; i<s; ++i)
    {
        for(int j=0; j<s; ++j)
        {
            out[i][j] = rand()%2;
        }
    }
    return out;
}


//Srednia magnetyzacja
double magnetization(int (**spin), int s)
{
    int M = 0;
    for (int i = 0; i < s; i++)
    {
        for (int j = 0; j < s; j++)
        {
        if (spin[i][j]) M++;
        else M--;
        }
    }
    return 1.0 * M / (s * s);
}

//Energia spinu
double site_energy(int i, int j, int(**spin), double J, int s)
{
    double E = 0.0;
    i = i % s;
    j = j % s;
    if (spin[i][j] == spin[(i+1)%s][j]) E++;
    else E--;
    if (spin[i][j] == spin[i][(j+1)%s]) E++;
    else E--;
    if (spin[i][j] == spin[(i-1+N)%s][j]) E++;
    else E--;
    if (spin[i][j] == spin[i][(j-1+s)%s]) E++;
    else E--;
    return -J * E;
}

//Energia calkowita
double total_energy(int(**spin), double J, int s)
{
    double E = 0.0;
    for (int i = 0; i < s; i++)
    for (int j = 0; j < s; j++)
        E += 0.5 * site_energy(i, j, spin, J, s);
    return E;
}

void step(int i, int j, int(**spin), double J, double kT, int s)
{
    i = i % s;
    j = j % s;
    double deltaE = -2 * site_energy(i, j, spin, J, s);
    if ((deltaE <= 0) || ((double)rand()/RAND_MAX <= exp(-deltaE / kT)))
    {
        if(spin[i][j]==0) spin[i][j]=1;
        else spin[i][j]=0;
    }
}

void phase(int(**spin), double J, double kT, int s)
{
    for (int steps = 0; steps < 1000000; steps++)
        {
        int i = (int) ((double)rand()/RAND_MAX * s);
        int j = (int) ((double)rand()/RAND_MAX * s);
        step(i, j, spin, J, kT, s);
        }
}

void testN6()
{
    double J=1.0;
    double kT=0.0;
    double mag[26];
    double x,y;
    int **spin;
    spin=spinGen(N);
    char *output_filename="dataN6.txt";
    FILE *output_file;
    output_file = fopen(output_filename, "w");
    if (output_file  == NULL)
    {
        fprintf(stderr, "Nie moge otworzyc %s\n", output_filename);
        getchar();
        return 1;
    }
    else
     {
        fprintf(output_file,"%s %s      %s \n","#","x","y");
        y = 1.0;
        x = 0.0;
        fprintf(output_file,"% 6.2f % 6.2f \n",x,y);
        for(int i=0;i<26;i++)
        {
            phase(spin, J, kT, N);
            mag[i] = magnetization(spin, N);
            if (mag[i]<0) mag[i]=mag[i]*(-1.0);
            x = kT;
            y = mag[i];
            fprintf(output_file,"% 6.2f % 6.2f \n",x,y);
            kT = kT + 0.2;
        }
    }
    fclose(output_file);

    fprintf(stderr,"file saved");
    getchar();
}

void testN10()
{
    double J=1.0;
    double kT=0.0;
    double mag[26];
    double x,y;
    int **spin;
    spin=spinGen(O);
    char *output_filename="dataN10.txt";
    FILE *output_file;
    output_file = fopen(output_filename, "w");
    if (output_file  == NULL)
    {
        fprintf(stderr, "Nie moge otworzyc %s\n", output_filename);
        getchar();
        return 1;
    }
    else
     {
        fprintf(output_file,"%s %s      %s \n","#","x","y");
        y = 1.0;
        x = 0.0;
        fprintf(output_file,"% 6.2f % 6.2f \n",x,y);
        for(int i=0;i<26;i++)
        {
            phase(spin, J, kT, O);
            mag[i] = magnetization(spin, O);
            if (mag[i]<0) mag[i]=mag[i]*(-1.0);
            x = kT;
            y = mag[i];
            fprintf(output_file,"% 6.2f % 6.2f \n",x,y);
            kT = kT + 0.2;
        }
    }
    fclose(output_file);

    fprintf(stderr,"file saved");
    getchar();
}

void testN20()
{
    double J=1.0;
    double kT=0.0;
    double mag[26];
    double x,y;
    int **spin;
    spin=spinGen(P);
    char *output_filename="dataN20.txt";
    FILE *output_file;
    output_file = fopen(output_filename, "w");
    if (output_file  == NULL)
    {
        fprintf(stderr, "Nie moge otworzyc %s\n", output_filename);
        getchar();
        return 1;
    }
    else
     {
        fprintf(output_file,"%s %s      %s \n","#","x","y");
        y = 1.0;
        x = 0.0;
        fprintf(output_file,"% 6.2f % 6.2f \n",x,y);
        for(int i=0;i<26;i++)
        {
            phase(spin, J, kT, P);
            mag[i] = magnetization(spin, P);
            if (mag[i]<0) mag[i]=mag[i]*(-1.0);
            x = kT;
            y = mag[i];
            fprintf(output_file,"% 6.2f % 6.2f \n",x,y);
            kT = kT + 0.2;
        }
    }
    fclose(output_file);

    fprintf(stderr,"file saved");
    getchar();
}



void dependence(int (**spin), double J, int s)
{
    double mag[18];
    //double kT[]={1.0,1.2,1.4,1.6,1.8,1.85,1.9,1.95,2.0,2.05,2.1,2.15,2.2,2.25,2.3,2.5,2.7};
    double kT[]={1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7};
    for(int k=0; k<18; k++)
    {
        printf("\nkT=%f\t", kT[k]);
        phase(spin, J, kT[k], s);
        mag[k] = magnetization(spin, s);
        printf("M=%f", mag[k]);
    }
}


void test1()
{
    int **spinN, **spinO, **spinP, **spinR, **spinS;
    double J=1.0;
    spinN = spinGen(N);
    spinO = spinGen(O);
    spinP = spinGen(P);
    spinR = spinGen(R);
    spinS = spinGen(S);
    printf("6x6");
    dependence(spinN, J, N);
    printf("\n\n10x10");
    dependence(spinO, J, O);
    printf("\n\n20x20");
    dependence(spinP, J, P);
    printf("\n\n40x40");
    dependence(spinR, J, R);
    printf("\n\n100x100");
    dependence(spinS, J, S);
}

void test2()
{
    double J=1.0;
    double kT[]={1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7};
    double magN[18], magO[18], magP[18], magR[18], magS[18];
    int spinNtab[18], spinOtab[18], spinPtab[18], spinRtab[18], spinStab[18];

    for(int i=0; i<18; i++)
    {
        spinNtab[i]=spinGen(N);
        spinOtab[i]=spinGen(O);
        spinPtab[i]=spinGen(P);
        spinRtab[i]=spinGen(R);
        spinStab[i]=spinGen(S);
        printf("\nkT=%f\t", kT[i]);
        phase(spinNtab[i], J, kT[i], N);
        phase(spinOtab[i], J, kT[i], O);
        phase(spinPtab[i], J, kT[i], P);
        phase(spinRtab[i], J, kT[i], R);
        phase(spinStab[i], J, kT[i], S);
        magN[i] = magnetization(spinNtab[i], N);
        magO[i] = magnetization(spinOtab[i], O);
        magP[i] = magnetization(spinPtab[i], P);
        magR[i] = magnetization(spinRtab[i], R);
        magS[i] = magnetization(spinStab[i], S);
        printf("M(6)=%f\tM(10)=%f\tM(20)=%f\tM(40)=%f\tM(100)=%f\t", magN[i], magO[i], magP[i], magR[i], magS[i]);
    }
}
