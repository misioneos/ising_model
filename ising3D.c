#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ising3D.h"

int *** spinGen3D(int s){
    int ***out;
    out = malloc(sizeof(int**) * s);
    for (int i = 0; i < s; ++i) {
        out[i] = malloc(sizeof(int*) * s);
        for(int j = 0; j < s; j++)
            out[i][j] =  malloc(sizeof(int) * s);
    }
    for(int i=0; i<s; ++i)
    {
        for(int j=0; j<s; ++j)
        {
            for(int k = 0; k < s; ++k)
                out[i][j][k] = rand()%2;
        }
    }
    return out;
}

//Srednia magnetyzacja
double magnetization3D(int (***spin), int s)
{
    int M = 0;
    for (int i = 0; i < s; i++)
    {
        for (int j = 0; j < s; j++)
        {
            for (int k = 0; k < s; k++)
            {
                if (spin[i][j][k]) M++;
                else M--;
            }
        }
    }
    return 1.0 * M / (s * s * s);
}

//Energia spinu
double site_energy3D(int i, int j, int k, int(***spin), double J, int s)
{
    double E = 0.0;
    i = i % s;
    j = j % s;
    k = k % s;
    if (spin[i][j][k] == spin[(i+1)%s][j][k]) E++;
    else E--;
    if (spin[i][j][k] == spin[i][(j+1)%s][k]) E++;
    else E--;
    if (spin[i][j][k] == spin[i][j][(k+1)%s]) E++;
    else E--;
    if (spin[i][j][k] == spin[(i-1+N)%s][j][k]) E++;
    else E--;
    if (spin[i][j][k] == spin[i][(j-1+s)%s][k]) E++;
    else E--;
    if (spin[i][j][k] == spin[i][j][(k-1+s)%s]) E++;
    else E--;
    return -J * E;
}

//Energia calkowita
double total_energy3D(int(***spin), double J, int s)
{
    double E = 0.0;
    for (int i = 0; i < s; i++)
    for (int j = 0; j < s; j++)
    for (int k = 0; k < s; k++)
        E += 0.5 * site_energy3D(i, j, k, spin, J, s);
    return E;
}

void step3D(int i, int j, int k, int(***spin), double J, double kT, int s)
{
    i = i % s;
    j = j % s;
    k = k % s;
    double deltaE = -2 * site_energy3D(i, j, k, spin, J, s);
    if ((deltaE <= 0) || ((double)rand()/RAND_MAX <= exp(-deltaE / kT)))
    {
        if(spin[i][j][k]==0) spin[i][j][k]=1;
        else spin[i][j][k]=0;
    }
}

void phase3D(int(***spin), double J, double kT, int s)
{
    for (int steps = 0; steps < 5000000; steps++)
        {
        int i = (int) ((double)rand()/RAND_MAX * s);
        int j = (int) ((double)rand()/RAND_MAX * s);
        int k = (int) ((double)rand()/RAND_MAX * s);
        step3D(i, j, k, spin, J, kT, s);
        }
}

void testN63D()
{
    double J=1.0;
    double kT=0.0;
    double mag[26];
    double x,y;
    int ***spin;
    spin=spinGen3D(N);
    char *output_filename="dataN63D.txt";
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
            phase3D(spin, J, kT, N);
            mag[i] = magnetization3D(spin, N);
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

void testN103D()
{
    double J=1.0;
    double kT=0.0;
    double mag[26];
    double x,y;
    int ***spin;
    spin=spinGen3D(O);
    char *output_filename="dataN103D.txt";
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
            phase3D(spin, J, kT, O);
            mag[i] = magnetization3D(spin, O);
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

void testN203D()
{
    double J=1.0;
    double kT=0.0;
    double mag[26];
    double x,y;
    int ***spin;
    spin=spinGen3D(P);
    char *output_filename="dataN203D.txt";
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
            phase3D(spin, J, kT, P);
            mag[i] = magnetization3D(spin, P);
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
