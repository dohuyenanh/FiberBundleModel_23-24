#include <stdio.h>
#include <stdlib.h>
#include <string.h> //for memset

#define N 10000    // number of fibers
#define E 1.0      // Young's modulus
#define e_y 0.35   // mean value of yielding thresholds
#define e_b 3.5    // mean value of breaking thresholds
#define de_y 0.2   // radius of the uniform distribution of e_y
#define de_b_1 1.1 // radius of the uniform distribution of e_b
#define de_b_2 1.2
#define de_b_3 1.3
#define alpha_1 0.2 // ratio of the elastic modulus of yielding fibers to E
#define alpha_2 0.5
#define alpha_3 0.8
#define limit 5 // the upper limit of the strain

// double elastic[N];    // elastic fibers
double yielding[N]; // yielding thresholds of fibers
double broken[N];   // breaking thresholds of fibers
int i, j;           // counters

double rand01();

int CmpFunc(const void *_a, const void *_b);

void uniform(double array[], double mean, double radius, int sizeOfArray);
// initialize and sort arrays of uniform distribution

void constitutive(double e_yielding, double de_yielding, double e_breaking, double de_breaking, double alpha, FILE *outputFile);
// function to export output data of the constitutive relation as one of the parameter change

int main()
{
    char filename[50];

    double alpha[3] = {alpha_1, alpha_2, alpha_3};
    double de_breaking[3] = {de_b_1, de_b_2, de_b_3};

    for (int index = 0; index < 3; index++)
    {
        // intialize yielding[] and broken[]
        uniform(yielding, e_y, de_y, N);
        uniform(broken, e_b, de_b_1, N);

        // name for output file
        sprintf(filename, "../data/constit(alpha-%d).txt", index + 1);
        FILE *constit;
        constit = fopen(filename, "w");

        constitutive(e_y, de_y, e_b, de_breaking[0], alpha[index], constit);
    }

    for (int index = 0; index < 3; index++)
    {
        // intialize yielding[] and broken[]
        uniform(yielding, e_y, de_y, N);
        uniform(broken, e_b, de_breaking[index], N);

        // name for output file
        sprintf(filename, "../data/constit(de_b-%d).txt", index + 1);
        FILE *constit;
        constit = fopen(filename, "w");
        
        constitutive(e_y, de_y, e_b, de_breaking[index], alpha[0], constit);
    }

    return 0;
}

double rand01()
{
    return (double)rand() / (double)RAND_MAX;
}

int CmpFunc(const void *_a, const void *_b)
{
    const double *a = (const double *)_a;
    const double *b = (const double *)_b;

    if (*a > *b)
        return 1; // first item is bigger than the second one -> return 1
    else if (*a == *b)
        return 0; // equality -> return 0
    else
        return -1; // second item is bigger than the first one -> return -1
}

void uniform(double array[], double mean, double radius, int sizeOfArray)
{
    for (int i = 0; i < sizeOfArray; i++)
    {
        double r = rand01();
        array[i] = r * 2 * radius + mean - radius;
    }
    qsort(array, sizeOfArray, sizeof(array[0]), CmpFunc);
}

void constitutive(double e_yielding, double de_yielding, double e_breaking, double de_breaking, double alpha, FILE *outputFile)
{
    double eps = 0.0;   // strain of the whole bundle
    double force = 0.0; // force on the whole bundle
    double sigma = 0.0; // stress of the whole bundle

    // elastic-only regime ----------
    for (eps = 0; eps < yielding[0]; eps += yielding[0] / (N * 1.0 / 3)) // shouldn't it be eps < (e_y - de_y)?
    {
        sigma = E * eps;
        fprintf(outputFile, "%lf\t%lf\n", eps, sigma);
    }

    // yielding + elastic regime ----------------
    for (i = 0; i < N; i++)
    {
        printf("%d\n", i);
        sigma = 0.0;       // stress of the whole bundle
        force = 0.0;       // force on the whole bundle, assuming each fiber is of a unit cross-section
        eps = yielding[i]; // let strain jump to the smallest yielding threshold
        for (j = 0; j < i; j++)
        {
            force += E * yielding[j] + alpha * E * (eps - yielding[j]);
            // force on the yielding fibers
        }
        force += E * eps * (N - i); // force on the remaining elastic fibers
        sigma = force / N;          // stress of the whole bundle
        fprintf(outputFile, "%lf\t%lf\n", eps, sigma);
    }

    // yielding-only regime -----------
    for (eps = yielding[N - 1]; eps < broken[0]; eps += (broken[0] - yielding[N - 1]) / (N / 2.5)) // shouldn't it be eps < (e_b - de_b)?
    {
        force = 0.0;
        for (i = 0; i < N; i++)
            force += E * yielding[i] + alpha * E * (eps - yielding[i]);
        sigma = force / N;
        fprintf(outputFile, "%lf\t%lf\n", eps, sigma);
    }

    // yielding+broken regime ---------
    for (i = 0; i < N; i++)
    {
        force = 0.0;
        eps = broken[i];
        for (j = i; j < N; j++)
            force += E * yielding[j] + alpha * E * (eps - yielding[j]);
        sigma = force / N;
        fprintf(outputFile, "%lf\t%lf\n", eps, sigma);
    }

    // broken-only regime ------------
    for (eps = broken[N - 1]; eps < limit; eps += (limit - broken[N - 1]) / (N * 1.0 / 4.5))
    {
        sigma = 0.0;
        fprintf(outputFile, "%lf\t%lf\n", eps, sigma);
    }
    fclose(outputFile);
}
