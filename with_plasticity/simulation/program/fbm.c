#include <stdio.h>
#include <stdlib.h>

#define N 10000   // number of fibers
#define E 1.0     // Young's modulus
#define e_y 0.5   // mean value of yielding thresholds
#define e_b 1.7   // mean value of breaking thresholds
#define de_y 0.3  // radius of the uniform distribution of e_y
#define de_b 0.5  // radius of the uniform distribution of e_b
#define alpha 0.2 // ratio of the elastic modulus of yielding fibers to E

double elastic[N];  // intact fibers
double yielding[N]; // yielding thresholds of fibers
double broken[N];   // breaking thresholds of fibers

double rand01();
int CmpFunc(const void *_a, const void *_b);
void uniform(double array[], double mean, double radius, int sizeOfArray);
// initialize and sort arrays of uniform distribution

int main()
{
    FILE *constit;
    constit = fopen("../data/constit.txt", "w");

    // intialize yielding[] and broken[]
    uniform(yielding, e_y, de_y, N);
    uniform(broken, e_b, de_b, N);

    //-------------------- implement the constitutive model -------------------
    double eps = 0.0; // strain of the whole bundle
    for (int i = 0; i < N; i++)
    {
        printf("%d\n", i);
        double sigma = 0.0; // stress of the whole bundle
        double force = 0.0; // force on the whole bundle
        eps = yielding[i];  // let eps jump to the smallest yielding threshold
        for (int j = 0; j < i; j++)
        {
            force += E * yielding[j] + alpha * E * (eps - yielding[j]);
            // force on the yielding fibers
        }
        force += E * eps * (N - i); // force on the remaining intact fibers
        sigma = force / N;          // stress of the whole bundle
        fprintf(constit, "%lf\t%lf\n", eps, sigma);
    }
    //-------------------------------------------------------------------------
    fclose(constit);

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
