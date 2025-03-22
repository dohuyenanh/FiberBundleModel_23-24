#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 10000
#define K 1000
#define E 1.0
#define alpha 0.095
#define e_y 0.45
#define de_y 0.1
#define e_b 5.9
#define de_b 4.0
#define limit 11
#define binAval 40 // number of bins for avalanche size distribution (log)

double yielding[N];     // yielding thresholds of fibers
double broken[N];       // breaking thresholds of fibers
double s[N];            // stress of the whole bundle
double sigmaC[K];       // critical/largest stress where catastrophic avalanche does not occur after that
double epsC[K];         // critical/largest strain where catastrophic avalanche does not occur after that
int linAvalSizeDist[N]; // avalanche size distribution (linear binning)
int logAvalSizeDist[N]; // avalanche size distribution (log binning)
int catasAval[K];       // catastrophic avalanches
int maxAval[K];         // largest avalanche size (excluding catastrophic) in each bundle
int i, j, k;            // counters
// int n;                  // number of output points in an interval
int noOfAvalanches[K]; // number of avalanches

double rand01();

int CmpFunc(const void *_a, const void *_b);

void uniform(double array[], double mean, double radius, int sizeOfArray);
// initialize and sort arrays of uniform distribution

void loga_bin(int num_bins, int max, int min, int *array_pointer, FILE *file_pointer); // logarithmic binning

int main()
{
    double aveMaxAval = 0;        // average size of the largest avalanches
    double aveCatasAval = 0;      // average size of the catastrophic avalanches
    double aveNoOfAvalanches = 0; // average number of avalanches
    double aveSigmaC = 0;         // average largest stress where catastrophic avalanche does not occur after that
    double aveEpsC = 0;           // average largest strain where catastrophic avalanche does not occur after that
    double s_i, s_j, TotDamage, aveTotDamage = 0.0, aveTotAvalSize = 0.0, TotAvalSize;  // note: what are these new variables?
    int i_avalstart;
    char filename[50];     // name of the output file
    int avalSize, subAval; // avalanche size and sub-avalanche size

    // Initialize the arrays
    memset(linAvalSizeDist, 0, N * sizeof(int));
    memset(logAvalSizeDist, 0, N * sizeof(int));
    memset(catasAval, 0, K * sizeof(int));
    memset(maxAval, 0, K * sizeof(int));
    memset(noOfAvalanches, 0, K * sizeof(int));

    FILE *constit, *aveParamAlpha, *linAvalSizeDistFile, *logAvalSizeDistFile;
    int intpart = (int)alpha;
    int fracpart = (int)round((alpha - intpart) * 10000);
    snprintf(filename, sizeof(filename), "constit_%d_%dp%04d.txt", N, intpart, fracpart);
    constit = fopen(filename, "w");
    snprintf(filename, sizeof(filename), "aveParamAlpha_%d_%dp%04d.txt", N, intpart, fracpart);
    aveParamAlpha = fopen(filename, "w");
    snprintf(filename, sizeof(filename), "linAvalSizeDist_%d_%dp%04d.txt", N, intpart, fracpart);
    linAvalSizeDistFile = fopen(filename, "w");
    snprintf(filename, sizeof(filename), "logAvalSizeDist_%d_%dp%04d.txt", N, intpart, fracpart);
    logAvalSizeDistFile = fopen(filename, "w");

    srand(time(NULL));  // note: why is this needed?

    for (k = 0; k < K; k++)
    {
        // printf("%d\n", k);
        uniform(yielding, e_y, de_y, N);
        uniform(broken, e_b, de_b, N);

        // attention: Constitutive Relation
        double eps;   // strain of the whole bundle
        double force; // force on the whole bundle
        double sigma; // stress of the whole bundle

        // elastic-only regime ----------
        // since we don't need too many data points, we can reduce the number of points to faster the simulation
        for (i = 0; i < N / 25 - K; i++)
        {
            eps = yielding[0] * i / (N / 25 - K - 1);
            sigma = E * eps;
            if (K == 1)
                fprintf(constit, "%lf\t%lf\n", eps, sigma);
        }

        // yielding + elastic regime ----------------
        for (i = 0; i < N; i++)
        {
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
            if (K == 1 && i % 7 == 0)   // print out fewer data points
                fprintf(constit, "%lf\t%lf\n", eps, sigma);
        }

        // yielding-only regime -----------
        for (i = 0; i < N; i++)
        {
            eps = yielding[N - 1] + (broken[0] - yielding[N - 1]) * i / (N - 1);
            force = 0.0;
            for (j = 0; j < N; j++)
                force += E * yielding[j] + alpha * E * (eps - yielding[j]);
            sigma = force / N;
            if (K == 1 && i % 100 == 0) // print out fewer data points
                fprintf(constit, "%lf\t%lf\n", eps, sigma);
        }

        // yielding+broken regime ---------
        for (i = 0; i < N; i++)
        {
            force = 0.0;
            eps = broken[i];
            for (j = i; j < N; j++)
                force += E * yielding[j] + alpha * E * (eps - yielding[j]);
            sigma = force / N;
            s[i] = sigma;
            if (K == 1)
                fprintf(constit, "%lf\t%lf\n", eps, sigma);
        }

        // broken-only regime ------------
        eps = broken[N - 1];
        // since we don't need to many data points, we can reduce the number of points to faster the simulation
        for (i = 0; i < N / 100; i++)
        {
            sigma = 0.0;
            if (K == 1)
                fprintf(constit, "%lf\t%lf\n", eps, sigma);
            eps += (limit - broken[N - 1]) / (N / 100 - 1);
        }

        // New avalanche ------------------------

        TotDamage = 0.0;
        TotAvalSize = 0.0;
        i = 0;
        while ((i + 1) < N)
        {
            s_i = s[i];
            i_avalstart = i;
            j = 0;
            do
            {
                j++;
                s_j = s[i + j];
            } while (s_j < s_i && (i + j + 1) < N);
            avalSize = j;
            i = i + j;
            if ((i + 1) < N) // excluding catastrophic avalanches
            {
                TotDamage += avalSize;
                TotAvalSize += avalSize;

                linAvalSizeDist[avalSize]++;
                if (avalSize > maxAval[k])
                    maxAval[k] = avalSize;
                noOfAvalanches[k]++;
                sigmaC[k] = s[i_avalstart];
                epsC[k] = broken[i_avalstart];
            }
        }
        // This part is needed because if already the first avalanche is catastrophic the critical strain and stress would be automatically zero, which is not correct
        if (noOfAvalanches[k] == 0)
        {
            sigmaC[k] = s[i_avalstart];
            epsC[k] = broken[i_avalstart];
        }

        //----------------------------------------

        // // attention: Breaking Avalanche Size Distribution
        // i = -1;
        // while ((i + 1) < N)
        // {
        //     i++;
        //     avalSize = 1;
        //     do
        //     {
        //         subAval = 0;
        //         while (s[i + 1] < s[i] && (i + 1) < N)
        //         {
        //             i++;
        //             subAval++;
        //         }
        //         avalSize += subAval;
        //     } while (subAval > 0 && (i + 1) < N);
        //     if ((i + 1) < N) // excluding catastrophic avalanches
        //     {
        //         linAvalSizeDist[avalSize]++;
        //         if (avalSize > maxAval[k] && avalSize < N)
        //         {
        //             maxAval[k] = avalSize;
        //         }
        //         noOfAvalanches[k]++;
        //         sigmaC[k] = s[i];
        //         epsC[k] = broken[i];
        //     }
        // }
        // printf("sigmaC[%d]: %lf\n", k, sigmaC[k]);
        // printf("maxAval[%d]: %d\n", k, maxAval[k]);
        aveTotDamage += TotDamage;
        if (noOfAvalanches[k] > 0)
            aveTotAvalSize += TotAvalSize / noOfAvalanches[k];
        catasAval[k] = avalSize;
        // if(catasAval[k] == maxAval[k])
        //     maxAval[k] = 0;
    } // attention: end of K loop

    // attention: Parameters dependent on alpha
    i = 0; // to exclude maxAval = catasAval, i.e., catastrophic avalanche happens at the start
    for (k = 0; k < K; k++)
    {
        aveCatasAval += catasAval[k];
        aveNoOfAvalanches += noOfAvalanches[k];
        aveSigmaC += sigmaC[k];
        aveEpsC += epsC[k];
        if (maxAval[k] > 0)
        {
            aveMaxAval += maxAval[k];
            i++;
        }
    }
    aveTotDamage /= K;
    aveTotAvalSize /= K;
    aveCatasAval /= K;
    if (i == 0)
        aveMaxAval = 0;
    else
        aveMaxAval /= i;
    aveNoOfAvalanches /= K;
    aveSigmaC /= K;
    aveEpsC /= K;

    printf("Average size of catastrophic avalanches: %lf\n", aveCatasAval);
    printf("Average size of the largest avalanches: %lf\n", aveMaxAval);
    printf("Average number of avalanches: %lf\n", aveNoOfAvalanches);
    printf("Average largest stress where catastrophic avalanche does not occur after that: %lf\n", aveSigmaC);
    printf("Average largest strain where catastrophic avalanche does not occur after that: %lf\n", aveEpsC);

    fprintf(aveParamAlpha, "%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", alpha, aveCatasAval, aveMaxAval, aveNoOfAvalanches, aveSigmaC, aveEpsC, aveTotDamage, aveTotAvalSize);

    // attention: Output the linear and logarithmic avalanche size distributions
    for (i = 0; i < N; i++)
    {
        if (linAvalSizeDist[i] > 0)
            fprintf(linAvalSizeDistFile, "%d\t%d\n", i, linAvalSizeDist[i]);
    }
    loga_bin(binAval, N, 1, linAvalSizeDist, logAvalSizeDistFile);

    fclose(constit);
    fclose(linAvalSizeDistFile);
    fclose(logAvalSizeDistFile);
    return 0;
}

void loga_bin(int num_bins, int max, int min, int *array_pointer, FILE *file_pointer)
{
    double bin_size, x_i, x_ip1, x, y, dist_inc[num_bins];
    int i, j;

    memset(dist_inc, 0, sizeof(dist_inc));

    bin_size = (log(max) - log(min)) / (double)num_bins;

    for (i = min; i < max; i++)
    {
        j = floor((log(i) - log(min)) / bin_size);
        dist_inc[j] += array_pointer[i];
    }

    for (i = 0; i < num_bins; i++)
    {
        if (dist_inc[i] > 0)
        {
            x_i = min * exp(i * bin_size);
            x_ip1 = min * exp((i + 1) * bin_size);
            x = sqrt(x_i * x_ip1);
            y = dist_inc[i] / x;
            fprintf(file_pointer, "%lf\t%lf\n", x, y);
        }
    }
    // fclose(file_pointer);
    return;
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
