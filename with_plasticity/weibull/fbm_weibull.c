#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 10000
#define K 1
#define M_YIELD 5.0
#define L_YIELD 1.0
#define M_BREAK 15.0
#define L_BREAK 3.0
#define E 1.0
#define alpha 0.2
#define binAval 40 // number of bins for avalanche size distribution (log)

typedef struct
{
    unsigned int id;      // fiber ID: the i-th and (i+N)-th elements are the same fiber
    double threshold;     // threshold of the fiber (yielding or breaking): if id = i then threshold = yielding threshold, if id = i+N then threshold = breaking threshold
    unsigned int plastic; // if the fiber has e_y < e_b then plastic = 1, otherwise plastic = 0
    unsigned int broken;  // if the fiber is already broken then broken = 1, otherwise broken = 0
} FIBER;

double yielding[N];     // yielding thresholds of fibers
double broken[N];       // breaking thresholds of fibers
double s[N];            // stress of the whole bundle
double sigmaC[K];       // critical/largest stress where catastrophic avalanche does not occur after that
double epsC[K];         // critical/largest strain where catastrophic avalanche does not occur after that
int linAvalSizeDist[N]; // avalanche size distribution (linear binning)
int logAvalSizeDist[N]; // avalanche size distribution (log binning)
int catasAval[K];       // catastrophic avalanches
int maxAval[K];         // largest avalanche size (excluding catastrophic) in each bundle
int noOfAvalanches[K];  // number of avalanches
FIBER bundle[2 * N];    // bundle of fibers
int i, j, k;            // counters

double rand01();
int CmpFunc(const void *_a, const void *_b);
int CmpFiberFunc(const void *_a, const void *_b);
void uniform(double array[], double mean, double radius, int sizeOfArray);
// initialize and sort arrays of uniform distribution
double weibull(double m, double lambda);
// generate Weibull distribution
void loga_bin(int num_bins, int max, int min, int *array_pointer, FILE *file_pointer); // logarithmic binning
void init_bundle(FIBER array[]);                                                       // initialize the bundle of fibers

int main()
{
    double aveMaxAval = 0;        // average size of the largest avalanches
    double aveCatasAval = 0;      // average size of the catastrophic avalanches
    double aveNoOfAvalanches = 0; // average number of avalanches
    double aveSigmaC = 0;         // average largest stress where catastrophic avalanche does not occur after that
    double aveEpsC = 0;           // average largest strain where catastrophic avalanche does not occur after that
    double s_i, s_j, TotDamage, aveTotDamage = 0.0, aveTotAvalSize = 0.0, TotAvalSize;
    int i_avalstart;
    char filename[50];     // name of the output file
    int avalSize, subAval; // avalanche size and sub-avalanche size
    int i_sigma;           // counter for the s[N] array

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

    srand(time(NULL));

    for (k = 0; k < K; k++)
    {
        printf("%d\n", k);
        init_bundle(bundle);

        // attention: Constitutive Relation
        double eps;    // strain of the whole bundle
        double force;  // force on the whole bundle
        double sigma;  // stress of the whole bundle
        int intact;    // number of intact fibers
        int noOfYield; // number of yielding fibers

        i = 0;
        i_sigma = 0;
        while (i < (2 * N))
        {
            intact = N;
            noOfYield = 0;
            force = 0.0;

            while (bundle[i].broken == 1)
            {
                i++;
            }

            eps = bundle[i].threshold; // let strain jump to the smallest yielding/breaking threshold

            for (j = 0; j < i; j++)
            {
                if (bundle[j].broken == 0 && bundle[j].plastic == 1)
                {
                    force += E * bundle[j].threshold + alpha * E * (eps - bundle[j].threshold); // force on yielding fibers
                    noOfYield++;
                }
                else if (bundle[j].broken == 1)
                {
                    if (bundle[j].plastic == 1 && bundle[j].id < N)
                        intact -= 1;
                    // bundle[j].id < N is used to avoid double counting the broken fibers hence bundle[j].id >= N will do the same trick
                    else if (bundle[j].plastic == 0 && bundle[j].id >= N)
                        intact -= 1;
                }
            }
            force += E * eps * (intact - noOfYield); // force on fibers still elastic
            if(intact == 0)
            {
                sigma = 0;
                break;
            }
            else
                sigma = force / intact;
            s[i_sigma] = sigma;
            if (bundle[i].plastic == 1 && bundle[i].id >= N)
            {
                bundle[i].broken = 1;
                j = 0;
                while (j < i) // switch its pair element to be broken too
                {
                    if (bundle[i].id == (bundle[j].id + N)) // finding its pair
                    {
                        bundle[j].broken = 1;
                        break;
                    }
                    j++;
                }
            }
            else if (bundle[i].plastic == 0 && bundle[i].id >= N)
            {
                bundle[i].broken = 1;
                j = i + 1;
                while (j < (2 * N))
                {
                    if (bundle[i].id == (bundle[j].id + N)) // finding its pair
                    {
                        bundle[j].broken = 1;
                        break;
                    }
                    j++;
                }
            }
            i++;
            i_sigma++;
            if (K == 1)
                fprintf(constit, "%lf\t%lf\n", eps, sigma);
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

        aveTotDamage += TotDamage;
        if (noOfAvalanches[k] > 0)
            aveTotAvalSize += TotAvalSize / noOfAvalanches[k];
        catasAval[k] = avalSize;
        // if(catasAval[k] == maxAval[k])
        //     maxAval[k] = 0;
    } // attention: end of K loop

    /*     // attention: Parameters dependent on alpha
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
     */
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

double weibull(double m, double lambda)
{
    double x, r;
    r = rand01();
    x = pow(pow(lambda, m) * (-log(1.0 - r)), 1.0 / m);
    return x;
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

void init_bundle(FIBER fiber_array[])
{
    for (int i = 0; i < N; i++)
    {
        fiber_array[i].id = i;
        fiber_array[i].broken = 0; // at first, no fiber is broken
        fiber_array[i].threshold = weibull(M_YIELD, L_YIELD);
        fiber_array[i + N].threshold = weibull(M_BREAK, L_BREAK);
        if (fiber_array[i].threshold < fiber_array[i + N].threshold)
        {
            fiber_array[i].plastic = 1;
            fiber_array[i + N].plastic = 1;
        }
        else
        {
            fiber_array[i].plastic = 0;
            fiber_array[i + N].plastic = 0;
        }
    }
    qsort(fiber_array, 2 * N, sizeof(FIBER), CmpFiberFunc);
}

int CmpFiberFunc(const void *a, const void *b)
{
    FIBER *fiberA = (FIBER *)a;
    FIBER *fiberB = (FIBER *)b;

    if (fiberA->threshold > fiberB->threshold)
        return 1;
    else if (fiberA->threshold < fiberB->threshold)
        return -1;
    else
        return 0;
}