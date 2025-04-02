#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 10000
#define K 50
#define E 1.0
#define alpha 0.8
#define e_y 0.45
#define de_y 0.1
#define e_b 5.9
#define de_b 4.0
#define limit 11
#define binAval 40     // number of bins for avalanche size distribution (log)
#define binEnergy 5000 // number of bins for dissipated energy distribution (lin)
#define logBinEnergy 40 // number of bins for dissipated energy distribution (log)
#define dE 0.025       // bin size of energy

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
int noOfAvalanches[K];           // number of avalanches
double E_max;                    // maximum energy (safe guess)
double avalEnergy[N];            // array storing dissipated energy in each avalanche
double linEnergyDist[binEnergy]; // dissipated energy distribution (lin binning)
double logEnergyDist[binEnergy]; // dissipated energy distribution (log binning)

double rand01();

int CmpFunc(const void *_a, const void *_b);

void uniform(double array[], double mean, double radius, int sizeOfArray);
// initialize and sort arrays of uniform distribution

void loga_bin_int(int num_bins, int max, int min, int *array_pointer, FILE *file_pointer); // logarithmic binning for array of integers
void loga_bin_double(int logNoOfBins, int maxValue, int minValue, double *linDistArray, FILE *outputFile); // logarithmic binning for array of doubles


int main()
{
    double aveMaxAval = 0;        // average size of the largest avalanches
    double aveCatasAval = 0;      // average size of the catastrophic avalanches
    double aveNoOfAvalanches = 0; // average number of avalanches
    double aveSigmaC = 0;         // average largest stress where catastrophic avalanche does not occur after that
    double aveEpsC = 0;           // average largest strain where catastrophic avalanche does not occur after that
    double s_i, s_j, totDamage, aveTotDamage = 0.0, aveTotAvalSize = 0.0, totAvalSize;
    int i_avalstart;
    char filename[50];     // name of the output file
    int avalSize, subAval; // avalanche size and sub-avalanche size
    double dlog=0.2, logmin=log(0.000001); // bin size of the log bin histogram, log of the smallest energy

    // Initialize the arrays
    memset(linAvalSizeDist, 0, N * sizeof(int));
    memset(logAvalSizeDist, 0, N * sizeof(int));
    memset(catasAval, 0, K * sizeof(int));
    memset(maxAval, 0, K * sizeof(int));
    memset(noOfAvalanches, 0, K * sizeof(int));
    memset(linEnergyDist, 0, binEnergy * sizeof(double));
    memset(logEnergyDist, 0, binEnergy * sizeof(double));
    memset(avalEnergy, 0, N * sizeof(double));

    FILE *constit, *aveParamAlpha, *linAvalSizeDistFile, *logAvalSizeDistFile, *linEnergyDistFile, *logEnergyDistFile;
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
    snprintf(filename, sizeof(filename), "linEnergyDist_%d_%dp%04d.txt", N, intpart, fracpart);
    linEnergyDistFile = fopen(filename, "w");
    snprintf(filename, sizeof(filename), "logEnergyDist_%d_%dp%04d.txt", N, intpart, fracpart);
    logEnergyDistFile = fopen(filename, "w");

    srand(time(NULL));

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
            {
                force += E * yielding[j] + alpha * E * (eps - yielding[j]);
            }
            sigma = force / N;
            s[i] = sigma; // stress-controlled loading
            if (K == 1)
                fprintf(constit, "%lf\t%lf\n", eps, sigma);
        }

        double min_sigma = s[0];
        for(i = 0; i < N; i++)
        {
            if (s[i] < min_sigma)
                min_sigma = s[i];
        }
        printf("Minimum stress: %lf\n", min_sigma);

        // broken-only regime ------------
        eps = broken[N - 1];
        // since we don't need too many data points, we can reduce the number of points to faster the simulation
        for (i = 0; i < N / 100; i++)
        {
            sigma = 0.0;
            if (K == 1)
                fprintf(constit, "%lf\t%lf\n", eps, sigma);
            eps += (limit - broken[N - 1]) / (N / 100 - 1);
        }

        // New avalanche ------------------------

        totDamage = 0.0;
        totAvalSize = 0.0;
        i = 0;
        int avalIndex = 0; // index for avalEnergy[]
        double energy;     // dissipated energy in avalanche
        while ((i + 1) < N)
        {
            energy = 0.0;
            s_i = s[i];
            i_avalstart = i;
            j = 0;
            do
            {
                energy += 0.5 * alpha * E * s[i + j] * s[i + j]; // a fiber's energy
                j++;
                s_j = s[i + j]; // stress-controled loading
            } while (s_j < s_i && (i + j + 1) < N);
            avalSize = j;
            i = i + j;
            if ((i + 1) < N) // excluding catastrophic avalanches
            {
                totDamage += avalSize;
                totAvalSize += avalSize;

                linAvalSizeDist[avalSize]++;
                if (avalSize > maxAval[k])
                    maxAval[k] = avalSize;
                noOfAvalanches[k]++;
                sigmaC[k] = s[i_avalstart];
                epsC[k] = broken[i_avalstart];

                avalEnergy[avalIndex] = energy;
                avalIndex++;
            }
        } // attention: end of avalanche loop
        aveTotDamage += totDamage;

        // This part is needed because if already the first avalanche is catastrophic the critical strain and stress would be automatically zero, which is not correct
        if (noOfAvalanches[k] == 0)
        {
            sigmaC[k] = s[i_avalstart];
            epsC[k] = broken[i_avalstart];
        }

        if (noOfAvalanches[k] > 0)
            aveTotAvalSize += totAvalSize / noOfAvalanches[k]; // the right hand side is the average avalanche size of each bundle

        catasAval[k] = avalSize; // catastrophic avalanches size
        // if(catasAval[k] == maxAval[k])
        //     maxAval[k] = 0;

        /*  Finding the safe guess E_max:
        E_max = avalEnergy[0];
        for (int ii = 0; ii < avalIndex; ii++)
            E_max = (avalEnergy[ii] > E_max) ? avalEnergy[ii] : E_max;
        printf("%lf\n", E_max);
         */

        for (int ii = 0; ii < avalIndex; ii++)
        {
            int jj = floor(avalEnergy[ii] / dE); // note: do we not count energy that is smaller than dE?
            jj = (jj > binEnergy-1) ? binEnergy-1 : jj;
            linEnergyDist[jj]++;
            jj = floor((log(avalEnergy[ii])-logmin) / dlog); // note: do we not count energy that is smaller than dE?
            jj = (jj > binEnergy-1) ? binEnergy-1 : jj;
            logEnergyDist[jj]++;
        }

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
    aveSigmaC /= K;
    aveEpsC /= K;

    for (i = 0; i < binEnergy; i++)
    {
        if (linEnergyDist[i] > 0)
        {
            linEnergyDist[i] /= aveNoOfAvalanches;
            fprintf(linEnergyDistFile, "%lf %.10lf\n", (i + 0.5) * dE, linEnergyDist[i]);
        }
        if(logEnergyDist[i] > 0)
        {
            double xx1, xx2, xx;
            xx1 = exp(logmin) * exp(((double)i) * dlog);
            xx2 = exp(logmin) * exp(((double)(i+1)) * dlog);
            xx = sqrt(xx1*xx2);
            fprintf(logEnergyDistFile,"%lf %.10lf\n", xx, logEnergyDist[i]/(xx*aveNoOfAvalanches));
        }
    }
     // loga_bin_double(logBinEnergy, binEnergy, 1, linEnergyDist, logEnergyDistFile);

    aveNoOfAvalanches /= K;

    printf("Average size of catastrophic avalanches: %lf\n", aveCatasAval);
    printf("Average size of the largest avalanches: %lf\n", aveMaxAval);
    printf("Average number of avalanches: %lf\n", aveNoOfAvalanches);
    printf("Average largest stress where catastrophic avalanche does not occur after that: %lf\n", aveSigmaC);
    printf("Average largest strain where catastrophic avalanche does not occur after that: %lf\n", aveEpsC);
    printf("Average total damage (number of fibers fail during avalanches, i.e. fiber that break with only itself and does not cause avalanche is not counted as damage): %lf\n", aveTotDamage);
    printf("Average of 'average avalanche size of each bundle': %lf\n", aveTotAvalSize);

    fprintf(aveParamAlpha, "%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", alpha, aveCatasAval, aveMaxAval, aveNoOfAvalanches, aveSigmaC, aveEpsC, aveTotDamage, aveTotAvalSize);

    // attention: Output the linear and logarithmic avalanche size distributions
    for (i = 0; i < N; i++)
    {
        if (linAvalSizeDist[i] > 0)
            fprintf(linAvalSizeDistFile, "%d\t%d\n", i, linAvalSizeDist[i]);
    }
    loga_bin_int(binAval, N, 1, linAvalSizeDist, logAvalSizeDistFile);

    // attention: Output the logarithmic dissipated energy distribution

    fclose(constit);
    fclose(linAvalSizeDistFile);
    fclose(logAvalSizeDistFile);
    fclose(linEnergyDistFile);
    fclose(logEnergyDistFile);
    fclose(aveParamAlpha);
    return 0;
}

void loga_bin_int(int num_bins, int max, int min, int *array_pointer, FILE *file_pointer)
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

void loga_bin_double(int logNoOfBins, int maxValue, int minValue, double *linDistArray, FILE *outputFile)
{
    double logBinWidth, lowerBoundary, upperBoundary, binCenter, density;
    double logBinnedValues[logNoOfBins];
    int index, logBinIndex;

    memset(logBinnedValues, 0, sizeof(logBinnedValues));

    logBinWidth = (log(maxValue) - log(minValue)) / (double)logNoOfBins;

    for (index = minValue; index < maxValue; index++)
    {
        if (linDistArray[index] > 0)
        {
            logBinIndex = floor((log(index) - log(minValue)) / logBinWidth);
            logBinnedValues[logBinIndex] += linDistArray[index];
        }
    }

    for (index = 0; index < logNoOfBins; index++)
    {
        if (logBinnedValues[index] > 0)
        {
            lowerBoundary = minValue * exp(index * logBinWidth);
            upperBoundary = minValue * exp((index + 1) * logBinWidth);
            binCenter = sqrt(lowerBoundary * upperBoundary);
            density = logBinnedValues[index] / binCenter;   // normalization
            fprintf(outputFile, "%lf\t%.10lf\n", binCenter, density);
        }
    }
    return;
}