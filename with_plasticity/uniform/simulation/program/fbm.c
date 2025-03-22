/* To compile:
 clang -o fbm fbm.c -I/opt/homebrew/Cellar/cjson/1.7.18/include -L/opt/homebrew/opt/cjson/lib -lcjson
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <cjson/cJSON.h>

// Declare global variables for the parameters
double alpha_1 = 0.2, alpha_2 = 0.5, alpha_3 = 0.8;
double e_b = 5.9, de_b_1 = 4.5, de_b_2 = 4.0, de_b_3 = 3.5 /* , de_b_4 */; 
double limit = 9, E = 1.0, e_y = 0.45, de_y = 0.1;
int K = 1000, N = 100000 /*, no_of_alpha, no_of_de_b, no_of_de_y, no_of_e_b, no_of_e_y */;
int *maxAval; // largest avalanche size (excluding catastrophic) in each bundle and their average
double aveMaxAval;
int binAval = 40;          // number of bins for avalanche size distribution (log)
double *yielding;     // yielding thresholds of fibers
double *broken;       // breaking thresholds of fibers
double *s;            // stress of the whole bundle
int *linAvalSizeDist; // avalanche size distribution (linear binning)
int *logAvalSizeDist; // avalanche size distribution (log binning)
int *catasAval;       // catastrophic avalanches and their average
double aveCatasAval;
int i, j, k; // counters
int n;       // number of non-zero elements in arrays

double rand01();

int CmpFunc(const void *_a, const void *_b);

void uniform(double array[], double mean, double radius, int sizeOfArray);
// initialize and sort arrays of uniform distribution

void constitutive(double e_yielding, double de_yielding, double e_breaking, double de_breaking, double alpha, FILE *outputFile);
// function to export output data of the constitutive relation as one of the parameter change

//void read_config(const char *filename);
// Function to read the configuration file

void loga_bin(int num_bins, int max, int min, int *array_pointer, FILE *file_pointer);
// logarithmic binning

int main()
{
    // Read the configuration file
    //read_config("parameters.json");

    // Initialize yielding and broken arrays dynamically based on the parameters
    yielding = malloc(N * sizeof(double));
    broken = malloc(N * sizeof(double));
    s = malloc(N * sizeof(double));
    linAvalSizeDist = malloc(N * sizeof(int));
    logAvalSizeDist = malloc(N * sizeof(int));
    catasAval = malloc(K * sizeof(int));
    maxAval = malloc(K * sizeof(int));

    // Initialize the arrays
    memset(linAvalSizeDist, 0, N * sizeof(int));
    memset(logAvalSizeDist, 0, N * sizeof(int));
    memset(catasAval, 0, K * sizeof(int));
    memset(maxAval, 0, K * sizeof(int));

    char filename[50];
    double alpha[3] = {alpha_1, alpha_2, alpha_3}; // if not in the varying alpha case, always use alpha[0] for a better looking plot
    double de_breaking[3] = {de_b_1, de_b_2, de_b_3 /* , de_b_4 */};

    /*
        //------------- Constitutive relation for varying parameters --------------

        // Varying alpha
        for (int index = 0; index < no_of_alpha; index++)
        {
            // intialize yielding[] and broken[]
            uniform(yielding, e_y, de_y, N);
            uniform(broken, e_b, de_b_2, N);

            // name for output file
            sprintf(filename, "../data/constit(alpha-%d).txt", index + 1);
            FILE *constit;
            constit = fopen(filename, "w");

            constitutive(e_y, de_y, e_b, de_breaking[2], alpha[index], constit);
        }

        // Varying de_b
        for (int index = 0; index < no_of_de_b; index++)
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

     */
    // attention:
    /*---------------- Breaking Avalanche Size Distribution -------------------
    -> take into account of the yielding+broken regime */
/*     for (k = 0; k < K; k++)
    {
        printf("Bundle no.%d\n", k + 1);
        // intialize yielding[] and broken[]
        uniform(yielding, e_y, de_y, N);
        uniform(broken, e_b, de_b_2, N);

        for (i = 0; i < N; i++)
        {
            double sigma = 0.0;
            double force = 0.0;
            double eps = broken[i];
            for (j = i; j < N; j++)
                force += E * yielding[j] + (alpha[0]) * E * (eps - yielding[j]);
            sigma = force / N;
            s[i] = sigma;
        }

        int avalSize, subAval;
        i = -1;
        while ((i + 1) < N)
        {
            i++;
            avalSize = 1;
            do
            {
                subAval = 0;
                while (s[i + 1] < s[i] && (i + 1) < N)
                {
                    i++;
                    subAval++;
                }
                avalSize += subAval;
            } while (subAval > 0 && (i + 1) < N);
            if ((i + 1) < N) // excluding catastrophic avalanches
            {
                linAvalSizeDist[avalSize]++;
                if (avalSize > maxAval[k])
                {
                    maxAval[k] = avalSize;
                }
            }
        }
        printf("maxAval[%d]: %d\n", k, maxAval[k]);
        catasAval[k] = avalSize;
    }
    aveCatasAval = 0;
    for (k = 0; k < K; k++)
    {
        aveCatasAval += catasAval[k];
    }
    aveCatasAval /= K;
    printf("Average size of catastrophic avalanches: %lf\n", aveCatasAval);

    aveMaxAval = 0;
    for (k = 0; k < K; k++)
    {
        aveMaxAval += maxAval[k];
    }
    aveMaxAval /= K;
    printf("Average size of the largest avalanches: %lf\n", aveMaxAval);
    FILE *aveParamAlpha;
    sprintf(filename, "../data/aveParamAlpha_%.2f.txt", alpha[0]);
    aveParamAlpha = fopen(filename, "a");
    fprintf(aveParamAlpha, "%f\t%lf\t%lf\n", alpha[0], aveCatasAval, aveMaxAval);

    FILE *linAvalSizeDistFile;
    sprintf(filename, "../data/linAvalSizeDist_%.2f.txt", alpha[0]);
    linAvalSizeDistFile = fopen(filename, "w");
    for (i = 0; i < N; i++)
    {
        if (linAvalSizeDist[i] > 0)
            fprintf(linAvalSizeDistFile, "%d\t%d\n", i, linAvalSizeDist[i]);
    }
    fclose(linAvalSizeDistFile);

    FILE *logAvalSizeDistFile;
    sprintf(filename, "../data/logAvalSizeDist_%.2f.txt", alpha[0]);
    logAvalSizeDistFile = fopen(filename, "w");
    loga_bin(binAval, N, 1, linAvalSizeDist, logAvalSizeDistFile);
 */

    FILE *constit;
    constit = fopen("../data/constit.txt", "w");
    uniform(yielding, e_y, de_y, N);
    uniform(broken, e_b, de_b_2, N);
    constitutive(e_y, de_y, e_b, de_b_2, alpha[0], constit);

    // Free allocated memory
    free(yielding);
    free(broken);
    free(s);
    free(linAvalSizeDist);
    free(logAvalSizeDist);
    free(catasAval);
    free(maxAval);

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
// attention:
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
    fclose(file_pointer);
    return;
}
// attention:
/* 
void read_config(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);
    char *data = malloc(length + 1);
    fread(data, 1, length, file);
    fclose(file);

    cJSON *json = cJSON_Parse(data);
    if (!json)
    {
        fprintf(stderr, "Error parsing JSON\n");
        exit(EXIT_FAILURE);
    }

    // Read parameters from the first object
    cJSON *values = cJSON_GetObjectItem(json, "values");
    alpha_1 = cJSON_GetObjectItem(values, "alpha_1")->valuedouble;
    alpha_2 = cJSON_GetObjectItem(values, "alpha_2")->valuedouble;
    alpha_3 = cJSON_GetObjectItem(values, "alpha_3")->valuedouble;
    e_b = cJSON_GetObjectItem(values, "e_b")->valuedouble;
    e_y = cJSON_GetObjectItem(values, "e_y")->valuedouble;
    de_b_1 = cJSON_GetObjectItem(values, "de_b_1")->valuedouble;
    de_b_2 = cJSON_GetObjectItem(values, "de_b_2")->valuedouble;
    de_b_3 = cJSON_GetObjectItem(values, "de_b_3")->valuedouble;
    // de_b_4 = cJSON_GetObjectItem(values, "de_b_4")->valuedouble;
    de_y = cJSON_GetObjectItem(values, "de_y")->valuedouble;
    limit = cJSON_GetObjectItem(values, "limit")->valuedouble;
    N = cJSON_GetObjectItem(values, "N")->valueint;
    E = cJSON_GetObjectItem(values, "E")->valuedouble;
    binAval = cJSON_GetObjectItem(values, "binAval")->valueint;
    K = cJSON_GetObjectItem(values, "K")->valueint;

    // Read parameters from the second object
    cJSON *number_of_values = cJSON_GetObjectItem(json, "number of values");
    no_of_e_y = cJSON_GetObjectItem(number_of_values, "no. of e_y")->valueint;
    no_of_e_b = cJSON_GetObjectItem(number_of_values, "no. of e_b")->valueint;
    no_of_de_y = cJSON_GetObjectItem(number_of_values, "no. of de_y")->valueint;
    no_of_de_b = cJSON_GetObjectItem(number_of_values, "no. of de_b")->valueint;
    no_of_alpha = cJSON_GetObjectItem(number_of_values, "no. of alpha")->valueint;

    // Initialize yielding and broken arrays dynamically based on the parameters
    yielding = malloc(N * sizeof(double));
    broken = malloc(N * sizeof(double));
    s = malloc(N * sizeof(double));
    linAvalSizeDist = malloc(N * sizeof(int));
    logAvalSizeDist = malloc(N * sizeof(int));
    catasAval = malloc(K * sizeof(int));
    maxAval = malloc(K * sizeof(int));

    // Initialize the arrays
    memset(linAvalSizeDist, 0, N * sizeof(int));
    memset(logAvalSizeDist, 0, N * sizeof(int));
    memset(catasAval, 0, K * sizeof(int));
    memset(maxAval, 0, K * sizeof(int));

    // Generate arrays dynamically using uniform distribution
    // uniform(yielding, e_y, de_y, N); // generate yielding array
    // uniform(broken, e_b, de_b_2, N); // generate broken array using de_b_2 (switch to other de_b_x as needed) 

    cJSON_Delete(json);
    free(data);
}
*/