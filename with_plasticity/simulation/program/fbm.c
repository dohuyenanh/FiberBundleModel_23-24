#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cjson/cJSON.h>

// Declare global variables for the parameters
double alpha_1, alpha_2, alpha_3, e_b, e_y, de_b_1, de_b_2, de_b_3, de_b_4, de_y, limit, E;
int N, no_of_alpha, no_of_de_b, no_of_de_y, no_of_e_b, no_of_e_y;
double *yielding; // yielding thresholds of fibers
double *broken;   // breaking thresholds of fibers
int i, j;         // counters

double rand01();

int CmpFunc(const void *_a, const void *_b);

void uniform(double array[], double mean, double radius, int sizeOfArray);
// initialize and sort arrays of uniform distribution

void constitutive(double e_yielding, double de_yielding, double e_breaking, double de_breaking, double alpha, FILE *outputFile);
// function to export output data of the constitutive relation as one of the parameter change

void read_config(const char *filename);
// Function to read the configuration file

int main()
{
    // Read the configuration file
    read_config("parameters.json");
    char filename[50];

    double alpha[3] = {alpha_1, alpha_2, alpha_3};
    double de_breaking[4] = {de_b_1, de_b_2, de_b_3, de_b_4};

    for (int index = 0; index < no_of_alpha; index++)
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

    // Free allocated memory
    free(yielding);
    free(broken);

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

void read_config(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
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
    if (!json) {
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
    de_b_4 = cJSON_GetObjectItem(values, "de_b_4")->valuedouble;
    de_y = cJSON_GetObjectItem(values, "de_y")->valuedouble;
    limit = cJSON_GetObjectItem(values, "limit")->valuedouble;
    N = cJSON_GetObjectItem(values, "N")->valueint;
    E = cJSON_GetObjectItem(values, "E")->valuedouble;

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

    // Generate arrays dynamically using uniform distribution
    uniform(yielding, e_y, de_y, N);  // generate yielding array
    uniform(broken, e_b, de_b_1, N);  // generate broken array using de_b_1 (you can switch to other de_b_x as needed)

    cJSON_Delete(json);
    free(data);
}