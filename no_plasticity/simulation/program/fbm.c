#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

//#define WEIBULL

#define N     100000
#define KMAX  1000
#define T     N*KMAX
#define M     40
#define MC    80
#define H     60    //number of bins for heights distribution of profiles
#define W1    10    //duration of profile
#define W2    12
#define W3    14    //difficult to use #define for an array
#define W4    16    //suppose W5 is always the largest number among W-s
#define W5    18    //if number of W-s changes, change size of l[], num_w[], ave_proflie[][]
//********************************************************

double strength[N], load_increment[N], dist_inc[M], dist_cataval[MC];
double m = 2.0, lambda = 1.0;
int avalanche_dist[N], duration_dist[N], duration, profile[N], h[T], wait_dist[N], aval_dur[2][N];

int CmpFunc(const void* _a, const void* _b);
double weibull(double ,double );
double uniform();
void initial_sample();
void loga_bin(int num_bins, int max, int min, int *array_pointer, FILE *file_pointer);

int main()
{
    int i, j, k, lav, lavloc, l[5], avalanche, avalloc, cataval_min, cataval_max, cataval[KMAX], ave_profile[5][W5], num_w[5], h_min, h_max, p, p_max, dist_hei[H], max_dur = 0, t, max_avaldist = 1;
    double sigma, s, eps, wait;
    double dm, max_inc, x_i, x_ip1, x, y, dh;

    int nbroken, nintact;
    FILE *constit, *inc, *liaval, *loaval, *cat, *cat_dist, *lidur, *pro5, *lihei, *w8t, *lodur, *lohei, *avadur;
    constit = fopen("../data/constit.txt", "w");
    inc = fopen("../data/load_increment.txt", "w");
    liaval = fopen("../data/lin_aval_size.txt", "w");
    loaval = fopen("../data/log_aval_size.txt", "w");
    // cat = fopen("../data/cat_aval.txt", "w");
    cat_dist = fopen("../data/dist_cat_aval.txt", "w");
    lidur = fopen("../data/lin_dur_dist.txt", "w");
    pro5 = fopen("../data/average_profile.txt", "w");
    lihei = fopen("../data/lin_hei_dist.txt", "w");
    w8t = fopen("../data/wait_time_dist.txt", "w");
    lodur = fopen("../data/log_dur_dist.txt", "w");
    lohei = fopen("../data/log_hei_dist.txt", "w");
    avadur = fopen("../data/dur_vs_avg_aval.txt", "w");

    memset(avalanche_dist, 0, sizeof(avalanche_dist));
    memset(duration_dist, 0, sizeof(duration_dist));
    memset(profile, 0, sizeof(profile));
    memset(l, 0, sizeof(l));
    memset(ave_profile, 0, sizeof(ave_profile));
    memset(cataval, 0, sizeof(cataval));
    memset(num_w, 0, sizeof(num_w));
    memset(dist_hei, 0, sizeof(dist_hei));
    memset(h, 0, sizeof(h)); // array to store profiles' heights
    memset(wait_dist, 0, sizeof(wait_dist));
    memset(aval_dur, 0, sizeof(aval_dur));
    memset(wait_dist, 0, sizeof(wait_dist));

    for (i = 0; i < 2; i++)
        for (j = 0; j < N; j++)
            aval_dur[i][j] = 0; // 1st row counts the number of avalanches of duration w and 2nd row is the sum of all these avalanches

    p = 0; // counter for profiles

    for (k = 0; k < KMAX; k++)
    {
        printf("%d\n", k);
        // intial conditions:
        initial_sample();
        nbroken = 0;
        nintact = N;
        i = -1;
        sigma = 0.0;
        // end of initial conditions

        while (i + 1 < N)
        {
            i++; // first round: i=0
            if (i != 0)
            { 
                wait = strength[i] - sigma; // strength is always between 0 and 1
                wait = floor(wait * N);     // to distinguish between different wait time
                wait_dist[(int)wait]++;     // frequency of a waiting time
            }
            sigma = strength[i]; // first round: strength[0] (the smallest)
            avalanche = 1;
            duration = 0;
            s = sigma * (double)nintact / (double)N;
            fprintf(constit, "%lf\t%lf\n", sigma, s);
            profile[duration] = 1;

            do
            {
                sigma = sigma * (double)nintact / (double)(N - i - 1);
                nbroken = i + 1;
                nintact = N - nbroken;
                avalloc = 0;

                while (sigma > strength[i + 1] && i + 1 < N)
                {
                    i++;
                    avalloc++; // sub-avalanche
                }

                avalanche += avalloc;
                duration++;
                profile[duration] = avalloc;

            } while (avalloc > 0 && i + 1 < N);

            if ((i + 1) < N) // excluding catastrophic avalanches
            {
                avalanche_dist[avalanche]++;                         // frequency of avalanche of a specific number of fibers
                duration_dist[duration]++;                           // frequency of specific duration
                max_dur = (max_dur > duration) ? max_dur : duration; // to find the maximum duration
                aval_dur[0][duration]++;                             // 1st row: same as duration_dist[]
                aval_dur[1][duration] += avalanche;

                // finding height of a profile
                for (j = 0; j < duration; j++)
                {
                    if (h[p] < profile[j])
                        h[p] = profile[j];
                }
                p++;
                // done

                int durations[] = {W1, W2, W3, W4, W5};

                for (int i = 0; i < sizeof(durations) / sizeof(durations[0]); i++)
                {
                    if (duration == durations[i])
                    {
                        num_w[i]++;
                        for (j = 0; j < durations[i]; j++)
                        {
                            ave_profile[i][j] += profile[j];
                        }
                        break;
                    }
                }
            } // done with excluding catastrophic avalanches
        } // while(i+1 < N) ------------------------

        cataval[k] = avalanche;
        //      fprintf(cat, "%d\t%d\n", k, cataval[k]);
    } // for (k=0; k<KMAX; k++) ------------------------------------
    fclose(constit);
    // fclose(cat);

    for (i = 0; i <= max_dur; i++)
        if (aval_dur[0][i] > 0)
            fprintf(avadur, "%5d\t%lf\n", i, (double)aval_dur[1][i] / (double)aval_dur[0][i]);
    fclose(avadur);

    // waiting time distribution
    for (j = 0; j < N; j++)
    {
        if (wait_dist[j] != 0)
            fprintf(w8t, "%d\t%d\n", j, wait_dist[j]);
    }
    fclose(w8t);
    // done

    p_max = p;
    // finding the biggest and smallest height
    h_max = h[0];
    h_min = h[0];
    for (p = 0; p < p_max; p++)
    {
        if (h_max < h[p])
            h_max = h[p];
        if (h_min > h[p])
            h_min = h[p];
    }
    // finding the size of a bin
    dh = (h_max - h_min) / H;
    // put h[p] in the right bin
    for (p = 0; p < p_max; p++)
    {
        j = floor((double)(h[p] - h_min) / dh);
        dist_hei[j]++; // bin
    }
    // printf("%d %d", h_min, h_max);
    for (j = 0; j < H; j++)
        if (dist_hei[j] > 0)
            fprintf(lihei, "%lf\t%20.15lf\n", h_min + j * dh, (double)dist_hei[j] / ((double)p_max * dh));
    fclose(lihei);

    loga_bin(50, H, 1, dist_hei, lohei);

    for (j = 0; j < W5; j++)
    {
        for (i = 4; i >= 0; i--)
        {

            if (ave_profile[i][j] != 0)
            {
                fprintf(pro5, "%d\t%7.6lf\t", j, (double)ave_profile[i][j] / (double)num_w[i]);
                if (l[i] == 0)
                    if (ave_profile[i][j + 1] == 0)
                        l[i] = j;
            }

            if (ave_profile[i][j] == 0)
            {
                ave_profile[i][j] = ave_profile[i][j - 1];
                fprintf(pro5, "%d\t%7.6lf\t", l[i], (double)ave_profile[i][j] / (double)num_w[i]);
            }

            if (i == 0)
            {
                fprintf(pro5, "\n");
            }
        }
    }
    fclose(pro5);

    for (i = 0; i < N; i++)
    {
        if (avalanche_dist[i] > 0)
        {
            fprintf(liaval, "%d\t%d\n", i, avalanche_dist[i]);
            if (max_avaldist < avalanche_dist[i])
                max_avaldist = avalanche_dist[i];
        }
        if (duration_dist[i] > 0)
            fprintf(lidur, "%d\t%d\n", i, duration_dist[i]);
    }
    fclose(liaval);
    fclose(lidur);

    loga_bin(60, max_dur, 1, duration_dist, lodur);

    /* // logarithmic binning for distribution of avalanche' size
    dm = (log(N) - log(1)) / M;     // bin size
    for(i=1; i<N; i++)
    {
        j = floor((log(i) - log(1)) / dm);
        // log(0) is undefined & avalanche >= 1
        dist_inc[j] += avalanche_dist[i];
    }

    for(i=0; i<M; i++)
    {
        if(dist_inc[i]>0)
        {
            x_i = 1 * exp(i*dm);
            x_ip1 = 1 * exp((i+1)*dm);
            x = sqrt(x_i * x_ip1);
            y = dist_inc[i] / x;
            fprintf(loaval, "%lf\t%lf\n", x, y);
        }
    }
    fclose(loaval);
    //done */

    loga_bin(M, N, 1, avalanche_dist, loaval);

    cataval_min = cataval[0];
    for (k = 0; k < KMAX; k++)
    {
        if (cataval[k] <= cataval_min)
            cataval_min = cataval[k];
    }
    cataval_max = cataval[0];
    for (k = 0; k < KMAX; k++)
    {
        if (cataval[k] >= cataval_min)
            cataval_max = cataval[k];
    }

    dm = (cataval_max - cataval_min) / (double)MC;

    for (i = 0; i < MC; i++)
        dist_cataval[i] = 0;

    for (i = 0; i < KMAX; i++)
    {
        j = floor((cataval[i] - cataval_min) / dm);
        dist_cataval[j]++;
    }

    for (i = 0; i < MC; i++)
    {
        dist_cataval[i] /= (double)KMAX * dm;
        fprintf(cat_dist, "%lf\t%lf\n", cataval_min + i * dm, dist_cataval[i]);
    }
    fclose(cat_dist);

    return 0;
}
//*************************************************************************

/*logarithmic binning function: avalsize and duration*/
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

double uniform()
{
    return (double)rand() / (double)RAND_MAX;
}

double weibull(double m, double lambda)
{
    double x, r;
    r = (double)rand() / (double)RAND_MAX;
    x = pow(pow(lambda, m) * (-log(1.0 - r)), 1.0 / m);

    return x;
}

void initial_sample()
{

    int i, k;
    for (i = 0; i < N; i++)
    {
#ifdef WEIBULL
        strength[i] = weibull(m, lambda);
#else
        strength[i] = uniform();
#endif
    }

    //  for (i=0; i<MAX; i++)
    //    printf("\n %f ", strength[i]);
    //  printf("\n ********************* \n");

    qsort(strength, N, sizeof(strength[0]), CmpFunc);

    /*  for (i=0; i<MAX; i++)
        printf("\n %f ", strength[i]);*/

    return;
}
