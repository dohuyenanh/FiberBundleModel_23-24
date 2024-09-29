#include<stdio.h>
#include<math.h>
#include<string.h>	//for memset

#define	X_MAX	140		//the self-estimated value of x_max
#define	STEPS	60		//number of steps that x_min needs to increase

void alpha(FILE *indata, int x_max, int x_min, int steps, FILE *outdata);

int main()
{
	FILE *linaval, *alp_aval;
	linaval = fopen("../data/lin_aval_size.txt", "r");
	alp_aval = fopen("../data/alpha(x_min).txt", "w");

	alpha(linaval, X_MAX, 1, STEPS, alp_aval);

	fclose(linaval);
	fclose(alp_aval);

	return 0;
}

void alpha(FILE *indata, int x_max, int x_min, int steps, FILE *outdata)
// scaling parameter(read data, x_max, first (smallest) x_min, number of steps that x_min increases, file contains x_min-s and their corresponding alpha)
{
	// Print values for debugging
    printf("x_max: %d, x_min: %d, steps: %d\n", x_max, x_min, steps);

    if (x_max < x_min || steps <= 0) {
        fprintf(stderr, "Invalid parameters: x_max must be >= x_min and steps must be > 0\n");
        return;
    }

	int d[x_max - x_min + 1][2];								// array to store data which is in the range, first column is value and second column is its corresponding frequency
	int i, j;													// ask prof. why we should always(?) declare counters outside the for loop
	int dxmin = floor((double)(x_max - x_min) / (double)steps); // value of a increasing step
	// int n[steps+1];	//array stores number of values x>=x_min where (steps+1) is the number of x_min-s we will deal with
	int n;		// number of values x>=x_min
	double a;	// alpha
	double sum; // the sum inside the inverse function of the alpha formula
	int new_xmin;

	memset(d, 0, sizeof(d)); // initialize d to 0

	for (i = 0; i < (x_max - x_min + 1); i++)
	{
		if (i == 0)
			fscanf(indata, "%d%d\n", &d[i][0], &d[i][1]);
		else
		{
			if (d[i - 1][0] < x_max)
				fscanf(indata, "%d%d\n", &d[i][0], &d[i][1]);
			else
			{
				d[i][0] = 0;
				d[i][1] = 0;
			}
		}
	}

	for (i = 0; i < (steps + 1); i++) // loop for varying alpha-s
	{
		printf("%d\n", i);
		new_xmin = x_min + i * dxmin;
		n = 0;
		sum = 0.0;
		for (j = 0; j < (x_max - x_min + 1); j++)
		{
			if (d[j][1] >= new_xmin)
			{
				sum += d[j][1] * log((double)d[j][0] / (double)new_xmin);
				n++;
			}
		}
		sum = 1 / sum;
		a = 1 + n * sum;
		fprintf(outdata, "%d\t%.12lf\n", new_xmin, a);
	}

	return;
}
