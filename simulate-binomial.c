// Simple mtDNA bottleneck calculations
// simulate-binomial.c

// Simulate a binomial-distribution model for mtDNA heteroplasmy

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NSAMP 10000
#define RND drand48()

// we don't use these factorial functions as it's simpler to simulate binomial events numerically
double fact(int n)
{
  if(n <= 0) return 1;
  return n*fact(n-1);
}

double binom(int n, int k)
{
  return fact(n) / (fact(k)*(fact(n-k)));
}

int main(void)
{
  double hist[100][101];
  int i, j;
  double x;
  double n;
  int score;
  int h;
  int k;
  FILE *fp;

  fp = fopen("simulate-binomial.txt", "w");
  j = -1;

  // the binomial distribution is simulated numerically
  // loop through bottleneck size
  for(n = 1; n < 1000; n *= 2)
    {
      j++;
      // loop through independent samples
      for(i = 0; i < NSAMP; i++)
	{
	  score = 0;
	  // do one sampling and record heteroplamy
	  for(k = 0; k < n; k++)
	    {
	      if(RND < 0.5) score++;
	    }
	  h = (100.*(double)score/n);
	  hist[j][(int)h]++;
	}
    }
  // output results of samples
  for(x = 0; x <= 1.01; x+= 0.01)
    {
      fprintf(fp, "%.3f ", x);
      for(i = 0; i < j; i++)
	{
	  fprintf(fp, "%.4f ", hist[i][(int)(x*100)]/NSAMP);
	}
      fprintf(fp, "\n");
    }
  fclose(fp);

  return 0;
}
      
