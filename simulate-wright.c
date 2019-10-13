// Simple mtDNA bottleneck calculations
// simulate-wright.c

// Simulates the process modelled by the Wright equation -- a series of binomial selections of (mtDNA) molecules from a population of size n
// Many (NSAMP) instances of this are run independently and the resulting mean and variance statistics are output

#include <stdio.h>
#include <stdlib.h>

#define NSAMP 10000
#define RND drand48()

int main(void)
{
  int n;
  double h[NSAMP];
  int i, j;
  int t;
  int count;
  double mean, var;
  FILE *fp;

  fp = fopen("simulate-wright.txt", "w");
  
  // loop through different bottleneck sizes
  for(n = 1; n < 100; n *= 2)
    {
      // initialise a set of samples
      for(i = 0; i < NSAMP; i++)
	{
	  h[i] = 0.5;
	}
      // loop through number of cell divisions
      for(t = 0; t < 100; t++)
	{
	  // loop through individual samples
	  for(i = 0; i < NSAMP; i++)
	    {
	      
	      // simulate a cell division that binomially chooses mtDNAs from a population of size n
	      count = 0;
	      for(j = 0; j < n; j++)
		{
		  if(RND < h[i]) count++;
		}
	      h[i] = (double) count/n;
	    }
	  // compute heteroplasmy statistics across samples
	  mean = var = 0;
	  for(i = 0; i < NSAMP; i++)
	    mean += h[i];
	  mean /= NSAMP;
	  for(i = 0; i < NSAMP; i++)
	    var += (mean-h[i])*(mean-h[i]);
	  var /= NSAMP-1.;

	  // output
	  fprintf(fp, "%i %i %f %f\n", n, t, mean, var);
	}
    }
  fclose(fp);
  
  return 0;
}
	
