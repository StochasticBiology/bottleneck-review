// Simple mtDNA bottleneck calculations
// bneck-errors.c

// Under one condition, the standard error associated with a variance measurement is V * sqrt(2 / (n-1)).
// This code assumes V = 0.01 and a mean heteroplasmy of 0.5 (hence h(1-h) = 0.25), and produces confidence intervals on an estimate of N = h(1-h)/V for different n. 

#include <stdio.h>
#include <math.h>

int main(void)
{
  int n;
  double vhp;
  double vh, hi, lo;
  FILE *fp;

  fp = fopen("bneck-errors.txt", "w");
  for(n = 2; n < 100; n++)
    {
      // mean, +sem and -sem estimates of V
      vh = 0.01;
      hi = 0.01+0.01*sqrt(2./(n-1));
      lo = 0.01-0.01*sqrt(2./(n-1));

      // avoid divide-by-zero problems (this is out of range on the resulting plot)
      if(lo <= 0) lo = 0.0001;

      // output values for h(1-h)/V
      fprintf(fp, "%i %f %f %f\n", n, 0.25/vh, 0.25/hi, 0.25/lo);
    }
  fclose(fp);

  return 0;
}
