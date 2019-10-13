// Simple mtDNA bottleneck calculations
// simulate-kimura.c

// Simulate a Kimura-distribution model for mtDNA heteroplasmy
// This uses the recursion relation from Hoang-Binh, D. (2005). A program to compute exact hydrogenic radial integrals oscillator strengths, and Einstein coefficients, for principal quantum numbers up to n approximate to 1000. Comput. Phys. Commun. 166, 191–196, which is cited in Wonnapinij et al. Am J Hum Genet 83 582 (2008)
// This could be sped up with dynamic programming but as the code only takes a minute anyway this is not implemented

#include <stdio.h>
#include <math.h>

#define NTERM 20

// recursive function used to compute Kimura pdf
double F(int a, double b, double c, double z, int orig)
{
  if(a == -1) return 1.-(b*z/c);
  if(a == 0) return 1;

  return ( (a+1)*(1-z)*F(a+1, b, c, z, orig) - (a+1)*(1-z)*F(a+2, b, c, z, orig) + (a+1 + b*z - c)*F(a+1, b, c, z, orig) ) / (a+1-c);
}

// probability mass at 0
double f0(double p0, double b)
{
  int i;
  double sum = 0;
  double delta;

  for(i = 1; i < NTERM; i++)
    {
      delta = (2*i+1)*p0*(1.-p0)*F(1-i, i+2, 2, 1-p0, 1-i)*pow(b, i*(i+1)/2.);    
	if(i % 2 == 1) delta *= -1; 
	sum += delta;
    }

  return (1.-p0) + sum;
}

// probability mass at 1
double f1(double p0, double b)
{
  int i;
  double sum = 0;
  double delta;

  for(i = 1; i < NTERM; i++)
    {
      delta = (2*i+1)*p0*(1.-p0)*F(1-i, i+2, 2, p0, 1-i)*pow(b, i*(i+1)/2.);    
	if(i % 2 == 1) delta *= -1; 
	sum += delta;
    }

  return p0 + sum;
}

// probability mass at x between 0 and 1 exclusive
double phi(double p0, double b, double x)
{
  int i;
  double sum = 0;
  double delta;

  if(fabs(x) < 0.001) return f0(p0, b);
  if(fabs(x-1) < 0.001) return f1(p0, b);
  for(i = 1; i < NTERM; i++)
    {
      delta = i*(i+1)*(2*i+1)*p0*(1-p0)*F(1-i, i+2, 2, x, 1-i)*F(1-i, i+2, 2, p0, 1-i)*pow(b, i*(i+1)/2);
      sum += delta;
    }

  return sum;
}

int main(void)
{
  double p0, b;
  double x;
  double r, total;
  double n;
  FILE *fp;

  fp = fopen("simulate-kimura.txt", "w");

  // loop through heteroplasmu
  for(x = 0; x <= 1.01; x += 0.05)
    {
      fprintf(fp, "%.3f ", x);
      // loop through bottleneck size
      for(n = 1; n < 1000; n *= 2)
	{
	  p0 = 0.5; b = 1. - 1./n;

	  // compute pdf and output
	  r = phi(p0, b, x);
	  if(x > 0.005 && x < 0.995) r *= 0.05;
	  fprintf(fp, "%.4f ", r);
	}
      fprintf(fp, "\n");
    }
  fclose(fp);

  return 0;
}
      
