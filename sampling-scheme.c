// Simple mtDNA bottleneck calculations
// sampling-scheme.c

// Simulates different sampling schemes for bottleneck calculations
// Single binomial sampling events (possibly with selection) are simulated in a large set of independent samples; statistics of resulting bottleneck estimates are reported as true bottleneck size, number of offspring, and level of selection vary

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()

int main(void)
{
  double n;
  int it;
  int i, j;
  double h[10000];
  double bneck1[10000], bneck2[10000], bneck3[10000];
  int nit = 1000;
  int noff = 20;
  double h0 = 0.501;
  double meanh, varh;
  double msd;
  double sel;
  double meanbneck1, meanbneck2, meanbneck3;
  double sdbneck1, sdbneck2, sdbneck3;
  double minbneck1, minbneck2, minbneck3;
  double maxbneck1, maxbneck2, maxbneck3;
  FILE *fp;
  
  srand48(122);

  fp = fopen("sampling-scheme.txt", "w");
  // loop through n, actual bottleneck size
  for(n = 1; n < 1000; n *= 1.2)
    {
      // loop through sel, level of selection
      for(sel = 1; sel >= 0.5; sel /= 1.1)
	{
	  // loop through noff, number of offspring in a litter
	  for(noff = 1; noff <= 100; noff *= 2)
	    {
	      // nit is number of iterations (independent samples)
	      for(it = 0; it < nit; it++)
		{
		  // loop through number of offspring and construct heteroplasmies
		  for(j = 0; j < noff; j++)
		    {
		      // construct new heteroplasmy by applying one round of binomial sampling, possibly with selection (if sel < 1 we preferentially lose "mutant" mtDNA)
		      h[j] = 0;
		      for(i = 0; i < n; i++)
			{
			  if(RND < h0*sel) h[j]++;
			}
		      h[j] /= n;
		    }

		  // calculate mean-square-difference and offspring stats 
		  msd = meanh = varh = 0;
		  for(j = 0; j < noff; j++)
		    {
		      msd += (h[j]-h0)*(h[j]-h0);
		      meanh += h[j];
		    }
		  meanh /= noff;
		  msd /= noff;
		  for(j = 0; j < noff; j++)
		    {
		      varh += (meanh-h[j])*(meanh-h[j]);
		    }
		  varh /= noff-1;

		  // estimates of bottleneck size using different calculations
		  bneck1[it] = h0*(1.-h0)/msd;
		  bneck2[it] = h0*(1.-h0)/varh;
		  bneck3[it] = meanh*(1.-meanh)/varh;
		}

	      // we've estimated bottleneck size for nit independent samples
	      // now compute statistics of bottleneck size estimates across this set
	      // not beautiful code...
	      meanbneck1 = meanbneck2 = meanbneck3 = sdbneck1 = sdbneck2 = sdbneck3 = 0;
	      maxbneck1 = maxbneck2 = maxbneck3 = 0; minbneck1 = minbneck2 = minbneck3 = 999999999;
	      for(it = 0; it < nit; it++)
		{
		  meanbneck1 += bneck1[it];
		  meanbneck2 += bneck2[it];
		  meanbneck3 += bneck3[it];
		  if(bneck1[it] < minbneck1) minbneck1 = bneck1[it];
		  if(bneck2[it] < minbneck2) minbneck2 = bneck2[it];
		  if(bneck3[it] < minbneck3) minbneck3 = bneck3[it];
		  if(bneck1[it] > maxbneck1) maxbneck1 = bneck1[it];
		  if(bneck2[it] > maxbneck2) maxbneck2 = bneck2[it];
		  if(bneck3[it] > maxbneck3) maxbneck3 = bneck3[it];
		}
	      meanbneck1 /= nit;
	      meanbneck2 /= nit;
	      meanbneck3 /= nit;
	      for(it = 0; it < nit; it++)
		{
		  sdbneck1 += (bneck1[it]-meanbneck1)*(bneck1[it]-meanbneck1);
		  sdbneck2 += (bneck2[it]-meanbneck2)*(bneck2[it]-meanbneck2);
		  sdbneck3 += (bneck3[it]-meanbneck3)*(bneck3[it]-meanbneck3);
		}
	      sdbneck1 = sqrt(sdbneck1/(nit-1));
	      sdbneck2 = sqrt(sdbneck2/(nit-1));
	      sdbneck3 = sqrt(sdbneck3/(nit-1));

	      // output
	      fprintf(fp, "%f %f %i %f %f %f %f %f %f %f %f %f %f %f %f\n", n, sel, noff, meanbneck1, sdbneck1, minbneck1, maxbneck1, meanbneck2, sdbneck2, minbneck2, maxbneck2, meanbneck3, sdbneck3, minbneck3, maxbneck3);
	    }
	}
    }
  fclose(fp);

  return 0;
}
