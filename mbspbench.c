/*

The MIT License (MIT)

Copyright (c) 2014 Marcelo Alaniz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


#include "mbsputil.h"
#include "mbsp-discover.h"

#define MANUAL_AFFINITY 2
multibsp_tree_node_t tn2;
multibsp_tree_node_t tnode;

#define NITERS 10     /* number of iterations */
#define MAXN 256     /* maximum length of DAXPY computation */
#define MAXH 2048       /* maximum h in h-relation */

#define MEGA 1000000.0


int P; /* number of processors requested */

void leastsquares(int h0, int h1, double *t, double *g, double *l){
    /* This function computes the parameters g and l of the 
       linear function T(h)= g*h+l that best fits
       the data points (h,t[h]) with h0 <= h <= h1. */

    double nh, sumt, sumth, sumh, sumhh, a;
    int h;

    nh= h1-h0+1;
    /* Compute sums:
        sumt  =  sum of t[h] over h0 <= h <= h1
        sumth =         t[h]*h
        sumh  =         h
        sumhh =         h*h     */
    sumt= sumth= 0.0;
    for (h=h0; h<=h1; h++){
        sumt  += t[h];
        sumth += t[h]*h;
    }
    sumh= (h1*h1-h0*h0+h1+h0)/2;  
    sumhh= ( h1*(h1+1)*(2*h1+1) - (h0-1)*h0*(2*h0-1))/6;

    /* Solve      nh*l +  sumh*g =  sumt 
                sumh*l + sumhh*g = sumth */
    if(fabs(nh)>fabs(sumh)){
        a= sumh/nh;
        /* subtract a times first eqn from second eqn */
        *g= (sumth-a*sumt)/(sumhh-a*sumh);
        *l= (sumt-sumh* *g)/nh;
    } else {
        a= nh/sumh;
        /* subtract a times second eqn from first eqn */
        *g= (sumt-a*sumth)/(sumh-a*sumhh);
        *l= (sumth-sumhh* *g)/sumh;
    }

} /* end leastsquares */


void bspbench(){
    void leastsquares(int h0, int h1, double *t, double *g, double *l);
    int p, s, s1, iter, i, n, h, destproc[MAXH], destindex[MAXH];
    double alpha, beta, x[MAXN], y[MAXN], z[MAXN], src[MAXH], *dest,
           time0, time1, time, *Time, mintime, maxtime,
           nflops, r, g0, l0, g, l, t[MAXH+1]; 
  
    size_t pin[100];

    // Determine p 
    // start: new code for pinning
    for (i=0; i< tnode->length; i++) pin[i] = tnode->sons[i]->index;
    mcbsp_set_pinning( pin, tnode->length );
    bsp_begin(tnode->length);
    // end: new code for pinning

    p= bsp_nprocs(); // p = number of processors obtained 
    s= bsp_pid();    // s = processor number

    Time= vecallocd(p); bsp_push_reg(Time,p*SZDBL);
    dest= vecallocd(2*(MAXH+p)); bsp_push_reg(dest,(2*(MAXH+p))*SZDBL);
    bsp_sync();

    // Determine r 

    for (n=1; n < MAXN; n *= 2){
        // Initialize scalars and vectors 
        alpha= 1.0/3.0;
        beta= 4.0/9.0;
        for (i=0; i<n; i++){
          z[i]= y[i]= x[i]= (double)i;
        }
        // Measure time of 2*NITERS DAXPY operations of length n 
        time0=bsp_time();
        for (iter=0; iter<NITERS; iter++){
          for (i=0; i<n; i++)
            y[i] += alpha*x[i];
          for (i=0; i<n; i++)
            z[i] -= beta*x[i];
        }
        time1= bsp_time(); 
        time= time1-time0; 
        bsp_put(0,&time,Time,s*SZDBL,SZDBL);
        bsp_sync();

        // Processor 0 determines minimum, maximum, average computing rate 
        if (s==0){
          mintime= maxtime= Time[0];
          for(s1=1; s1<p; s1++){
            mintime= MIN(mintime,Time[s1]);
            maxtime= MAX(maxtime,Time[s1]);
          }
          if (mintime>0.0){
            // Compute r = average computing rate in flop/s 
            nflops= 4*NITERS*n;
            r= 0.0;
            for(s1=0; s1<p; s1++)
              r += nflops/Time[s1];
            r /= p; 

            //printf("n= %5d min= %7.3lf max= %7.3lf av= %7.3lf Mflop/s ",
            //       n, nflops/(maxtime*MEGA),nflops/(mintime*MEGA), r/MEGA);
            //fflush(stdout);
            // Output for fooling benchmark-detecting compilers 
            printf( "", y[n-1]+z[n-1] );
          } 
        }
      }

      // Determine g and l 
      for (h=0; h<=MAXH; h++){
        // Initialize communication pattern 
        for (i=0; i<h; i++){
          src[i]= (double)i;
          if (p==1){
            destproc[i]=0;
            destindex[i]=i;
          } else {
            // destination processor is one of the p-1 others 
            destproc[i]= (s+1 + i%(p-1)) %p;
            // destination index is in my own part of dest 
            destindex[i]= s + (i/(p-1))*p;
          }
        }
        for (i=0; i<h; i++){
          src[i]= (double)i;
          if (p==1){
            destproc[i]=0;
            destindex[i]=i;
          } else {
            // destination processor is one of the p-1 others 
            destproc[i]= (s+1 + i%(p-1)) %p;
            // destination index is in my own part of dest 
            destindex[i]= s + (i/(p-1))*p;
          }
        }


        // Measure time of NITERS h-relations 
        bsp_sync(); 

        time0= bsp_time(); 
        for (iter=0; iter<NITERS; iter++){
          for (i=0; i<h; i++) {
            //bsp_get(0,  dest, destindex[i]*SZDBL, &src[i] , SZDBL);
            //bsp_get(destproc[i],  dest, destindex[i]*SZDBL, &src[i] , SZDBL);
            bsp_put(destproc[i],  &src[i] , dest              , destindex[i]*SZDBL, SZDBL);
          }

          //if (s == 0) 
          //  bsp_get(0,  dest, destindex[i]*SZDBL, &src[i] , SZDBL);

          bsp_sync(); 
          
        }

        time1= bsp_time();
        time= time1-time0;

        // Compute time of one h-relation 
        if (s==0){
          t[h]= (time*r)/NITERS;
//#define SEHLOC_BENCH_VERBOSE
#ifdef SEHLOC_BENCH_VERBOSE
          char strnodes[256];
          sprintf(strnodes, "");
          for (i=0; i<tnode->length; i++) {
            sprintf(strnodes, "%s %d", strnodes, tnode->sons[i]->index);
          }
          printf("SEH# Level%d %5d %lf %8.0lf\n", tnode->level, h, time/NITERS, t[h]); fflush(stdout);
#endif
        }
      }

      if (s==0){
        leastsquares(0,p,t,&g0,&l0); 
        printf("Range h=0 to p   : g= %.1lf, l= %.1lf\n",g0,l0);
        leastsquares(p,MAXH,t,&g,&l);
        g=(g>0)? g: g0*2;
        printf("Range h=p to HMAX: g= %.1lf, l= %.1lf\n",g,l);
        //printf("plot# %d %.1lf  %.1lf\n",tnode->level, g,l);
        printf("The bottom line for this MultiBSP component is:\n");
        printf("<p= %d, r= %.3lf Mflop/s, g= %.1lf, l= %.1lf>\n",
               p,r/MEGA,g,l);
        fflush(stdout);
      }
      bsp_pop_reg(dest); vecfreed(dest);
      bsp_pop_reg(Time); vecfreed(Time);
      bsp_end();

} /* end bspbench */


int main(int argc, char **argv){
    multibsp_tree_node_t  tree_root = NULL;
    int i;
    char* strp;


    tree_root = multibsp_discover_new();

    strp = multibsp_discover_print(tree_root);
    printf("%s\n", strp);
    mcbsp_set_affinity_mode( MANUAL_AFFINITY ); // go MANUAL affinity

    printf("\nMBSPDiscover Benchmark\n");
    //printf("G   \t\tLEVEL\tP  \tR(Mflop/s)\tL \n");
    for (tnode=tree_root; tnode->length != 0; tnode = tnode->sons[0]) {
      printf("\n***************** LEVEL %d *******************\n", tnode->level);
      bsp_init(bspbench, argc, argv);
      bspbench();
    }

    /*
    for (tnode=tree_root; tnode->length > 0; ) {
      for (i =0; i< tnode->length; i++) {
        tnode = tnode->sons[i];
        printf("***************** LEVEL %d *******************\n", tnode->level);
        bsp_init(bspbench, argc, argv);
        bspbench();
      }
    }
    */

    printf("\n");
    exit(0);
} /* end main */
