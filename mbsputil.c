#include "mbsputil.h"

double *vecallocd(int n)
{
    double *pd;

    if (n==0){
        pd= NULL;
    } else {
        pd= (double *)malloc(n*SZDBL);
        if (pd==NULL)
            bsp_abort("vecallocd: not enough memory");
    }
    return pd;
} 

int *vecalloci(int n)
{
    int *pi;

    if (n==0){
        pi= NULL; 
    } else { 
        pi= (int *)malloc(n*SZINT);
        if (pi==NULL)
            bsp_abort("vecalloci: not enough memory");
    }
    return pi;
} 

double **matallocd(int m, int n)
{
    int i;
    double *pd, **ppd;

    if (m==0){
        ppd= NULL;  
    } else { 
        ppd= (double **)malloc(m*sizeof(double *));
        if (ppd==NULL)
            bsp_abort("matallocd: not enough memory");
        if (n==0){
            for (i=0; i<m; i++)
                ppd[i]= NULL;
        } else {  
            pd= (double *)malloc(m*n*SZDBL); 
            if (pd==NULL)
                bsp_abort("matallocd: not enough memory");
            ppd[0]=pd;
            for (i=1; i<m; i++)
                ppd[i]= ppd[i-1]+n;
        }
    }
    return ppd;
} 

void vecfreed(double *pd)
{

    if (pd!=NULL)
        free(pd);

}

void vecfreei(int *pi)
{
    if (pi!=NULL)
        free(pi);

} 

void matfreed(double **ppd)
{
    if (ppd!=NULL){
        if (ppd[0]!=NULL)
            free(ppd[0]);
        free(ppd);
    }
}
