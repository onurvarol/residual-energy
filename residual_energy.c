/*******************************************************************************
* Residual energy contribution calculation.
* Onur Varol - 21.03.2012
*
* residual_energy -nc <Nresidue> -f <Nframe> -cmap <cmap file> -dr <drfile> -out <outputfile>
* nc and f are number of row and column drfile
* cmap file shold have 0-1 file
*******************************************************************************/

#include<stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double **cmap;
double **dr;
double **U;

double** Make2DDoubleArray(int arraySizeX, int arraySizeY) 
{  
    int i;     
    double** theArray;  
    theArray = (double**) malloc(arraySizeX*sizeof(double*));  
    for (i = 0; i < arraySizeX; i++)
    {  
       theArray[i] = (double*) malloc(arraySizeY*sizeof(double));  
    }
       return theArray;  
} 

main(int argc, char **argv)
{
    int i, j, k, l, t,
    count, Nframe, Nresidue;  
    double temp, ttemp, dr_kl;
    char* cmapfile = "cmap_temp.pr";
    char* drfile = "dr_all.pr";
    char* outfile = "U_all.pr";

    for(i=1; i<argc; i++)
    {
        if(i+1 != argc)
        {
           if(!strcmp(argv[i], "-nc")){
                   Nresidue = atoi(argv[i+1]);}
           if(!strcmp(argv[i], "-f")){
                   Nframe = atoi(argv[i+1]);}
           if(!strcmp(argv[i], "-cmap")){
                   cmapfile = argv[i+1];} 
           if(!strcmp(argv[i], "-dr")){
                   drfile = argv[i+1];}   
           if(!strcmp(argv[i], "-out")){
                   outfile = argv[i+1];} 
        }     
    }
    printf("%d Ca - %d frame\n",Nresidue,Nframe);
    //fprintf(stdout,"%d Ca - %d frame\n",Nresidue,Nframe);

    double *dr_ij = (double*)malloc(Nframe*sizeof(double));
    
    cmap = Make2DDoubleArray(Nresidue,Nresidue);
    dr = Make2DDoubleArray(3*Nresidue,Nframe);
    U = Make2DDoubleArray(Nresidue,Nresidue);
    printf("Matrices created\n");
	
    /*----- Read Matrices -----*/
    FILE *fdr;
    FILE *fcmap;
    if(!(fdr = fopen(drfile ,"r"))){
          fprintf(stderr,"file:dr not found\n");
          exit(-1);
    }
    for(i=0; i<Nresidue; i++)
    {
        for(j=0; j<Nframe; j++)
        {
            fscanf(fdr, "%lf ", &dr[i][j]);     
        }
    }
    fclose(fdr);
    
    if(!(fcmap = fopen(cmapfile ,"r"))){
          fprintf(stderr,"file:cmap not found\n");
          exit(-1);
    }
    for(i=0; i<3*Nresidue; i++)
    {
        for(j=0; j<Nresidue; j++)
        {
            fscanf(fcmap, "%lf", &cmap[i][j]);
        }
    }
    fclose(fcmap);         
         
    printf("Files read finished\n");
    //fprintf(stdout,"Files read finished\n");
    fflush(stdout);	

    #pragma omp parallel
    for(i=0; i<Nresidue; i++)
    {
       for(k=i; k<Nresidue; k++)
       {
           temp = 0;
           for(j=0; j<Nresidue; j++)
           {
              if(cmap[i][j] > 0)
              {
                 for(t=0; t<Nframe; t++)
                 {
                     dr_ij[t] = (dr[i][t] - dr[j][t]) * (dr[i][t] - dr[j][t]) + 
                                (dr[i + Nresidue][t] - dr[j + Nresidue][t]) * (dr[i + Nresidue][t] - dr[j + Nresidue][t]) + 
                                (dr[i + 2 * Nresidue][t] - dr[j + 2 * Nresidue][t]) * (dr[i + 2 * Nresidue][t] - dr[j + 2 * Nresidue][t]);     
                 }           
              }
              
              for(l=0; l<Nresidue; l++)
              {
                 ttemp = 0;
                 if(cmap[i][j]*cmap[k][l] > 0)
                 {
                    dr_kl = 0;
                    for(t=0; t<Nframe; t++)
                    {
                        dr_kl = (dr[k][t] - dr[l][t]) * (dr[k][t] - dr[l][t]) + 
                                (dr[k + Nresidue][t] - dr[l + Nresidue][t]) * (dr[k + Nresidue][t] - dr[l + Nresidue][t]) + 
                                (dr[k + 2 * Nresidue][t] - dr[l + 2 * Nresidue][t]) * (dr[k + 2 * Nresidue][t] - dr[l + 2 * Nresidue][t]);
                        ttemp += dr_ij[t] * dr_kl;
                    }
                    temp += ttemp / (Nframe*1.0);
                 }             
              }     
           }
           U[i][k] = temp;
           printf(".");     
       }    
       printf("%d-%d \n",i,Nresidue);
	   fflush(stdout);
       //fprintf(stdout,"%d-%d \n",i,Nresidue); 
    }
    printf("Finished, saving file\n");
    FILE *u;
    if(!(u = fopen(outfile ,"w+"))){
          fprintf(stderr,"file:output not found\n");
          exit(-1);
    }
    for(i=0; i<Nresidue; i++)
    {
        for(j=i; j<Nresidue; j++)
        {
            fprintf(u,"%d %d %f\n",i,j,U[i][j]);         
        }     
        //fprintf(u,"\n");
    }
    fclose(u);
    free(cmap);
    free(U);
    free(dr);
    return(0);
}
