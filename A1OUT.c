/* Define the following four macros. */
#define DISCR_NAME  "A1OUT"
#define NAMEGEN(A) A1OUT ## A
#define NISPECS 1
#define NFSPECS 0

/* Do not alter the next 5 lines. */
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "parm_red_structs.h"
#include "readtoken.h"
/* Make any further definitions or includes that are necessary. */
#include <math.h>
#define S (discr_set->ispecs[0])

/* The following are the A1OUT (alpha on outside two points only) constants 
 * for beta, for step  sizes 1 through 17. Since the beta are symmetric, only the
 * first half are listed. (The first row is not used.) */
const static double A1OUT_b[18][9] = {
    {0.0},
    
    {0.5},
    
    {1.0/3.0, 4.0/3.0},
    
    {0.375, 1.125},
    
    {14.0/45.0, 64.0/45.0, 8.0/15.0},
    
    {95.0/288.0, 125.0/96.0, 125.0/144.0},
    
    {41.0/140.0, 54.0/35.0, 27.0/140.0, 68.0/35.0},
    
    {5257.0/17280.0, 25039.0/17280.0, 343.0/640.0, 20923.0/17280.0},
    
    {9836.0/33075.0, 1856.0/1225.0, 1184.0/4725.0, 9152.0/4725.0, 0.0},
    
    {25713.0/89600.0, 141669.0/89600.0, 243.0/2240.0, 10881.0/5600.0, 
        26001.0/44800.0},
    
    {27535.0/95256.0, 149375.0/95256.0, 3875.0/31752.0, 16375.0/7938.0, 
        0.0, 4345.0/2268.0},
    
    {11746361.0/40824000.0, 6379483.0/4032000.0, 815903.0/13063680.0,
        118459.0/54432.0, 0.0, 75732569.0/54432000.0},
    
    {653.0/2310.0, 15552.0/9625.0, 0.0, 1088.0/525.0, 243.0/350.0, 0.0,
        2336.0/875.0},
    
    {4881098833.0/17244057600.0, 3972529717.0/2463436800.0, 0.0,
        4659246007.0/2155507200.0, 844282933.0/3448811520.0, 79577537.0/70963200.0,
        773535139.0/718502400.0},
    
    {216835969.0/764478000.0, 55932373.0/34749000.0, 0.0, 894103.0/408375.0,
        548849.0/2673000.0, 2582447.0/2673000.0, 3111157.0/1782000.0, 0.0},
    
    {20597515.0/72855552.0, 9807885.0/6071296.0, 0.0, 7209025.0/3469312.0, 
        2364795.0/3469312.0, 0.0, 630445.0/236544.0, 390825.0/2207744.0},
    
    {46236616.0/164189025.0, 45051392.0/27785835.0, 0.0, 66059264.0/32837805.0,
        117107392.0/138929175.0, 0.0, 25996288.0/12629925.0, 14955008.0/12629925.0,
        0.0},
    
    {35662824992161.0/126931591987200.0, 441821701104101.0/271996268544000.0, 0.0,
        52647387834529.0/26444081664000.0, 9901044586361.0/11333177856000.0, 0.0,
        210391762153.0/96864768000.0, 1059993717659.0/1830744115200.0,
        3414363256259.0/3487131648000.0}
};

/* The only lines in the following function that need to be editted are the
 * lines between the two ************** marks. */
void NAMEGEN(_discr_set_init)(struct s_discr_set *discr_set, FILE *fptr) {
    
    int i,num;
    char sep;
    struct s_discr *discr;
    
    discr_set->name = DISCR_NAME;
    discr_set->n_ispecs = NISPECS;
    discr_set->ispecs = (int *) malloc(NISPECS*sizeof(int));
    discr_set->n_fspecs = NFSPECS;
    discr_set->fspecs = (double *) malloc(NFSPECS*sizeof(double));

    /* This adjusts the set-up depending on whether or not this discretization is the primary one chosen */
    if (discr_set->flag[0] == 0) {
        if (readtoken(fptr,NISPECS,"discretization_ispecs","%d",discr_set->ispecs,&sep) != NISPECS) {
            fprintf(stderr,"parm_reduction: error reading discretization_ispecs (A1OUT).\n");
            exit(1);
        }
        if (readtoken(fptr,NFSPECS,"discretization_fspecs","%lf",discr_set->fspecs,&sep) != NFSPECS) {
            fprintf(stderr,"parm_reduction: error reading discretization_fspecs.\n");
            exit(1);
        }
//        discr_set->ispecs[0] = discr_set->input_size; /* delete this later (only with output4) */
        discr_set->ispecs[0]--;
    }
    else {
        discr_set->ispecs[0] = discr_set->flag[0];
    }
    
    if (S > 17 || S < 1) {
        fprintf(stderr,"parm_reduction: error discr_set.ispec[0] must be between 2 and 18.\n");
        exit(1);
    }
    
    /* ***************** Edit the stuff below. */
    discr_set->n_discr = S;
    if ((discr_set->discr = (struct s_discr *) malloc(discr_set->n_discr*sizeof(struct s_discr))) == NULL) {
        fprintf(stderr,"parm_reduction: error allocating space for discr\n");
        exit(1);
    }

    for (num=1; num<=discr_set->n_discr; num++) {
        discr = &discr_set->discr[num-1];
        discr->discrtype = parm_bc;
        discr->winsize = num + 1;
        if ((discr->stencil = (enum e_stencil *) malloc((discr_set->n_discr+1+discr_set->flag[1])*
                                                        sizeof(enum e_stencil))) == NULL) {
            fprintf(stderr,"parm_reduction: error allocating space for stencil\n");
            exit(1);
        }
        /* Window stencils and constants: */
        if (((discr->alpha = (double *) malloc((discr_set->n_discr+1+discr_set->flag[1])*sizeof(double))) == NULL) ||
            ((discr->beta = (double *) malloc((discr_set->n_discr+1+discr_set->flag[1])*sizeof(double))) == NULL)) {
            fprintf(stderr,"parm_reduction: error allocating space for alpha/beta\n");
            exit(1);
        }
        discr->stencil[0] = var_vec;
        discr->alpha[0] = -1.0;
        discr->beta[0] = A1OUT_b[num][0];
        for (i=1; i<=(num)/2; i++) {
            discr->alpha[i] = 0.0;
            discr->alpha[num-i] = 0.0;
            discr->beta[i] = A1OUT_b[num][i];
            discr->beta[num-i] = A1OUT_b[num][i];
            if (fabs(A1OUT_b[num][i]) < DBL_EPSILON) {
                discr->stencil[i] = not_used;
                discr->stencil[num-i] = not_used;
            }
            else {
                discr->stencil[i] = vec;
                discr->stencil[num-i] = vec;
            }
        }
        discr->stencil[num] = var_vec;
        discr->alpha[num] = 1.0;
        discr->beta[num] = A1OUT_b[num][0];
        
    }
    /* ***************** Edit the stuff above. */
    return;
}
