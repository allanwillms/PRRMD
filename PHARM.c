/* This file contains a template for the model-specific information and code.
 * The ??? symbols should be replaced with appropriate code and then the file 
 * s_parm_reduction.h should be editted to add this model to the list in the
 * variable all_models. Then the program should be recompiled. */

/* First define a prefix for the function names which is unique to this model.
 * Put this in quotes in the definition of MODEL_NAME, and put it without quotes
 * in the definition of NAMEGEN. */
#define MODEL_NAME  "PHARM"
#define NAMEGEN(A) PHARM ## A

/* Next define values for the model: */
#define DEDIM  4   /* Enter the dimension of the Differential Equation. */
#define OFDIM  0   /* Enter the dimension of other functions. */
#define AEDIM  0   /* Enter the dimension of the Algebraic State Equations */
#define NPARMS  7  /* Enter the number of parameters in the DE. */
#define NISPECS  0  /* Enter number of integer specs for the model. */
#define NFSPECS  1  /* Enter number of floating point specs for the model. */

/* The user must code two external variables and two functions for the MODEL.
 * The two variables simply define whether each equation in the model involves 
 * each parameter or variable.
 * The two functions are:
 * 1) A dependence function which determines whether the vector field is monotonic
 *    with respect to the parameters.
 * 2) A function which computes f + ay where f is the vector field.
 */

/* TO DO in this file
 * 1) Unmeasured monotonicity for variables 1,2,3,4
 * 2) Beta hull consistency for unmeasured variables 
 * 3) Fix equation 0 (used wrong version) */

/* Do not change the next 4 lines. */
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "parm_red_structs.h"
#include "readtoken.h"
/* Add any other include files or defines needed for the coded functions. */
#include <math.h>
#include "interval_math.h"
#define TOTALDRUG (model->fspecs[0])
#define VB p[0]
#define k1 p[1]
#define K2 p[2]
#define kM p[3]
#define muD p[4]
#define muM p[5]
#define ET p[6]
#define DB y[0]
#define MB y[1]
#define DU y[2]
#define MU y[3]


/* PHARM:  pharmacological model of drug taken orally and removed through urine
 * The overall model is
 *
 * DB' = (I'(t) - (k1 + muD + kM)*DB + k2*DT)/VB
 * DT' = (k1*DB - k2*DT)/VT
 * MB' = (kM*DB - muM*MB)/VB
 * DU' = muD*DB
 * MU' = muM*MB
 *
 * where:
 * DB is the drug in the blood (mg/litre)
 * DT is the drug in the tissue (mg/litre)
 * DU is the cumulative drug in the urine (mg)
 * MB is the metabolite in the blood (mg/litre)
 * MU is the cumulative metabolite in the urine (mg)
 * time is in hours
 * VB is the effective blood volume (litres)
 * VT is the effective tissue volume (litres)
 * k1 and k2 are the rates of transfer for the drug between the blood
 *    and the tissue  (litres/hour)
 * kM is the rate of metabolism of drug in the blood (litres/hour)
 * muD and muM are the rates of output of the drug from the bloood
 *    to the urine (litres/hour)
 * I'(t) is the input rate: TOTALDRUG/ET  if 0<=t<ET,  0 if t>=ET, so that
 *    the total drug administered is TOTALDRUG by time t=ET.
 *
 * The redundancy in the model is
 * (VB*(DB+MB) + VT*DT + DU + MU)' = I'(t)
 * so that
 * VB*(DB+MB) + VT*DT + DU + MU = I(t) = [ TOTALDRUG*t/ET  if 0<=t<ET
 *                                       [ TOTALDRUG       if t>=ET
 * We use this to remove VT*DT from the model (since DT is not measured) by
 * defining K2 = k2/VT.
 * Thus the resulting 4-dimensional model is:
 *
 * DB' = (I'(t) - (k1 + muD + kM)*DB + K2*(I(t) - DU - MU))/VB - K2*(DB + MB)
 * MB' = (kM*DB - muM*MB)/VB
 * DU' = muD*DB
 * MU' = muM*MB
 *   */

/* Specify whether each equation involves the parameters and variables and other
 * functions.  1 indicates this equation depends on the parameter/variable, 
 * 0 indicates no. */
static int PARMDEP[DEDIM+AEDIM][NPARMS] = {
    {1, 1, 1, 1, 1, 0, 1},
    {1, 0, 0, 1, 0, 1, 0},
    {0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 1, 0}
};

static int VARDEP[DEDIM+AEDIM][DEDIM] = {
    {1, 1, 1, 1},
    {1, 1, 0, 0},
    {1, 0, 0, 0},
    {0, 1, 0, 0}
};

static int VARPARTIALDEP[DEDIM+AEDIM][DEDIM+1] = {
    {1, 1, 1, 1, 1},
    {1, 1, 0, 0, 0},
	{0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
};

static int OFDEP[DEDIM+AEDIM][OFDIM] = {
    {},
    {},
    {},
    {}
};

static int OFPARTIALDEP[DEDIM+AEDIM][OFDIM] = {
    {},
    {},
	{},
    {}
};

/* Specify any a priori extreme bounds on the parameters.  Use +-DBL_MAX to
 * represent +- infinity, if there are no bounds, and +-DBL_MIN for the numbers
 * closest to zero. */
static double P_EXTREME[NPARMS][2] = {
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX}
};

static double Y_EXTREME[DEDIM][2] = {
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX}
};

void NAMEGEN(_determine_dependence)(int eqn, Range *p, Range t, Range *y, Range *u, Range partial,
                                    struct s_model *model, int *dep, struct s_discr *discr, int parm) {
    
    Range It,Itprime,z,w;
    /* dep[i]=1 if equation eqn of the vector field is monotonically nondecreasing
     * with p[i] for each range of t, y, and u values.  -1 if it is monotonically 
     * nonincreasing, 0 if it does not depend on p[i], and -2 if it is not monotonic. 
     * Return a value of 0 if any of dep entries are -2, else return 1. */
    
    switch (eqn) {
        case 0:
            /* Programming It and Itprime for use in eqn 0 */
            if (t[0] >= ET[1]) {
                It[0] = TOTALDRUG;
                It[1] = TOTALDRUG;
                Itprime[0] = 0.0;
                Itprime[1] = 0.0;
            }
            else if (t[1] < ET[0]) {
                It[0] = TOTALDRUG*t[0]/ET[1];
                It[1] = TOTALDRUG*t[1]/ET[0];
                Itprime[0] = TOTALDRUG/ET[1];
                Itprime[1] = TOTALDRUG/ET[0];
            }
            else {
                It[0] = TOTALDRUG*t[0]/ET[1];
                It[1] = TOTALDRUG;
                Itprime[0] = 0.0;
                Itprime[1] = TOTALDRUG/ET[0];
            }
            
            /* Now determine monotonicities */
            switch (parm) {
                case 0:
                    z[0] = Itprime[0] - (k1[1] + muD[1] + kM[1])*DB[1] + K2[0]*(It[0] - DU[1] - MU[1]);
                    z[1] = Itprime[1] - (k1[0] + muD[0] + kM[0])*DB[0] + K2[1]*(It[1] - DU[0] - MU[0]);
                    int_multiply(VB,VB,partial);
                    int_divide(z,partial,partial);
                    multiply(-1.0,partial);
                    if (ISZERO_INTERIOR(z)) dep[parm] = -2;
                    else if (z[0] + z[1] > 0.0) dep[parm] = -1;
                    else dep[parm] = 1;
                    break;
                case 1:
                    dep[parm] = -1;
                    break;
                case 2:
                    z[0] = (It[0] - DU[1] - MU[1])/VB[1] - DB[1] - MB[1];
                    z[1] = (It[1] - DU[0] - MU[0])/VB[0] - DB[0] - MB[0];
                    int_put(z,partial);
                    if (ISZERO_INTERIOR(z)) dep[parm] = -2;
                    else if (z[0] + z[1] > 0.0) dep[parm] = 1;
                    else dep[parm] = -1;
                    break;
                case 3:
                    dep[parm] = -1;
                    break;
                case 4:
                    dep[parm] = -1;
                    break;
                case 6:
//                    dep[parm] = -2;
                    if (t[0] >= ET[1]) dep[parm] = -2;
                    else dep[parm] = -1;
                    break;
                default: /* variable */
                    if (model->var == 0) { 
                        int_add(k1,muD,partial);
                        int_add(partial,kM,partial);
                        int_divide(partial,VB,partial);
                        multiply(-1.0,partial);
                        int_subtract(partial,K2,partial);
                    }
                    else if (model->var == 1) { 
                        int_put(K2,partial);
                        multiply(-1.0,partial);
                    } 
                    else { /* either w.r.t. DU or MU */
                        int_divide(K2,VB,partial);
                        multiply(-1.0,partial);
                    }
            }
            break;
        case 1: /* eqn == 1 */
            switch (parm) {
                case 0:
                    int_multiply(kM,DB,z);
                    int_multiply(muM,MB,w);
                    int_subtract(z,w,z);
                    int_multiply(VB,VB,w);
                    int_divide(partial,w,partial);
                    multiply(-1.0,partial);
                    if (ISZERO_INTERIOR(z)) dep[parm] = -2;
                    else if (z[0] + z[1] > 0.0) dep[parm] = -1;
                    else dep[parm] = 1;
                    break;
                case 3:
                    dep[parm] = 1;
                    break;
                case 5: 
                    dep[parm] = -1;
                    break;
                default: /* variable */
                    if (model->var == 0) { 
                        int_divide(kM,VB,partial);
                    }
                    else {
                        int_divide(muM,VB,partial);
                        multiply(-1.0,partial);
                    }
            }
            break;
        case 2: /* eqn == 2 */
            switch (parm) {
                case 4:
                    dep[parm] = 1;
                    break;
                default: /* variable */
                    if (model->var == 0) { 
                        int_put(muD,partial);
                    }
                    else {
                        partial[0] = 0.0;
                        partial[1] = 0.0;
                    }
            }
            break;
        default: /* eqn == 3 */
            switch (parm) {
                case 5:
                    dep[parm] = 1;
                    break;
                default: /* variable */
                    if (model->var == 1) { 
                        int_put(muM,partial);
                    }
                    else {
                        partial[0] = 0.0;
                        partial[1] = 0.0;
                    }
            }
            break;
    }

    return;
}

double NAMEGEN(_vec_field)(int cor, int eqn, double *p, double a, double t, 
                           Range *y, Range *u, struct s_model *model) {
    /* This function returns either the lower or upper value (depending on whether
     * cor is 0 or 1) of the eqn'th entry of f + ay, where f is the vector
     * field.
     * */
        
    int j,opp_cor;
    double It,Itprime,g,temp;
    
    opp_cor = (cor == 0) ? 1 : 0;
        
    switch (eqn) {
        case 0:
            if (t >= ET) {
                It = TOTALDRUG;
                Itprime = 0.0;
            }
            else {
                It = TOTALDRUG*t/ET;
                Itprime = TOTALDRUG/ET;
            }
            temp = -(k1 + muD + kM)/VB - K2 + a;
            j = (temp >= 0.0) ? cor : opp_cor;
            g = (Itprime + K2*(It - DU[opp_cor] - MU[opp_cor]))/VB - K2*MB[opp_cor] + temp*DB[j];
            break;
        case 1:
            temp = -muM/VB + a;
            j = (temp >= 0.0) ? cor : opp_cor;
            g = kM*DB[cor]/VB + temp*MB[j];
            break;
        case 2:
            j = (a > 0.0) ? cor : opp_cor;
            g = muD*DB[cor] + a*DU[j];
            break;
        default: /* eqn == 3 */
            j = (a > 0.0) ? cor : opp_cor;
            g = muM*MB[cor] + a*MU[j];
            break;
    }
    return g;
}

double NAMEGEN(_interval_vec_field)(int eqn, int cor, Range *p, double *t, Range **y, Range **u, 
                                    struct s_model *model, struct s_discr *discr) {
    /* This function returns either the lower or upper value (depending on whether
     * cor is 0 or 1) of the eqn'th equation of f + ay, where f is the vector
     * field.
     *
     * The parameters p are assumed nonnegative. */
    /* g is where the final discretization value goes */

    int i;
    Range *z, *It, *Itprime;;
    if (((z = (Range *) malloc(discr->sub_window[1]*sizeof(Range))) == NULL) ||
        ((It = (Range *) malloc(discr->sub_window[1]*sizeof(Range))) == NULL) ||
        ((Itprime = (Range *) malloc(discr->sub_window[1]*sizeof(Range))) == NULL)) {
        fprintf(stderr,"parm_reduction: error allocating memory for z,It,Itprime\n"); exit(1);
    }
    Range g,w,v;  
    
    switch (eqn) {
        case 0:
            for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++){
                if (t[i] < ET[0]){
                    Itprime[i][0] = TOTALDRUG/ET[1];
                    Itprime[i][1] = TOTALDRUG/ET[0];
                }
                else if (t[i] < ET[1]) {
                    Itprime[i][0] = 0.0;
                    Itprime[i][1] = TOTALDRUG/ET[0];
                }
                else {
                    Itprime[i][0] = 0.0;
                    Itprime[i][1] = 0.0;
                }
            }
            int_vfield_sum(Itprime,w,discr);
            int_put(w,g);
            
            int_vfield_sum(DB,w,discr);
            int_add(k1,muD,v);
            int_add(v,kM,v);
            multiply(-1.0,v);
            int_multiply(v,w,w);
            int_add(g,w,g);

            for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++){
                if (t[i] < ET[0]){
                    It[i][0] = TOTALDRUG*t[i]/ET[1];
                    It[i][1] = TOTALDRUG*t[i]/ET[0];
                }
                else if (t[i] < ET[1]) {
                    It[i][0] = TOTALDRUG*t[i]/ET[1];
                    It[i][1] = TOTALDRUG;
                }
                else {
                    It[i][0] = TOTALDRUG;
                    It[i][1] = TOTALDRUG;
                }
                int_subtract(It[i],DU[i],z[i]);
                int_subtract(z[i],MU[i],z[i]);
            }
            int_vfield_sum(z,w,discr);
            int_multiply(K2,w,w);
            int_add(g,w,g);
            
            int_divide(g,VB,g);
            
            for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) int_add(DB[i],MB[i],z[i]);
            int_vfield_sum(z,w,discr);
            int_multiply(K2,w,w);
            multiply(-1.0,w);
            int_add(g,w,g);
            break;
        case 1:
            int_vfield_sum(DB,w,discr);
            int_multiply(kM,w,w);
            int_put(w,g);
            
            int_vfield_sum(MB,w,discr);
            int_multiply(muM,w,w);
            multiply(-1.0,w);
            int_add(g,w,g);
            
            int_divide(g,VB,g);
            break;
        case 2:
            int_vfield_sum(DB,w,discr);
            int_multiply(muD,w,w);
            int_put(w,g);
            break;
        default: /* eqn == 3 */
            int_vfield_sum(MB,w,discr);
            int_multiply(muM,w,w);
            int_put(w,g);
            break;
    }
    
    free(z);
    free(It);
    free(Itprime);
    return g[cor];
}

int NAMEGEN(_invert_vec_field)(int eqn, int i, Range F, struct s_model *model, Range *p, double t, Range *y, Range *u) {
    
    /* in this function, we are solving for parameter i in equation eqn.
     * If we start with f^s(p_i,**) = F (f^s is the vector field at the s-th time step, 
     * and F is the input to the function then we must invert f^s, 
     * and find the expression p_i = (f^s)^(-1)(F,**) */ 
    
    int ok = 0;
    Range hull,K,H,It,Itprime;
    
    switch (eqn) {
        case 0:
            /* Programming It and Itprime */
            if (t >= ET[1]) {
                It[0] = TOTALDRUG;
                It[1] = TOTALDRUG;
                Itprime[0] = 0.0;
                Itprime[1] = 0.0;
            }
            else if (t < ET[0]) {
                It[0] = TOTALDRUG*t/ET[1];
                It[1] = TOTALDRUG*t/ET[0];
                Itprime[0] = TOTALDRUG/ET[1];
                Itprime[1] = TOTALDRUG/ET[0];
            }
            else {
                It[0] = TOTALDRUG*t/ET[1];
                It[1] = TOTALDRUG;
                Itprime[0] = 0.0;
                Itprime[1] = TOTALDRUG/ET[0];
            }
            
            switch (i) {
                case 0: /* parameter VB */
                    int_add(DB,MB,K);
                    int_multiply(K2,K,K);
                    int_add(F,K,hull);
                    
                    int_add(k1,muD,K);
                    int_add(K,kM,K);
                    int_multiply(K,DB,K);
                    int_subtract(Itprime,K,K);
                    
                    int_subtract(It,DU,H);
                    int_subtract(H,MU,H);
                    int_multiply(K2,H,H);
                    int_add(K,H,K);
                                        
                    int_divide(hull,K,hull);
                    ok = 1;
                    break;
                case 1: /* parameter k1 */
                    int_add(DB,MB,K);
                    int_multiply(K2,K,K);
                    int_add(F,K,hull);
                
                    int_multiply(hull,VB,hull);
                    
                    int_subtract(hull,Itprime,hull);
                    
                    int_subtract(It,DU,K);
                    int_subtract(K,MU,K);
                    int_multiply(K,K2,K);
                    int_subtract(hull,K,hull);
                    
                    int_divide(hull,DB,hull);
                    multiply(-1.0,hull);
                    
                    int_subtract(hull,muD,hull);
                    int_subtract(hull,kM,hull);
                    ok = 1;
                    break;
                case 2: /* parameter K2 */
                    int_multiply(F,VB,hull);
                    
                    int_subtract(hull,It,hull);
                    
                    int_add(k1,muD,K);
                    int_add(K,kM,K);
                    int_multiply(K,DB,K);
                    int_add(hull,K,hull);
                    
                    int_subtract(It,DU,K);
                    int_subtract(K,MU,K);
                    int_divide(K,VB,K);
                    int_subtract(K,DB,K);
                    int_subtract(K,MB,K);
                    int_divide(hull,K,hull);
                    ok = 1;
                    break;
                case 3: /* parameter kM */
                    int_add(DB,MB,K);
                    int_multiply(K2,K,K);
                    int_add(F,K,hull);
                    
                    int_multiply(hull,VB,hull);
                    
                    int_subtract(hull,Itprime,hull);
                    
                    int_subtract(It,DU,K);
                    int_subtract(K,MU,K);
                    int_multiply(K,K2,K);
                    int_subtract(hull,K,hull);
                    
                    int_divide(hull,DB,hull);
                    multiply(-1.0,hull);
                    
                    int_subtract(hull,muD,hull);
                    int_subtract(hull,k1,hull);
                    ok = 1;
                    break;
                case 4: /* parameter muD */
                    int_add(DB,MB,K);
                    int_multiply(K2,K,K);
                    int_add(F,K,hull);
                    
                    int_multiply(hull,VB,hull);
                    
                    int_subtract(hull,Itprime,hull);
                    
                    int_subtract(It,DU,K);
                    int_subtract(K,MU,K);
                    int_multiply(K,K2,K);
                    int_subtract(hull,K,hull);
                    
                    int_divide(hull,DB,hull);
                    multiply(-1.0,hull);
                    
                    int_subtract(hull,kM,hull);
                    int_subtract(hull,k1,hull);
                    ok = 1;
                    break;
            }
            break;
        case 1:
            switch (i) {
                case 0:
                    int_multiply(kM,DB,hull);
                    int_multiply(muM,MB,K);
                    int_subtract(hull,K,hull);
                    
                    int_divide(hull,F,hull);
                    ok = 1;
                    break;
                case 3:
                    int_multiply(F,VB,hull);
                    int_multiply(muM,MB,K);
                    int_add(hull,K,hull);
                    int_divide(hull,DB,hull);                    
                    ok = 1;
                    break;
                case 5:
                    int_multiply(F,VB,hull);
                    int_multiply(kM,DB,K);
                    int_subtract(hull,K,hull);
                    int_divide(hull,MB,hull);
                    multiply(-1.0,hull);
                    ok = 1;
                    break;
            }
            break;
        case 2:
            switch (i) {
                case 4:
                    int_divide(F,DB,hull);
                    ok = 1;
                    break;
            }
            break;
        case 3:
            switch (i) {
                case 5:
                    int_divide(F,MB,hull);
                    ok = 1;
                    break;
            }
            break;
    }
    
    /* now compare parameter to the interval hull */
    if (ok) ok = intersect(p[i],hull);
        
    return ok;
}


/* Do not edit below this line. */
/* This function initializes the model structure. */
void NAMEGEN(_model_init)(struct s_model *model, FILE *fptr) {
    
    char sep;
    int i,j;
    
    model->name = MODEL_NAME;
    model->DEdim = DEDIM;
    model->OFdim = OFDIM;
    model->AEdim = AEDIM;
    model->TEdim = AEDIM+DEDIM;
    model->nparms = NPARMS;
    model->p_extreme = (Range *) malloc(NPARMS*sizeof(Range));
    for (i=0; i<NPARMS; i++) {
        for (j=0; j<2; j++) model->p_extreme[i][j] = P_EXTREME[i][j];
    }
    model->y_extreme = (Range *) malloc(DEDIM*sizeof(Range));
    for (i=0; i<DEDIM; i++) {
        for (j=0; j<2; j++) model->y_extreme[i][j] = Y_EXTREME[i][j];
    }
    model->n_ispecs = NISPECS;
    model->ispecs = (int *) malloc(NISPECS*sizeof(int));
    if (readtoken(fptr,NISPECS,"model_ispecs","%d",model->ispecs,&sep)
	    != NISPECS) {
        fprintf(stderr,"parm_reduction: error reading model_ispecs.\n");
        exit(1);
    }
    model->n_fspecs = NFSPECS;
    model->fspecs = (double *) malloc(NFSPECS*sizeof(double));
    if (readtoken(fptr,NFSPECS,"model_fspecs","%lf",model->fspecs,&sep)
	    != NFSPECS) {
        fprintf(stderr,"parm_reduction: error reading model_fspecs.\n");
        exit(1);
    }
    model->determine_dependence = &NAMEGEN(_determine_dependence);
    model->vec_field = &NAMEGEN(_vec_field);
    model->interval_vec_field = &NAMEGEN(_interval_vec_field);
    model->invert_vec_field = &NAMEGEN(_invert_vec_field);
    model->parmdep = (int **) malloc((DEDIM+AEDIM)*sizeof(int *));
    model->parmdep_flag = (int *) calloc(DEDIM+AEDIM,sizeof(int));
    model->vardep = (int **) malloc((DEDIM+AEDIM)*sizeof(int *));
    model->ofdep = (int **) malloc((DEDIM+AEDIM)*sizeof(int *));
    model->varpartialdep = (int **) malloc((DEDIM+AEDIM)*sizeof(int *));
    model->ofpartialdep = (int **) malloc((DEDIM+AEDIM)*sizeof(int *));
    for (i=0; i<DEDIM+AEDIM; i++) {
        model->parmdep[i] = PARMDEP[i];
        model->vardep[i] = VARDEP[i];
        model->ofdep[i] = OFDEP[i];
        model->varpartialdep[i] = VARPARTIALDEP[i];
        model->ofpartialdep[i] = OFPARTIALDEP[i];
        for (j=0; j<NPARMS; j++) {
            if (model->parmdep[i][j]) { model->parmdep_flag[i] = 1; break; }
        }
    }
}
