/* This file contains a template for the model-specific information and code.
 * The ??? symbols should be replaced with appropriate code and then the file 
 * s_parm_reduction.h should be editted to add this model to the list in the
 * variable all_models. Then the program should be recompiled. */

/* First define a prefix for the function names which is unique to this model.
 * Put this in quotes in the definition of MODEL_NAME, and put it without quotes
 * in the definition of NAMEGEN. */
#define MODEL_NAME  "NLP"
#define NAMEGEN(A) NLP ## A

/* Next define values for the model: */
#define DEDIM  2   /* Enter the dimension of the Differential Equation. */
#define OFDIM  0   /* Enter the dimension of other functions. */
#define AEDIM  0   /* Enter the dimension of the Algebraic State Equations */
#define NPARMS  3  /* Enter the number of parameters in the DE. */
#define NISPECS  0  /* Enter number of integer specs for the model. */
#define NFSPECS  3  /* Enter number of floating point specs for the model. */

/* The user must code two external variables and two functions for the MODEL.
 * The two variables simply define whether each equation in the model involves 
 * each parameter or variable.
 * The two functions are:
 * 1) A dependence function which determines whether the vector field is monotonic
 *    with respect to the parameters.
 * 2) A function which computes f + ay where f is the vector field.
 */

/* Do not change the next 4 lines. */
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "parm_red_structs.h"
#include "readtoken.h"
/* Add any other include files or defines needed for the coded functions. */
#include <math.h>
#include "interval_math.h"
#define GRAV (model->fspecs[0])
#define BB (model->fspecs[1])
#define OMEGA (model->fspecs[2])

/* NLP: Forced, damped, nonlinear pendulum.
 * The vector field is
 * f = [ y[1]                                                    ]
 *     [ -GRAV*p[0]*sin(y[0]) - p[1]*y[1] + p[2]*BB*sin(OMEGA*t) ] 
 * where 
 *   GRAV is the the gravitational acceleration,
 *   BB is the amplitude of the forcing function,
 *   OMEGA is the frequency of the forcing function,
 *   p[0] is 1/L
 *   p[1] is a/m
 *   p[2] is 1/mL
 *   L is the length of the pendulum
 *   a is the damping coefficient 
 *   m is the mass.
 *   It is assumed grav, b, and omega are known precisely.
 *   */

/* Specify whether each equation involves the parameters and variables and other
 * functions.  1 indicates this equation depends on the parameter/variable, 
 * 0 indicates no. */
static int PARMDEP[DEDIM+AEDIM][NPARMS] = {
    {0, 0, 0},
    {1, 1, 1}
};

static int VARDEP[DEDIM+AEDIM][DEDIM] = {
    {0, 1},
    {1, 1}
};

static int VARPARTIALDEP[DEDIM+AEDIM][DEDIM+1] = {
    {0, 0, 0},
	{1, 1, 1},
};

static int OFDEP[DEDIM+AEDIM][OFDIM] = {
    {},
    {}
};

static int OFPARTIALDEP[DEDIM+AEDIM][OFDIM] = {
    {},
	{}
};

/* Specify any a priori extreme bounds on the parameters.  Use +-DBL_MAX to
 * represent +- infinity, if there are no bounds, and +-DBL_MIN for the numbers
 * closest to zero. */
static double P_EXTREME[NPARMS][2] = {
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX},
    {DBL_MIN, DBL_MAX}
};

static double Y_EXTREME[DEDIM][2] = {
    {-DBL_MAX, DBL_MAX},
    {-DBL_MAX, DBL_MAX}
};

void NAMEGEN(_determine_dependence)(int eqn, Range *p, Range t, Range *y, Range *u, Range partial,
                                   struct s_model *model, int *dep, struct s_discr *discr, int parm) {
    
    int i;
    Range z,t_omega;
    double temp;
    
    /* dep[i]=1 if equation eqn of the vector field is monotonically nondecreasing
     * with p[i] for each range of t, y, and u values.  -1 if it is monotonically 
     * nonincreasing, 0 if it does not depend on p[i], and -2 if it is not monotonic. 
     * Return a value of 0 if any of dep entries are -2, else return 1. */
    /* explain -3 setting */
    switch (eqn) {
        case 0:
            if (model->var == 0) { /* If y[0] (x) is unmeasured */
                partial[0] = 0.0;
                partial[1] = 0.0;
            }
            else if (model->var == 1) { /* If y[1] (y) is unmeasured */
                partial[0] = 1.0;
                partial[1] = 1.0;                
            }
            break;
        default: /* eqn == 1 */
            switch (parm) {
                case 0:
//                    dep[parm] = -2;

                    int_sine(y[0],z);
                    int_put(z,partial);
                    multiply(-GRAV,partial);
                    if (ISZERO_INTERIOR(z)) {
                        dep[parm] = -3;
//                        dep[parm] = -2;
                    }
                    else {
                        temp = z[0] + z[1];
                        if (temp > 0.0)
                            dep[parm] = -1;
                        else if (temp < 0.0)
                            dep[parm] = 1;
                        else
                            dep[parm] = 0;
                    }
                    break;
                case 1:
//                    dep[parm] = -2;

                    int_put(y[1],partial);
                    multiply(-1.0,partial);
                    if (ISZERO_INTERIOR(y[1])) {
                        dep[parm] = -3;
//                        dep[parm] = -2;
                    }
                    else {
                        temp = y[1][0] + y[1][1];
                        if (temp > 0.0)
                            dep[parm] = -1;
                        else if (temp < 0.0)
                            dep[parm] = 1;
                        else
                            dep[parm] = 0;
                    }
                    break;
                case 2:
//                    dep[parm] = -2;

                    for (i=0; i<2; i++) t_omega[i] = t[i]*OMEGA;
                    int_sine(t_omega,z);
                    int_put(z,partial);
                    multiply(BB,partial);
                    if (ISZERO_INTERIOR(z)) {
                        dep[parm] = -3;
//                        dep[parm] = -2;
                    }
                    else {
                        temp = z[0] + z[1];
                        if (temp > 0.0)
                            dep[parm] = 1;
                        else if (temp < 0.0)
                            dep[parm] = -1;
                        else
                            dep[parm] = 0;
                        break;
                    }
                    break;
                default: /* unmeasured variable */
                    if (model->var == 0) { /* If y[0] (x) is unmeasured */
                        int_cosine(y[0],partial);
                        int_multiply(p[0],partial,partial);
                        multiply(GRAV,partial);
                        multiply(-1.0,partial);
                    }
                    else { /* If y[1] (y) is unmeasured */
                        int_put(p[1],partial);
                        multiply(-1.0,partial);
                    } 
            }
    }

    return;
}

double NAMEGEN(_vec_field)(int cor, int eqn, double *p, double a, double t,
                           Range *y, Range *u, struct s_model *model) {
    /* This function returns either the lower or upper value (depending on whether
     * cor is 0 or 1) of the eqn'th equation of f + ay, where f is the vector
     * field.
     *
     * The parameters p are assumed nonnegative. */
        
    int j,opp_cor;
    Range z;
    double g,temp;
    
    if (cor == 0) opp_cor = 1;
    else opp_cor = 0;
        
    switch (eqn) {
        case 0:
            if (a >= 0.0) j = cor;
            else j = opp_cor;
            g = y[1][cor] + a*y[0][j];
            break;
        default: /* eqn == 1 */
            int_sine(y[0],z);
            g = -GRAV*p[0]*z[opp_cor];
            temp = -p[1] + a;
            if (temp >= 0.0) j = cor;
            else j = opp_cor;
            g += temp*y[1][j];
            temp = sin(OMEGA*t);
            if (temp >= 0.0) j = cor;
            else j = opp_cor;
            g += temp*p[2]*BB;
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
        
    int i;
    Range *z=model->window.z;
    Range g,w;
    
    switch (eqn) {
        case 0:
            /* no need to code, never non-monotonic */
            break;
        default: /* eqn == 1 */            
            for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) int_sine(y[0][i],z[i]);
            int_vfield_sum(z,w,discr);
            
            int_multiply(p[0],w,w);        
            multiply(-GRAV,w);
            int_put(w,g);
            
            int_vfield_sum(y[1],w,discr);
            int_multiply(p[1],w,w);
            multiply(-1.0,w);
            int_add(g,w,g);

            for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++){
                z[i][0] = sin(OMEGA*t[i]);  
                z[i][1] = z[i][0];
            }    
            int_vfield_sum(z,w,discr);
            int_multiply(p[2],w,w);
            multiply(BB,w);
            int_add(g,w,g);
            break;
    }
    
    return g[cor];
}

int NAMEGEN(_invert_vec_field)(int eqn, int i, Range F, struct s_model *model, Range *p, double t, Range *y, Range *u) {
    
    /* in this function, we are solving for parameter i in equation eqn.
     * If we start with f^s(p_i,**) = F (f^s is the vector field at the s-th time step, and F is the input to the function
     * then we must invert f^s, and find the expression p_i = (f^s)^(-1)(F,**) */ 
    
    int ok = 0;
    Range hull,K;
    
    switch (eqn) {
        case 0:
            int_put(F,hull);
            break;
        case 1:
            switch (i) {
                case 0:
                    int_multiply(p[1],y[1],K);
                    int_add(F,K,hull);
                    
                    int_put(p[2],K);
                    multiply(BB*sin(OMEGA*t),K);
                    int_subtract(hull,K,hull);
                    
                    int_sine(y[0],K);
                    multiply(-GRAV,K);
                    int_divide(hull,K,hull);
                    ok = 1;
                    break;
                case 1:
                    int_sine(y[0],K);
                    int_multiply(p[0],K,K);
                    multiply(GRAV,K);
                    int_add(F,K,hull);
                    
                    int_put(p[2],K);
                    multiply(BB*sin(OMEGA*t),K);
                    int_subtract(hull,K,hull);
                    
                    int_divide(hull,y[1],hull);
                    multiply(-1.0,hull);
                    ok = 1;
                    break;
                case 2:
                    int_sine(y[0],K);
                    int_multiply(p[0],K,K);
                    multiply(GRAV,K);
                    int_add(F,K,hull);
                    
                    int_multiply(p[1],y[1],K);
                    int_add(hull,K,hull);
                    
                    multiply(1/(BB*sin(OMEGA*t)),hull);
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
