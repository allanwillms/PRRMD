#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include "parm_red_structs.h"
#include "parm_red_discr_sets.h"
#include "parm_red_models.h"
#include "parm_reduction.h"
#include "read_data.h"
#include "readtoken.h"
#include "tasle.h"
#include "interval_math.h"
#include "discr_error_test.h"

#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define ISZERO_INTERIOR(a) (((a)[0] * (a)[1]) < 0.0 ? 1 : 0 )
#define db1 
//#define db2
//#define db3
//#define db5
//#define db_boxcut
//#define influence_split
//#define measured_variable_monotonicity
//#define db_hull_results
//#define db9
//#define tasleoutput
//#define unknown_output
//#define box_summary
//#define E_test
//#define output1    /* input: maxboxes - output: mu, p(c.o.m) (p-cols), distance from p(true), mu_hull, time  */
//#define output4    /* input: 1000*stepsize,s+1 - output: mu, p(c.o.m) (p-cols), distance from p(true), mu_hull, time */
//#define timing
//#define ppp_u_temp

//#define double_ctl_points
#define single_ctl_points
//--show-reachable=yes
//valgrind --leak-check=yes --dsymutil=yes --leak-check=full  ./parm_reduction sir2
// valgrind --tool=cachegrind --dsymutil=yes ./parm_reduction sir2
// cg_annotate cachegrind.out.90871 --auto=yes > profile_18.txt

// TO DO:   1. monotonicity of u
//          2. new splitting (of all unmeasured for option 2)

double *P_TRUE;

int main(int argc, char *argv[]) {
    
#ifdef timing
    time_t start = clock();
#endif
    char *USEAGE_TEXT[] = {
		"parm_reduction"
		"   Usage:  parm_reduction basefile",
		"      basefile: base file name.  The file basefile_in will be read as",
		"             input and the file basefile_out will be written as output.",
		" ",
		" parm_reduction reduces parameter boxes for an ODE model given some data.",
		" ",
		" Copyright 2010, Allan Willms,",
		" Dept. of Mathematics and Statistics,",
		" University of Guelph, Guelph, ON N1G 2W1, Canada"
    };
#define NUM_USEAGE_TEXT_LINES (sizeof USEAGE_TEXT / sizeof(char *))
	
    struct s_data data;
    struct s_parms parms;
    struct s_model model;
    struct s_discr_set_list discr_set_list;
    char fname[FILENAME_MAX];
    int i;
    
#ifdef output4
//    if (atof(argv[2])/1000*atoi(argv[3]) > 20.0) {
//        printf("%d\t%f\t0.01\t\t200\t\t100.005\t\t0.000001\t\t1000\t\t500.0000005\t\t1\t\t0\n",
//               atoi(argv[3]),atof(argv[2])/1000);
//        exit(1);
//    }
#endif
    
    /* initialize random generator */
    struct timeval t1;
    gettimeofday(&t1, NULL);
    srand(t1.tv_usec * t1.tv_sec);
//    srand(1323000);
    if (argc < 2 || strncmp(argv[1],"-h",2) == 0) {
		for (i=0; i<NUM_USEAGE_TEXT_LINES; i++) 
			fprintf(stdout,"%s\n",USEAGE_TEXT[i]);
		exit(0);
    }
    sprintf(fname,"%s_in",argv[1]);
#ifdef db1
	printf("Calling aquire_input...\n");
#endif
    aquire_input(&model,&discr_set_list,&data,&parms,fname,argv);
#ifdef db1
	printf("Calling reduce_all...\n");
	fflush(NULL);
#endif
    reduce_all(&model,&discr_set_list,&data,&parms);
#ifdef db1
	printf("Calling write_output...\n");
#endif
    sprintf(fname,"%s_out",argv[1]);
    write_output(&model,&discr_set_list,&data,&parms,fname);
    
#ifdef timing
    time_t stop = clock();
    printf("%f",difftime(stop,start)/CLOCKS_PER_SEC);
//    printf("Finished in about %f seconds. \n", difftime(stop, start)/CLOCKS_PER_SEC);
#endif
#if defined(output1) || defined(output4)
    printf("\n");
#endif
    exit(0);
}

void aquire_input(struct s_model *model, struct s_discr_set_list *discr_set_list,
				  struct s_data *data, struct s_parms *parms, char *fname, char *argv[]) {

    char *colspec[4],**filename,name[100],name2[100],token[30],sep;
    FILE *fptr;
    int base_winsize,col,file,flag,i,j,k,m_lo,m_hi,n,p,nfiles,ncols,*npts,ptruespec,*temp_use_var;
    double cpt_spread,***filedata,max_slope,temp,tempa[4];
    
    /* Open the input file. */
    if ((fptr = fopen(fname,"r")) == NULL) {
		fprintf(stderr,"parm_reduction: error opening %s.\n",fname); exit(1);
    }
	
    /* Get the discretization name. */
    if (readtoken(fptr,1,"discretization_set_name","%99s",name,&sep) != 1) {
		fprintf(stderr,"parm_reduction: error, could not find discretization_set_name.\n"); exit(1);
    }
    
    /* How many discretizations are needed? (more can be added manually if required)
     * the array index refers to -> 0: user_defined, 1: base, 2, 'alpha', 3 'beta', 4, algebraic */
    discr_set_list->n_discr_types = 4; 
    if (((discr_set_list->discr_set = (struct s_discr_set*) 
          malloc(discr_set_list->n_discr_types*sizeof(struct s_discr_set))) == NULL)){
        fprintf(stderr,"parm_reduction: error allocating memory for "
                " discr_set_list->discr_set\n"); exit(1);
    }
    /* Initialize the discretization structure inputted by the user. */
    i = 0;
    while (i<NUM_DSETS && strcmp(all_discr_sets[i].name,name) != 0) i++;
    if (i == NUM_DSETS) {
		fprintf(stderr,"parm_reduction: specified discretization set "
                "%s not valid.\n",name); exit(1);
    }
    /* set flag to zero to indicate default discretization */
#ifdef output4
    discr_set_list->discr_set[0].input_size = atoi(argv[2]);
#endif
    for (k=0; k<2; k++) discr_set_list->discr_set[0].flag[k] = 0; 
    (*all_discr_sets[i].init_fcn)(&discr_set_list->discr_set[0],fptr);
    discr_set_list->discr_set[0].flag[0] = 0;
    if(strchr(name,'C') != NULL) 
        base_winsize = discr_set_list->discr_set[0].ispecs[0] + discr_set_list->discr_set[0].ispecs[1] + 1;
    else
        base_winsize = discr_set_list->discr_set[0].ispecs[0] + 1;

    /* Get the discretization name for unmeasured subwindows */
    if (readtoken(fptr,1,"discretization_unmeasured_name","%99s",name2,&sep) != 1) 
        strcpy(name2,name);
    /* Initialize the discretization structure inputted by the user. */
    i = 0;
    while (i<NUM_DSETS && strcmp(all_discr_sets[i].name,name2) != 0) i++;
    if (i == NUM_DSETS) {
		fprintf(stderr,"parm_reduction: specified discretization set "
                "%s not valid.\n",name2); exit(1);
    }
    discr_set_list->discr_set[1].flag[0] = base_winsize-1;
    /* if the chosen discretization is a cumulative discretization, we will need to make enough
     * space in memory when initializing the non-chosen discretizations */
    if(strchr(name,'C') == NULL) discr_set_list->discr_set[1].flag[1] = 0;
    else discr_set_list->discr_set[1].flag[1] = discr_set_list->discr_set[0].ispecs[1];
    (*all_discr_sets[i].init_fcn)(&discr_set_list->discr_set[1],fptr); 

    /* Read in the sub-window parameters (maximum and minimum window lengths) */
    i = readtoken(fptr,2,"discretization_unmeasured_minmax","%d",&discr_set_list->sub_win_var,&sep);
    /* ensure that all entries are appropriate and fit within the full window size */
    j = 0; /* used as a flag */
    if (i == -1) {
        fprintf(stderr,"parm_reduction: Failed to located 'discretization_unmeasured_minmax'"
                " in input file.\n"); exit(1);
    }
    if (i == 0) /* if there are no sub-window parameters */
        for (j=0; j<2; j++) discr_set_list->sub_win_var[j] = -1;
    else if (i == 1)
        discr_set_list->sub_win_var[1] = discr_set_list->sub_win_var[0];
    
    if ((discr_set_list->sub_win_var[0] > discr_set_list->sub_win_var[1]) ||
        (discr_set_list->sub_win_var[0] < 2 && j != 2) || 
        (discr_set_list->sub_win_var[1] > base_winsize)) {
        fprintf(stderr,"parm_reduction: error in 'discretization_unmeasured_minmax':\n\t\t(1) "
                "subwindow minimum must be less than subwindow maximum.\n\t\t(2) subwindow minimum "
                "must be at least 2.\n\t\t(3) subwindow maximum must be no greater than %d.\n",
                base_winsize); exit(1);
    }
    
    /* Get the discretization name for hull consistency subwindows */
    if (readtoken(fptr,1,"discretization_hull_name","%99s",name2,&sep) != 1) 
        strcpy(name2,name);
    /* Initialize the discretization structure inputted by the user. */
    i = 0;
    while (i<NUM_DSETS && strcmp(all_discr_sets[i].name,name2) != 0) i++;
    if (i == NUM_DSETS) {
		fprintf(stderr,"parm_reduction: specified discretization set "
                "%s not valid.\n",name2); exit(1);
    }
    
    discr_set_list->discr_set[2].flag[0] = base_winsize-1;
    /* if the chosen discretization is a cumulative discretization, we will need to make enough
     * space in memory when initializing the non-chosen discretizations */
    if(strchr(name,'C') == NULL)  discr_set_list->discr_set[2].flag[1] = 0;
    else  discr_set_list->discr_set[2].flag[1] = discr_set_list->discr_set[0].ispecs[1];
    (*all_discr_sets[i].init_fcn)(&discr_set_list->discr_set[2],fptr); 

    /* Read in window sizes */
    i = readtoken(fptr,2,"discretization_hull_minmax","%d",&discr_set_list->sub_win_parm,&sep);
    /* ensure that all entries are appropriate and fit within the full window size */
    j = 0; /* used as a flag */
    if (i == -1) {
        fprintf(stderr,"parm_reduction: Failed to located 'discretization_hull_minmax' "
                "in input file.\n"); exit(1);
    }
    else if (i == 0) /* if there are no sub-window parameters */
        for (j=0; j<2; j++) discr_set_list->sub_win_parm[j] = -1;
    else if (i == 1)
        discr_set_list->sub_win_parm[1] = discr_set_list->sub_win_parm[0];
    
    if ((discr_set_list->sub_win_parm[0] > discr_set_list->sub_win_parm[1]) ||
        (discr_set_list->sub_win_parm[0] < 2 && j!=2) || 
        (discr_set_list->sub_win_parm[1] > base_winsize)) {
        fprintf(stderr,"parm_reduction: error in 'discretization_hull_minmax':\n\t\t(1) "
                "subwindow minimum must be less than subwindow maximum.\n\t\t(2) subwindow "
                "minimum must be at least 2.\n\t\t(3) subwindow maximum must be no greater "
                "than %d.\n",base_winsize); exit(1);
    }

    /* set-up the blank discretization for any potential algebraic equations */
    discr_set_list->discr_set[3].n_discr = 1;
    if (((discr_set_list->discr_set[3].discr = (struct s_discr *) 
         malloc(discr_set_list->discr_set[3].n_discr*sizeof(struct s_discr))) == NULL) ||
        ((discr_set_list->discr_set[3].discr[0].stencil = (enum e_stencil *) 
          malloc(base_winsize*sizeof(enum e_stencil))) == NULL) ||
        ((discr_set_list->discr_set[3].discr[0].alpha = (double *) 
          calloc(base_winsize,sizeof(double))) == NULL) ||
        ((discr_set_list->discr_set[3].discr[0].beta = (double *) 
          malloc(base_winsize*sizeof(double))) == NULL)) {
        fprintf(stderr,"parm_reduction: error allocating space for discr_set[3] (AEdim)\n"); 
        exit(1);
    }
    for (i=0; i<base_winsize; i++){ 
        discr_set_list->discr_set[3].discr[0].stencil[i] = vec;
        discr_set_list->discr_set[3].discr[0].beta[i] = 1.0;
    }
    discr_set_list->discr_set[3].discr[0].winsize = base_winsize;

    /* Get the discretization error. */
    if (readtoken(fptr,1,"discretization_error","%lf",&model->discr_error,&sep) != 1) {
		fprintf(stderr,"parm_reduction: error, could not find discretization_error.\n"); exit(1);
    }

    /* Get the model name. */
    if (readtoken(fptr,1,"model_name","%99s",name,&sep) != 1) {
		fprintf(stderr,"parm_reduction: error, could not find model_name.\n"); exit(1);
    }
    
    /* Initialize the model structure. */
    i = 0;
    while (i<NUM_MODELS && strcmp(all_models[i].name,name) != 0) i++;
    if (i == NUM_MODELS) {
		fprintf(stderr,"parm_reduction: specified model %s not found.\n",name); exit(1);
    }
    (*all_models[i].init_fcn)(model,fptr);
        
    /* Read in the maximum number of parameter boxes. */
    if (readtoken(fptr,1,"maxboxes","%d",&parms->maxboxes,&sep) != 1) {
		fprintf(stderr,"parm_reduction: error, could not find maxboxes.\n"); exit(1);
    }
#ifdef output4
//    parms->maxboxes = atoi(argv[2]);
//    printf("parms->maxboxes updated to %d\n",parms->maxboxes);
#endif
	
    if (((parms->p_gw = (double*) malloc(model->nparms*sizeof(double))) == NULL) ||
		((parms->pbox = (struct s_pbox*) malloc(parms->maxboxes*sizeof(struct s_pbox))) == NULL) ||
		((parms->pbox[0].p = (Range*) malloc(model->nparms*sizeof(Range))) == NULL) ||
        ((parms->pbox[0].inf = (struct s_inf**) malloc(model->TEdim*sizeof(struct s_inf*))) == NULL)) {
            fprintf(stderr,"parm_reduction: error allocating memory for "
                    "parms->p_gw,pbox,pbox[0].p,pbox[0].inf\n"); exit(1);
    }
    for (i=0; i<model->TEdim; i++) {
        if (((parms->pbox[0].inf[i] = (struct s_inf*) 
              calloc((model->nparms+1),sizeof(struct s_inf))) == NULL)) {
            fprintf(stderr,"parm_reduction: error allocating memory for inf[i]"); exit(1);
        }
    }

    parms->nboxes = 1;
    parms->pbox[0].status = putative;

    /* Read in the parameter ranges and goal widths. */
    parms->orig_vol = 1.0;
    ptruespec = 0;
    P_TRUE = NULL;
    for (i=0; i<model->nparms; i++) {
		sprintf(token,"p[%d]",i);
		if ((j = readtoken(fptr,4,token,"%lf",tempa,&sep)) < 3) {
			fprintf(stderr,"parm_reduction: error, could not read range for %s.\n",token); exit(1);
		}
		for (k=0; k<2; k++) parms->pbox[0].p[i][k] = tempa[k];
		parms->p_gw[i] = tempa[2];
		if (j == 4) {
			if (!ptruespec && i>0) {
				fprintf(stderr,"parm_reduction: error, true parameter value for p[0] not given.\n"); 
                exit(1);
			}
			else if (!ptruespec) {
				ptruespec = 1;
				if ((P_TRUE = (double *) malloc(model->nparms*sizeof(double))) == NULL) {
					fprintf(stderr,"parm_reduction: error allocating memory for P_TRUE\n"); exit(1);
				}
			}
			P_TRUE[i] = tempa[3];
		}
		else if (ptruespec) {
			fprintf(stderr,"parm_reduction: error, only some true parameter values were given.\n"); 
            exit(1);
		}
		temp = SPAN(parms->pbox[0].p[i]);
		parms->orig_vol *= temp;
		if (temp < 0.0) {
			fprintf(stderr,"parm_reduction: error, parameter range for %s must be specified as "
                    "low high.\n",token); exit(1);
		}
		else if (parms->pbox[0].p[i][0] < model->p_extreme[i][0] ||
				 parms->pbox[0].p[i][1] > model->p_extreme[i][1]) {
			fprintf(stderr,"parm_reduction: error, parameter range for %s is outside a priori "
                    "bounds.\n",token); exit(1);
		}
    }
    parms->pbox[0].volume = parms->orig_vol;

    /* Read in the minimum fraction of the goal width for a successful reduction. */
    if (readtoken(fptr,1,"p_minreduction_fraction","%lf",&parms->minred,&sep) != 1) {
		fprintf(stderr,"parm_reduction: error, could not find p_minreduction_fraction.\n"); exit(1);
    }
    
    /* Read in the window specs */
    i = readtoken(fptr,4,"loop_specs","%lf",&discr_set_list->loop_specs,&sep);
    if (i == -1) {  
        fprintf(stderr,"parm_reduction: error, could not find window_specs.\n"); exit(1); 
    }
    else if (i == 0) {  
        discr_set_list->loop_specs[0] = 0.0;
        discr_set_list->loop_specs[1] = 1.0;
        discr_set_list->loop_specs[2] = 1.0;
        discr_set_list->loop_specs[3] = 10.0;
    }
    else if (i != 4) {
        fprintf(stderr,"parm_reduction: error, must be 4 specs entered for window_specs.\n"); exit(1); 
    }
    if(discr_set_list->loop_specs[0] >= 1){
        fprintf(stderr,"parm_reduction: window overlap spec (p) must be less than 1.\n"); exit(1);
    }

    /* Read in the subbox options */
    i = readtoken(fptr,2,"subbox_specs","%d",&model->window.subbox_specs,&sep);
    
    if (i == 0) 
        model->window.subbox_specs[2] = 0;
    else if (i == 2) 
        model->window.subbox_specs[2] = 1 + (sep == ':');
    else {  
        fprintf(stderr,"parm_reduction: error, could not find subbox_specs.\n"); exit(1);  
    }
    
#ifdef output4
    model->window.subbox_specs[2] = atoi(argv[4]);
    model->window.subbox_specs[0] = atoi(argv[5]);
    model->window.subbox_specs[1] = atoi(argv[6]);
//    printf("subbox specs = %d,%d,%d\n",model->window.subbox_specs[0],model->window.subbox_specs[1],model->window.subbox_specs[2]); exit(1);
#endif    
    
    if (i > 0 && model->window.subbox_specs[0] > base_winsize) {
        printf("%d > %d\n",model->window.subbox_specs[0],base_winsize);
        fprintf(stderr,"parm_reduction: error in subbox_specs, m must be less than or equal "
                "to discretization width.\n"); exit(1);
    }
    if (model->window.subbox_specs[2] == 2 && model->window.subbox_specs[0] <= 0) {
        fprintf(stderr,"parm_reduction: error in subbox_specs, m must be greater than zero.\n"); 
        exit(1);
    }
    /* set up the splitting ratios required for subboxes */
    if (((model->window.split_ratios = (double *) 
          malloc((model->window.subbox_specs[1]+1)*sizeof(double))) == NULL) ||
        ((model->window.time_split = (int *) malloc(base_winsize*sizeof(int))) == NULL)) {
        fprintf(stderr,"parm_reduction: error in process_equation allocating split_ratios memory.\n"); 
        exit(1);
    }      
    if (model->window.subbox_specs[2] == 0)
        for (i=0; i<2; i++) model->window.subbox_specs[i] = 0;
    else{
        /* below divides by powers of 2 - non-uniform */
//        temp = pow(2,(int)(model->window.subbox_specs[1]/2 + 1));
//        for (i = 0; i <= (model->window.subbox_specs[1]-1)/2; i++) {
//            model->window.split_ratios[i] = i/temp;
//            model->window.split_ratios[model->window.subbox_specs[1]-i] = 1 - i/temp;
//        }
        /* below divides evenly */
        temp = model->window.subbox_specs[1];
        for (i = 0; i <= (model->window.subbox_specs[1]-1)/2; i++) {
            model->window.split_ratios[i] = i/temp;
            model->window.split_ratios[model->window.subbox_specs[1]-i] = 1 - i/temp;
        }
        
        /* if the number of boxes is even, then set the middle element to 0.5 */
        if (model->window.subbox_specs[1]%2==0) 
            model->window.split_ratios[(model->window.subbox_specs[1]/2)] = 0.5;
        
        /* choosing which time steps to split */
        if (model->window.subbox_specs[2] == 2) {
            j=0;
            for (i=0; i<base_winsize; i++) {
                if (base_winsize*j/model->window.subbox_specs[0] == i) { 
                    /* select this one as a splitter */
                    model->window.time_split[i] = 1; j++;
                }
                else 
                    model->window.time_split[i] = 0;
            }
        }
    }

    /* Read in the number of data sets. */
    if (readtoken(fptr,1,"num_data_sets","%d",&data->nsets,&sep) != 1) {
		fprintf(stderr,"parm_reduction: error, could not read number of data files.\n"); exit(1);
    }
    
    if ((data->set = (struct s_set *) malloc(data->nsets*sizeof(struct s_set))) == NULL) {
		fprintf(stderr,"parm_reduction: error allocating memory for data->set\n"); exit(1);
    }
    for (i=0; i<4; i++) {
		if ((colspec[i] = (char *) malloc(40*sizeof(char))) == NULL) {
			fprintf(stderr,"parm_reduction: error allocating memory for colspec[%d]\n",i); exit(1);
		}
    }
    
    /* For any unmeasured variables found later, we declare an struct which will contain control
     * points for each unmeasured variables. This is first indexed by set, then by variable */
    if ((parms->pbox[0].vp = (struct s_yu **) malloc(data->nsets*sizeof(struct s_yu *))) == NULL) {
		fprintf(stderr,"parm_reduction: error allocating memory for parms->pbox[0].vp\n"); exit(1);
    }
    for (i=0; i<data->nsets; i++) {
        if ((parms->pbox[0].vp[i] = (struct s_yu *) 
             malloc(model->DEdim*sizeof(struct s_yu))) == NULL) {
            fprintf(stderr,"parm_reduction: error allocating memory for "
                    "parms->pbox[0].vp[i]\n"); exit(1);
        }
    }
    
    /* Set up a summary of unmeasured/measured variables over all data sets */
    if (((data->y_unmsrd_max = (int *) calloc(model->DEdim,sizeof(int))) == NULL) ||
        ((temp_use_var = (int *) malloc((model->DEdim+model->OFdim)*sizeof(int))) == NULL)) {
		fprintf(stderr,"parm_reduction: error allocating memory for data->y_unmsrd_max,temp_use_var\n"); exit(1);
    }
    
    /* For each data set: */
    data->num_y_unmsrd_max = 0;
    for (i=0; i<data->nsets; i++) {
		/* Allocate space for data structures for this data set. */
		if (((data->set[i].y = (struct s_yu*) malloc(model->DEdim*sizeof(struct s_yu))) == NULL) ||
			((data->set[i].u = (struct s_yu*) malloc(model->OFdim*sizeof(struct s_yu))) == NULL) ||
            ((data->set[i].unmsrd.y_unmsrd = (int*) calloc(model->DEdim,sizeof(int))) == NULL) ||
            ((data->set[i].unmsrd.y_unmsrd_index = (int*) calloc(model->DEdim,sizeof(int))) == NULL) ||
            ((data->set[i].unmsrd.y_unmsrd_order = (int*) malloc(model->DEdim*sizeof(int))) == NULL)) {
                fprintf(stderr,"parm_reduction: error allocating memory for "
                        "data->set[%d].y,u,y_unmsrd,y_unmsrd_order\n",i); exit(1);
        }
        data->set[i].unmsrd.num_y_unmsrd = 0; /* number of unmeasured variables in each set */
        data->set[i].unmsrd.base_winsize = base_winsize;
        for (j=0; j<model->DEdim+model->OFdim; j++) temp_use_var[j] = 0;
        
		/* Read in the maximum discretization step size for this data set. */
        j = readtoken(fptr,2,"step_size","%lf",&data->set[i].step_size,&sep);
		if (j == 1)
            data->set[i].step_size[1] = data->set[i].step_size[0];
        else if (j != 2) {
			fprintf(stderr,"parm_reduction: error, could not read maximum step "
					"size for data set %d.\n",i); exit(1);
        }
#ifdef output4
        data->set[i].step_size[0] = atof(argv[3])/1000;
        data->set[i].step_size[1] = atof(argv[3])/1000;
#endif

		/* Read in the number of data files for this data set. */
		if (readtoken(fptr,1,"num_files_in_set","%d",&nfiles,&sep) != 1) {
			fprintf(stderr,"parm_reduction: error, could not read number of "
					"files in data set %d.\n",i); exit(1);
		}
		/* Allocate space for the filename and filedata arrays. */ 
		if (((filename = (char **) malloc(nfiles*sizeof(char *))) == NULL) ||
			((filedata = (double ***) malloc(nfiles*sizeof(double **))) == NULL) ||
			((npts = (int *) malloc(nfiles*sizeof(int))) == NULL)) {
			fprintf(stderr,"parm_reduction: error allocating memory for filename array "
					"for data set %d.\n",i); exit(1);
		}
		for (file=0; file<nfiles; file++) {
			/* Allocate space for this filename. */
			if ((filename[file] = (char *) malloc(100*sizeof(char))) == NULL) {
				fprintf(stderr,"parm_reduction: error allocating memory for filename "
						"%d of data set %d.\n",file,i); exit(1);
			}
			/* Read in the data file name. */
			sprintf(token,"data_file_name_%d",file);
			if (readtoken(fptr,1,token,"%99s",filename[file],&sep) != 1) {
				fprintf(stderr,"parm_reduction: error, could not read %dth data "
						"file name in data set %d.\n",file,i); exit(1);
			}
			/* Read in the data from this file. */
			filedata[file] = read_data(filename[file],&ncols,&npts[file]);
		}
		/* Initialize the max time range. */
		data->set[i].t_range[0] = DBL_MAX;
		data->set[i].t_range[1] = -DBL_MAX;
        
		/* Read in the indications of how the variables relate to the read 
		 * in data and assign the data structure appropriately. */
		for (j=0; j<model->DEdim; j++) {
			sprintf(token,"y[%d]",j);
			if (readtoken(fptr,4,token,"%39s",colspec,&sep) < 4) {
				fprintf(stderr,"parm_reduction: error reading specification for y[%d] "
                        "for data set %d.\n",j,i); exit(1);
			}
			/* Read in the time input for this variable. */
			if ((col = parse_colspec(colspec[0],&file,&temp)) == -1) {
				fprintf(stderr,"parm_reduction: error parsing time specification for y[%d] "
                        "in data set %d.\n",j,i); exit(1);
			}
			data->set[i].y[j].npts = npts[file];
			data->set[i].y[j].input_t = filedata[file][col];
			data->set[i].t_range[0] = MIN(data->set[i].t_range[0],filedata[file][col][0]);
			data->set[i].t_range[1] = MAX(data->set[i].t_range[1],filedata[file][col][npts[file]-1]);
            
            /* Parse/apply the goal width. Negative goal width means y ranges will not be reduced. */
            data->set[i].y[j].gw = parse_gw(colspec[3],tempa);

			/* Read or calculate the data values for this variable. */
			if (sep == ':') {
                
                temp_use_var[j] = 1; /* indicating that this variable used tasle */
				/* use tasle to get bounds on data */
#ifdef db1
				printf("doing tasle fit for set %d, var %d\n",i,j);
				fflush(NULL);
#endif
				if ((col = parse_colspec(colspec[1],&file,&temp)) == -1 ||
					sscanf(colspec[2],"%lf",&(data->set[i].y[j].tmin)) != 1) {
					fprintf(stderr,"parm_reduction: error parsing specification "
							"for y[%d] in data set %d.\n",j,i); exit(1);
				}
				if (data->set[i].y[j].tmin < 0.0) data->set[i].y[j].tmin = 0.0;
				data->set[i].y[j].input = filedata[file][col];
				if (data->set[i].y[j].npts != npts[file]) {
					fprintf(stderr,"parm_reduction: error in specification for y[%d] "
							" of data set %d.\nThe time specification has %d points, "
							"but the data has %d points.\n",j,i,data->set[i].y[j].npts,
							npts[file]); exit(1);
				}                
				data->set[i].y[j].ncpt = tasle(data->set[i].y[j].npts,data->set[i].y[j].input_t,
                                               data->set[i].y[j].input,data->set[i].y[j].tmin,
                                               data->set[i].y[j].gw,&data->set[i].y[j].t,
                                               &data->set[i].y[j].cpt);
			}
			else if (sep == '=') {
#ifdef db1
				printf("doing automatic coarse banding for future reduction for set %d, var %d\n",i,j);
				fflush(NULL);
#endif                
                /* (2) Here is where we read in missing data sets and put in constant heights */
                data->set[i].unmsrd.y_unmsrd[j] = 1;                
                data->set[i].unmsrd.y_unmsrd_index[j] = model->nparms+data->set[i].unmsrd.num_y_unmsrd*base_winsize;
                data->set[i].unmsrd.y_unmsrd_order[data->set[i].unmsrd.num_y_unmsrd++] = j;
                data->y_unmsrd_max[j] = 1;
                
//                cpt_spread = 10*data->set[i].step_size[0]; /* NLP */
//                cpt_spread = 1*data->set[i].step_size[0]; /* SIR */
//                cpt_spread = 2*data->set[i].step_size[0]; /* HIV */
//                cpt_spread = data->set[i].step_size[0]; /* PHARM */
                cpt_spread = 1.0; /* none */
                
				data->set[i].y[j].tmin = -1.0;
#ifdef single_ctl_points /* number of control points (non-doubled) */
                data->set[i].y[j].ncpt = cpt_spread*SPAN(data->set[i].t_range)/data->set[i].step_size[0]+1;
#endif
#ifdef double_ctl_points /* number of control points (doubled) */
                data->set[i].y[j].ncpt = 2*cpt_spread*SPAN(data->set[i].t_range)/data->set[i].step_size[0];
                if (data->set[i].y[j].ncpt%2) data->set[i].y[j].ncpt++;
#endif
                parms->pbox[0].vp[i][j].ncpt = data->set[i].y[j].ncpt;
				if (((data->set[i].y[j].t = (double*) 
                      malloc((data->set[i].y[j].ncpt)*sizeof(double))) == NULL) ||
					((data->set[i].y[j].cpt = (Range*) 
                      malloc((data->set[i].y[j].ncpt)*sizeof(Range))) == NULL) ||
                    ((parms->pbox[0].vp[i][j].t = (double*) 
                      malloc((data->set[i].y[j].ncpt)*sizeof(double))) == NULL) ||
                    ((parms->pbox[0].vp[i][j].cpt = (Range*) 
                      malloc((data->set[i].y[j].ncpt)*sizeof(Range))) == NULL)) {
                        fprintf(stderr,"parm_reduction: error allocating memory for "
                                "data->set[%d].y[%d].t/cpt, parms->pbox[0].sp\n",i,j); exit(1);
                }
#ifdef single_ctl_points
                /* set the t-values (non-doubled) */        
				for (n=0; n<data->set[i].y[j].ncpt; n++){
                    data->set[i].y[j].t[n] = data->set[i].t_range[0] + 
                                             n*data->set[i].step_size[0]/cpt_spread;
                    parms->pbox[0].vp[i][j].t[n] = data->set[i].y[j].t[n];
                }
                n = data->set[i].y[j].ncpt-1;
                data->set[i].y[j].t[n] = data->set[i].t_range[1];
                parms->pbox[0].vp[i][j].t[n] = data->set[i].y[j].t[n];
#endif 
#ifdef double_ctl_points
                /* set the t-values (doubled) */
                data->set[i].y[j].t[0] = data->set[i].t_range[0];
                parms->pbox[0].vp[i][j].t[0] = data->set[i].y[j].t[0];
                for (n=1; n<data->set[i].y[j].ncpt/2; n++){
                    data->set[i].y[j].t[2*n-1] = data->set[i].t_range[0]+ 
                                                 n*data->set[i].step_size[0]/cpt_spread;
                    parms->pbox[0].vp[i][j].t[2*n-1] = data->set[i].y[j].t[2*n-1];
                    data->set[i].y[j].t[2*n] = data->set[i].t_range[0] + 
                                                 n*data->set[i].step_size[0]/cpt_spread;
                    parms->pbox[0].vp[i][j].t[2*n] = data->set[i].y[j].t[2*n];
                }
                n = data->set[i].y[j].ncpt-1;
                data->set[i].y[j].t[n] = data->set[i].t_range[1];
                parms->pbox[0].vp[i][j].t[n] = data->set[i].y[j].t[n];
#endif
                
				for (k=0; k<2; k++) {
					col = parse_colspec(colspec[k+1],&file,&temp);
					if (data->set[i].y[j].npts != npts[file]) {
						fprintf(stderr,"parm_reduction: error in specification for y[%d][%d] " 
                                " set %d.\nThe time specification has %d points, but the data "
                                "has %d points.\n",j,k,i,data->set[i].y[j].npts,npts[file]); exit(1);
					}  
					if (col >= 0) {
                        /* this means a column was inputted, so we interpolate to get control points */                        
                        m_lo = 0; m_hi = 1;
						for (n=0; n<data->set[i].y[j].ncpt; n++){
                            if (m_hi < data->set[i].y[j].npts-1 && 
                                data->set[i].y[j].input_t[m_hi] < data->set[i].y[j].t[n]) {
                                m_lo++; m_hi++;
                            }
                            data->set[i].y[j].cpt[n][k] = filedata[file][col][m_lo] + 
                                                    (filedata[file][col][m_hi]-filedata[file][col][m_lo])/
                                                    (data->set[i].y[j].input_t[m_hi] - 
                                                     data->set[i].y[j].input_t[m_lo])*
                                                    (data->set[i].y[j].t[n]-data->set[i].y[j].input_t[m_lo]);
                            parms->pbox[0].vp[i][j].cpt[n][k] = data->set[i].y[j].cpt[n][k];
                        }   
                    }
					else {
						for (n=0; n<data->set[i].y[j].ncpt; n++){
							data->set[i].y[j].cpt[n][k] = temp;
                            parms->pbox[0].vp[i][j].cpt[n][k] = temp;
                        }
                    }
				}
			}
            else if (sep == ';') {
#ifdef db1
				printf("doing manual banding for future reduction for set %d, var %d\n",i,j);
				fflush(NULL);
#endif          
                /* (1) this is when we want to do our own manual banding */
                data->set[i].y[j].tmin = -1.0;
				data->set[i].y[j].ncpt = data->set[i].y[j].npts;
				if (((data->set[i].y[j].t = (double *) 
                      malloc(data->set[i].y[j].ncpt*sizeof(double))) == NULL) ||
					((data->set[i].y[j].cpt = (Range *) 
                      malloc(data->set[i].y[j].ncpt*sizeof(Range))) == NULL)) {
                        fprintf(stderr,"parm_reduction: error allocating memory for "
                                "data->set[%d].y[%d].t/cpt\n",i,j); exit(1);
                    }

				for (n=0; n<data->set[i].y[j].ncpt; n++)
					data->set[i].y[j].t[n] = data->set[i].y[j].input_t[n];
				for (k=0; k<2; k++) {
					col = parse_colspec(colspec[k+1],&file,&temp);
					if (data->set[i].y[j].npts != npts[file]) {
						fprintf(stderr,"parm_reduction: error in specification for y[%d][%d] "
								" of data set %d.\nThe time specification has %d points, "
								"but the data has %d points.\n",j,k,i,data->set[i].y[j].npts,
								npts[file]); exit(1);
					}
					if (col >= 0)
						for (n=0; n<data->set[i].y[j].ncpt; n++)
							data->set[i].y[j].cpt[n][k] = filedata[file][col][n] + temp;
					else
						for (n=0; n<data->set[i].y[j].ncpt; n++)
							data->set[i].y[j].cpt[n][k] = temp;
				}
                
            }
            
			/* Ensure ranges are valid. */
			for (n=0; n<data->set[i].y[j].ncpt; n++) {
				if (SPAN(data->set[i].y[j].cpt[n]) < 0.0) {
					fprintf(stderr,"parm_reduction: error in range of y[%d] at time point "
                            "%d (t=%f).\n",j,n,data->set[i].y[j].t[n]);
                    printf("[%f,%f]\n",data->set[i].y[j].cpt[n][0],data->set[i].y[j].cpt[n][1]); 
                    exit(1);
				}
                
                if (data->set[i].y[j].cpt[n][0] < model->y_extreme[j][0]) {
                    if (data->set[i].y[j].cpt[n][1] > model->y_extreme[j][1]) {
                        printf("ERROR: calculated band does not fall in a priori variable "
                               "ranges.\n"); exit(1);
                    }
                    data->set[i].y[j].cpt[n][0] = model->y_extreme[j][0];
                }
                if (data->set[i].y[j].cpt[n][1] > model->y_extreme[j][1]) {
                    if (data->set[i].y[j].cpt[n][0] < model->y_extreme[j][0]) {
                        printf("ERROR: calculated band does not fall in a priori variable "
                               "ranges.\n"); exit(1);
                    }
                    data->set[i].y[j].cpt[n][1] = model->y_extreme[j][1];
                }
			}
            /* Determine the maximum size needed for extended pboxes */
            data->num_y_unmsrd_max = MAX(data->num_y_unmsrd_max,data->set[i].unmsrd.num_y_unmsrd);
    
#ifdef tasleoutput
//            double average_height = 0.0;
			int pp;
			for(pp=0; pp<data->set[i].y[j].ncpt; pp++){
				fprintf(stderr,"%f\t%f\t%f\t%f\n",data->set[i].y[j].t[pp],data->set[i].y[j].cpt[pp][0],
                        data->set[i].y[j].t[pp],data->set[i].y[j].cpt[pp][1]);
//                average_height += SPAN(data->set[i].y[j].cpt[pp]);
			} 
//            printf("AVERAGE HEIGHT = %f\n",average_height/data->set[i].y[j].ncpt);
#endif			
		}
        
        /* Now that we have looped through all sets, extend the pbox->p to include unmeasured variables */
        Range *more_p;
        if ((more_p = (Range*) realloc(parms->pbox[0].p,(model->nparms+data->num_y_unmsrd_max*
                                       base_winsize)*sizeof(Range))) == NULL){
            fprintf(stderr,"parm_reduction: error reallocating memory for more_p\n"); exit(1);
        }
        parms->pbox[0].p = more_p;
                
		/* Do the same for the other functions, except there are no goal widths. */
		for (j=0; j<model->OFdim; j++) {            

			sprintf(token,"u[%d]",j);
			if (readtoken(fptr,4,token,"%39s",colspec,&sep) < 3) {
				fprintf(stderr,"parm_reduction: error reading specification "
						"for u[%d] for data set %d.\n",j,i); exit(1);
			}
			/* Read in the time input for this variable. */
			if ((col = parse_colspec(colspec[0],&file,&temp)) == -1) {
				fprintf(stderr,"parm_reduction: error parsing time specification "
						"for u[%d] in data set %d.\n",j,i); exit(1);
			}
			data->set[i].u[j].npts = npts[file];
			data->set[i].u[j].input_t = filedata[file][col];
			data->set[i].t_range[0] = MIN(data->set[i].t_range[0],filedata[file][col][0]);
			data->set[i].t_range[1] = MAX(data->set[i].t_range[1],filedata[file][col][npts[file]-1]);
            
            /* Parse and apply the goal width.  */
            data->set[i].u[j].gw = parse_gw(colspec[3],tempa);

			/* Read or calculate the data values for this variable. */
			if (sep == ':') {
                temp_use_var[j] = 1; /* indicating that this variable used tasle bounding */

				/* use tasle to get bounds on data */
#ifdef db1
				printf("doing tasle fit for set = %d, function = %d\n",i,j); fflush(NULL);
#endif
				if ((col = parse_colspec(colspec[1],&file,&temp)) == -1 ||
					sscanf(colspec[2],"%lf",&(data->set[i].u[j].tmin)) != 1) {
					fprintf(stderr,"parm_reduction: error parsing specification "
							"for u[%d] in data set %d.\n",j,i); exit(1);
				}
				if (data->set[i].u[j].tmin < 0.0) data->set[i].u[j].tmin = 0.0;
				data->set[i].u[j].input = filedata[file][col];
				if (data->set[i].u[j].npts != npts[file]) {
					fprintf(stderr,"parm_reduction: error in specification for u[%d] of data set "
                            "%d.\nThe time specification has %d points, but the data has %d points.\n",
                            j,i,data->set[i].u[j].npts,npts[file]); exit(1);
				}
				data->set[i].u[j].ncpt = tasle(data->set[i].u[j].npts,data->set[i].u[j].input_t,
                                               data->set[i].u[j].input,data->set[i].u[j].tmin,
                                               data->set[i].u[j].gw,&data->set[i].u[j].t,
                                               &data->set[i].u[j].cpt);
#ifdef tasleoutput
                int pp;
                for(pp=0; pp<data->set[0].u[j].ncpt; pp++){
                    fprintf(stderr,"%f\t%f\t%f\t%f\n",data->set[0].u[j].t[pp],
                            data->set[0].u[j].cpt[pp][0],data->set[0].u[j].t[pp],
                            data->set[0].u[j].cpt[pp][1]);
                } 
#endif		
			}
			else {
				data->set[i].u[j].tmin = -1.0;
				data->set[i].u[j].ncpt = data->set[i].u[j].npts;
				if (((data->set[i].u[j].t = (double *) 
                      malloc(data->set[i].u[j].ncpt*sizeof(double))) == NULL) ||
					((data->set[i].u[j].cpt = (Range *) 
                      malloc(data->set[i].u[j].ncpt*sizeof(Range))) == NULL)) {
                        fprintf(stderr,"parm_reduction: error allocating memory for "
                                "data->set[%d].u[%d].t/cpt\n",i,j); exit(1);
                    }
				for (n=0; n<data->set[i].u[j].ncpt; n++)
					data->set[i].u[j].t[n] = data->set[i].u[j].input_t[n];
				for (k=0; k<2; k++) {
					col = parse_colspec(colspec[k+1],&file,&temp);
					if (data->set[i].u[j].npts != npts[file]) {
						fprintf(stderr,"parm_reduction: error in specification for u[%d][%d] "
								" of data set %d.\nThe time specification has %d points, "
								"but the data has %d points.\n",j,k,i,data->set[i].u[j].npts,
								npts[file]); exit(1);
					}
					if (col >= 0)
						for (n=0; n<data->set[i].u[j].ncpt; n++)
							data->set[i].u[j].cpt[n][k] = filedata[file][col][n] + temp;
					else
						for (n=0; n<data->set[i].u[j].ncpt; n++)
							data->set[i].u[j].cpt[n][k] = temp;
				}
			}
			/* Ensure ranges are valid. */
			for (n=0; n<data->set[i].u[j].ncpt; n++) {
				if (SPAN(data->set[i].u[j].cpt[n]) < 0.0) {
					fprintf(stderr,"parm_reduction: error in range of u[%d] "
							"at time point %d.\n",j,n); exit(1);
				}
#ifdef ppp_u_temp
                if (data->set[i].u[j].cpt[n][0] <= 0.0) data->set[i].u[j].cpt[n][0] = 0.0;
#endif
			}
		}
        
		/* Set the max t_range and extend any data so that all t ranges are equal */
		for (j=0; j<model->DEdim; j++) {
			if (data->set[i].t_range[0] < data->set[i].y[j].t[0]) {
//				fprintf(stdout,"parm_reduction: warning, smallest time for y[%d] in data "
//						"set %d is %lf, which is larger than \nthat of other variables: %lf.\n"
//						"y[%d][0] will be extended linearly.\n",j,i,data->set[i].y[j].t[0],
//						data->set[i].t_range[0],j);
				for (k=0; k<2; k++) {
					data->set[i].y[j].cpt[0][k] = data->set[i].y[j].cpt[0][k] +
					(data->set[i].t_range[0] - data->set[i].y[j].t[0])*
					(data->set[i].y[j].cpt[1][k] - data->set[i].y[j].cpt[0][k])/
					(data->set[i].y[j].t[1] - data->set[i].y[j].t[0]);
				}
				data->set[i].y[j].t[0] = data->set[i].t_range[0];
                /* reset if extension violates variable bounds */
                data->set[i].y[j].cpt[0][0] = MAX(data->set[i].y[j].cpt[0][0],model->y_extreme[j][0]);
                data->set[i].y[j].cpt[0][1] = MAX(data->set[i].y[j].cpt[0][1],
                                                  model->y_extreme[j][0]+0.5*data->set[i].y[j].gw);
                
                data->set[i].y[j].cpt[0][1] = MIN(data->set[i].y[j].cpt[0][1],model->y_extreme[j][1]);
                data->set[i].y[j].cpt[0][0] = MIN(data->set[i].y[j].cpt[0][0],
                                                  model->y_extreme[j][1]-0.5*data->set[i].y[j].gw);
			}
			n = data->set[i].y[j].ncpt - 1;
			if (data->set[i].t_range[1] > data->set[i].y[j].t[n]) {
//				fprintf(stdout,"parm_reduction: warning, largest time for y[%d] in data "
//						"set %d is %lf, which is smaller than \nthat of other variables: %lf.\n"
//						"y[%d][%d] will be extended linearly.\n",j,i,data->set[i].y[j].t[n],
//						data->set[i].t_range[1],j,n);
				for (k=0; k<2; k++) {
					data->set[i].y[j].cpt[n][k] = data->set[i].y[j].cpt[n][k] +
					(data->set[i].t_range[1] - data->set[i].y[j].t[n])*
					(data->set[i].y[j].cpt[n-1][k] - data->set[i].y[j].cpt[n][k])/
					(data->set[i].y[j].t[n-1] - data->set[i].y[j].t[n]);
				}
				data->set[i].y[j].t[n] = data->set[i].t_range[1];
                /* reset if extension violates variable bounds */
                data->set[i].y[j].cpt[n][0] = MAX(data->set[i].y[j].cpt[n][0],model->y_extreme[j][0]);
                data->set[i].y[j].cpt[n][1] = MAX(data->set[i].y[j].cpt[n][1],
                                                  model->y_extreme[j][0]+0.5*data->set[i].y[j].gw);

                data->set[i].y[j].cpt[n][1] = MIN(data->set[i].y[j].cpt[n][1],model->y_extreme[j][1]);
                data->set[i].y[j].cpt[n][0] = MIN(data->set[i].y[j].cpt[n][0],
                                                  model->y_extreme[j][1]-0.5*data->set[i].y[j].gw);
            }
		}
        
        /* Allocate space for array for storing steepness */
        p = (int) SPAN(data->set[i].t_range)/data->set[i].step_size[0];
        if (((data->set[i].change = (enum e_change*) calloc(p+2,sizeof(enum e_change))) == NULL)) {
            fprintf(stderr,"parm_reduction: error allocating memory for data->set[i].change\n"); exit(1);
        }
        /* For all bands that were constructed using tasle */
        flag = 0;
        for (j=0; j<model->DEdim; j++) {
            if (temp_use_var[j] == 1) { /* this variable was constructed using tasle */
                flag = 1; max_slope = 0.0;
                for (n=0; n<data->set[i].y[j].ncpt-1; n++) {
                    max_slope = MAX(max_slope,fabs((data->set[i].y[j].cpt[n+1][1]-data->set[i].y[j].cpt[n][1])/
                                                   (data->set[i].y[j].t[n+1]-data->set[i].y[j].t[n])));
                    max_slope = MAX(max_slope,fabs((data->set[i].y[j].cpt[n+1][0]-data->set[i].y[j].cpt[n][0])/
                                                   (data->set[i].y[j].t[n+1]-data->set[i].y[j].t[n])));
                }
                /* Now go through each point and classify according to steepness */               
//                p = 0;
//                for (n=0; n<data->set[i].y[j].ncpt-1; n++) {
//                    if (data->set[i].t_range[0]+p*data->set[i].step_size[0] < data->set[i].y[j].t[n+1]) {
//                        temp = MAX(fabs((data->set[i].y[j].cpt[n+1][1]-data->set[i].y[j].cpt[n][1])/
//                                        (data->set[i].y[j].t[n+1]-data->set[i].y[j].t[n])),
//                                   fabs((data->set[i].y[j].cpt[n+1][0]-data->set[i].y[j].cpt[n][0])/
//                                        (data->set[i].y[j].t[n+1]-data->set[i].y[j].t[n])));
//                        if (temp > 0.25*max_slope) { /* set the new value */
//                            if (temp <= 0.5*max_slope) data->set[i].change[p] = MAX(1,data->set[i].change[p]);
//                            else if (temp <= 0.75*max_slope) data->set[i].change[p] = MAX(2,data->set[i].change[p]);
//                            else data->set[i].change[p] = MAX(3,data->set[i].change[p]);
//                        }
//                        while (data->set[i].t_range[0]+(++p)*data->set[i].step_size[0] < data->set[i].y[j].t[n+1])
//                            data->set[i].change[p] = MAX(data->set[i].change[p],data->set[i].change[p-1]);
//                    }
//                }
            }
        }
        /* If no variables were banded using tasle, then manually set the steepness */
        if (flag == 0) for (n=0; n<p; n++) data->set[i].change[p] = 2;       
    }
    free(temp_use_var);
    
    /* For all equations, check to see if they are dependent on an unmeasured variable
     * if so, then adjust the flag value for use in process_window */
    for (j=0; j<model->TEdim; j++) {
        for (i=0; i<model->DEdim; i++) {
            if (data->y_unmsrd_max[i] && (j==i || model->vardep[j][i])) {
                model->parmdep_flag[j] = 2; break;
            }
        }
    }
        
    /* We do not free filedata[file][j] because these addresses were copied to other variables */
    if (data->nsets > 0) {
		free(filedata);
		for (i=0; i<nfiles; i++) free(filename[i]);
		free(filename);
		free(npts);
    }
    for (i=0; i<4; i++) free(colspec[i]);
    /* Read in the minimum fraction of the goal width for a successful reduction. */
    if (readtoken(fptr,1,"y_minreduction_fraction","%lf",&data->minred,&sep) != 1) {
		fprintf(stderr,"parm_reduction: error, could not find y_minreduction_fraction.\n"); exit(1);
    }
    /* Read in the output variables option. */
    if (readtoken(fptr,1,"write_box_vars","%s",name,&sep) != 1) {
		fprintf(stderr,"parm_reduction: error, could not find write_box_vars.\n"); exit(1);
    }
    
    k = model->nparms+base_winsize*data->num_y_unmsrd_max;
    model->tparms = k;
    /* Initialize memory for each window to be used later in the program */
    if (((model->window.windep[0] = (int *) calloc(k,sizeof(int))) == NULL) ||
        ((model->window.windep[1] = (int *) calloc(k,sizeof(int))) == NULL) ||
        ((model->window.windep[2] = (int *) calloc(k,sizeof(int))) == NULL) ||
        ((model->window.windep[3] = (int *) calloc(k,sizeof(int))) == NULL) ||
        ((model->window.windep[4] = (int *) calloc(k,sizeof(int))) == NULL) ||
        ((model->window.windep[5] = (int *) calloc(k,sizeof(int))) == NULL) ||
        ((model->window.p_prev_mono = (Range *) malloc(model->nparms*sizeof(Range))) == NULL) ||
        ((model->window.windep_msrd[0] = (int **) malloc(model->DEdim*sizeof(int*))) == NULL) ||
        ((model->window.windep_msrd[1] = (int **) malloc(model->DEdim*sizeof(int*))) == NULL) ||
        ((model->window.windep_msrd[2] = (int **) malloc(model->DEdim*sizeof(int*))) == NULL) ||
        ((model->window.eqn_order = (int *) malloc(model->TEdim*sizeof(int))) == NULL) ||
        ((model->window.partial = (Range *) malloc(model->nparms*sizeof(Range))) == NULL) ||
        ((model->window.attempt = (intRange *) malloc(k*sizeof(intRange))) == NULL) ||
        ((model->window.p_hull = (Range *) malloc(k*sizeof(Range))) == NULL) ||
		((model->window.tuse = (double *) malloc(base_winsize*sizeof(double))) == NULL) ||
		((model->window.y = (Range **) malloc(base_winsize*sizeof(Range *))) == NULL) ||
		((model->window.u = (Range **) malloc(base_winsize*sizeof(Range *))) == NULL) ||
		((model->window.t_ind = (int *) calloc(model->DEdim,sizeof(int))) == NULL) ||
        ((model->window.yext = (Range *) malloc(model->DEdim*sizeof(Range))) == NULL) ||
		((model->window.uext = (Range *) malloc(model->OFdim*sizeof(Range))) == NULL) ||
        ((model->window.yext_sort = (intRange **) malloc(model->DEdim*sizeof(intRange*))) == NULL) ||
        ((model->window.uext_sort = (intRange **) malloc(model->OFdim*sizeof(intRange*))) == NULL) ||
        ((model->window.p_alt = (Range *) malloc(model->nparms*sizeof(Range))) == NULL) ||
        ((model->window.z = (Range *) malloc(base_winsize*sizeof(Range))) == NULL) ||
        ((model->window.y_alt = (Range **) malloc(model->DEdim*sizeof(Range*))) == NULL) ||
        ((model->window.y_alt1 = (Range *) malloc(model->DEdim*sizeof(Range))) == NULL) ||
        ((model->window.u_alt = (Range **) malloc(model->OFdim*sizeof(Range*))) == NULL) ||
        ((model->window.p_prev[0] = (Range *) malloc(k*sizeof(Range))) == NULL) ||
        ((model->window.p_prev[1] = (Range *) malloc(k*sizeof(Range))) == NULL) ||
        ((model->window.p_prev[2] = (Range *) malloc(k*sizeof(Range))) == NULL) ||
        ((model->window.p_ord0 = (double *) malloc(k*sizeof(double))) == NULL) ||
        ((model->window.p_ord1 = (double *) malloc(k*sizeof(double))) == NULL) ||
        ((model->window.p_ord2 = (double *) malloc(k*sizeof(double))) == NULL)) {
		fprintf(stderr,"parm_reduction: error in process_equation allocating memory for window.\n"); 
        exit(1);
    }
    for (j = 0; j < base_winsize; j++) {
		if (((model->window.y[j] = (Range *) malloc(model->DEdim*sizeof(Range))) == NULL) ||
			((model->window.u[j] = (Range *) malloc(model->OFdim*sizeof(Range))) == NULL)) {
			fprintf(stderr,"parm_reduction: error allocating memory for window."); exit(1);
		}
    }
    for(j = 0; j < model->DEdim; j++){
        if (((model->window.y_alt[j] = (Range *) malloc(base_winsize*sizeof(Range))) == NULL) ||
            ((model->window.yext_sort[j] = (intRange*) malloc(base_winsize*sizeof(intRange))) == NULL) ||
            ((model->window.windep_msrd[0][j] = (int*) malloc(base_winsize*sizeof(int))) == NULL) ||
            ((model->window.windep_msrd[1][j] = (int*) malloc(base_winsize*sizeof(int))) == NULL) ||
            ((model->window.windep_msrd[2][j] = (int*) malloc(base_winsize*sizeof(int))) == NULL)) {
            fprintf(stderr,"parm_reduction: error allocating memory for "
                    "model->(y_alt/yext_sort/windep_msrd[0,1,2])[j]\n"); exit(1);
        }
    }
    for(j = 0; j < model->OFdim; j++){
        if (((model->window.u_alt[j] = (Range *) malloc(base_winsize*sizeof(Range))) == NULL) ||
            ((model->window.uext_sort[j] = (intRange*) malloc(base_winsize*sizeof(intRange))) == NULL)) {
            fprintf(stderr,"parm_reduction: error allocating memory for model->(su/uext_sort)[j]\n"); 
            exit(1);
        }
    }
    if (data->num_y_unmsrd_max == 0) for (i=0; i<3; i++) model->window.subbox_specs[i] = 0;
    
    /* Set equation ordering */
    j = 0;
    for (i=0; i<model->AEdim; i++) model->window.eqn_order[(j++)] = model->DEdim + i;
    for (i=0; i<model->DEdim; i++)
        if (data->y_unmsrd_max[i]) model->window.eqn_order[(j++)] = i;
    for (i=0; i<model->DEdim; i++)
        if (!data->y_unmsrd_max[i]) model->window.eqn_order[(j++)] = i; 

    for (i=0; i<3; i++) name[i] = tolower(name[i]);
    if (strncmp(name,"yes",3) == 0) data->write_vars = 1;
    else data->write_vars = 0;
    if (fclose(fptr) == EOF) {
		fprintf(stderr,"parm_reduction: error closing %s.\n",fname); exit(1);
    }
    
    return;
}

int parse_colspec(char *expr, int *file, double *val) {
	
    /* This function parses expressions of the form
     *      [file%d][col%d][{+-}%lf]
     *      for example:   file2col3+0.1
     *  It returns the column number and puts the file number
     *  and the {+-}%lf value into *file and *val respectively.
     *  If [file%d] is not present then 0 is put in *file.
     *  If col%d is not present, then -1 is returned.
     *  If %f is not present, 0.0 is put into *val. */
    int col;
    char *oploc,*cloc,*floc;
	
    if ((floc = strstr(expr,"file")) == NULL) *file = 0;
    else if (sscanf(floc+4,"%d",file) != 1) *file = 0;
    if ((cloc = strstr(expr,"col")) == NULL) {
		col = -1;
		if (sscanf(expr,"%lf",val) != 1) *val = 0.0;
    }
    else {
		if (sscanf(cloc+3,"%d",&col) != 1) col = -1;
		if ((oploc = strpbrk(cloc+3,"+-")) == NULL) *val = 0.0;
		else {
			if (sscanf(oploc+1,"%lf",val) != 1) *val = 0.0;
			if (*oploc == '-') *val = -(*val);
		}
    }
    return col;
}

double parse_gw(char *expr, double mmm[3]) {
	
    /* This function parses expressions involving numbers, the four basic arithmetic
     * operators: +,-,*,/, and the words "min", "max", and "mean".  It computes the
     * expression left to right but does not obey order of operations rules and does
     * not support parentheses. */
    int i,len;
    double temp,val;
    char *mloc,*oploc[2];
    static char *ops = "+-*/";
	
    val = 0.0;
    oploc[0] = ops;  
    /* The above makes the first operator a + so, 0.0 + is effectively prepended 
     * to the expression. */
    len = (int) strlen(expr);
    i = 0;
    while (i<len) {
		/* Find the next operator. */
		if ((oploc[1] = strpbrk(expr+i,ops)) == NULL) oploc[1] = expr+len;
		/* Find the next location with "m" */
		if (((mloc = strchr(expr+i,'m')) != NULL) && (mloc < oploc[1])) {
			if (strncmp(mloc,"min",3) == 0)
				temp = mmm[0];
			else if (strncmp(mloc,"max",3) == 0)
				temp = mmm[1];
			else if (strncmp(mloc,"mean",4) == 0)
				temp = mmm[2];
			else
				temp = 0.0;  /* This is really a parse error, but we ignore. */
		}
		else if (oploc[1] > expr+i) {
			if (sscanf(expr+i,"%lf",&temp) != 1) temp = 0.0;
		}
		else
			temp = 0.0;
		
		/* Operate between the previous value and the value just read in. */
		switch (*oploc[0]) {
			case '+' : val += temp; break;
			case '-' : val -= temp; break;
			case '*' : val *= temp; break;
			case '/' : val /= temp; break;
		}
		oploc[0] = oploc[1];
		i = oploc[0] - expr + 1;
    }
    return val;
}

void write_output(struct s_model *model, struct s_discr_set_list *discr_set_list,
				  struct s_data *data, struct s_parms *parms, char *fname) {
	
    char datafname[100];
    int box,i,j,k,n,nvalid;
    FILE *fptr;
    double avgdist,avgarea,dist,totalvol,vol;
    	
    if ((fptr = fopen(fname,"w")) == NULL) {
		fprintf(stderr,"parm_reduction: error opening %s.\n",fname); exit(1);
    }
    fprintf(fptr,"DISCRETIZATION:\n");
    fprintf(fptr,"discretization_set_name = %s\ndiscretization_ispecs =",discr_set_list->discr_set[0].name);
    for (i=0; i<discr_set_list->discr_set[0].n_ispecs; i++) 
		fprintf(fptr," %d",discr_set_list->discr_set[0].ispecs[i]);
    fprintf(fptr,"\ndiscretization_fspecs =");
    for (i=0; i<discr_set_list->discr_set[0].n_fspecs; i++) 
		fprintf(fptr," %f",discr_set_list->discr_set[0].fspecs[i]);
    fprintf(fptr,"\ndiscretization.n_discr = %d\n",discr_set_list->discr_set[0].n_discr);
    for (i=0; i<discr_set_list->discr_set[0].n_discr; i++) {
		fprintf(fptr,"discretization.discr[%d]:\n",i);
		fprintf(fptr,"     discrtype = %d\n",discr_set_list->discr_set[0].discr[i].discrtype);
		fprintf(fptr,"     winsize = %d\n",discr_set_list->discr_set[0].discr[i].winsize);
		fprintf(fptr,"     stencil = ");
		for (j=0; j<discr_set_list->discr_set[0].discr[i].winsize; j++) 
			fprintf(fptr," %d",discr_set_list->discr_set[0].discr[i].stencil[j]);
		fprintf(fptr,"\n");
    }
    fprintf(fptr,"\nMODEL:\n");
    fprintf(fptr,"model_name = %s\nmodel_ispecs =",model->name);
    for (i=0; i<model->n_ispecs; i++) 
		fprintf(fptr," %d",model->ispecs[i]);
    fprintf(fptr,"\nmodel_fspecs =");
    for (i=0; i<model->n_fspecs; i++) 
		fprintf(fptr," %f",model->fspecs[i]);
    fprintf(fptr,"\nmodel.DEdim = %d\nmodel.OFdim = %d\nmodel.nparms = %d\n",
			model->DEdim,model->OFdim,model->nparms);
    fprintf(fptr,"model.parmdep =");
    for (j=0; j<model->nparms; j++) fprintf(fptr," %d",model->parmdep[0][j]);
    for (i=1; i<model->DEdim; i++) {
		fprintf(fptr,"\n               ");
		for (j=0; j<model->nparms; j++) fprintf(fptr," %d",model->parmdep[i][j]);
    }
    fprintf(fptr,"\nmodel.vardep =");
    for (j=0; j<model->DEdim; j++) fprintf(fptr," %d",model->vardep[0][j]);
    for (i=1; i<model->DEdim; i++) {
		fprintf(fptr,"\n              ");
		for (j=0; j<model->DEdim; j++) fprintf(fptr," %d",model->vardep[i][j]);
    }
    fprintf(fptr,"\n\nDATA:\n");
    fprintf(fptr,"number of sets (files) = %d\n",data->nsets);
    for (i=0; i<data->nsets; i++) {
		fprintf(fptr,"  data set %d:\n",i);
        fprintf(fptr,"    step size = %f to %f\n",data->set[i].step_size[0],data->set[i].step_size[1]);
		for (j=0; j<model->DEdim; j++) {
			if (data->set[i].y[j].tmin >= 0.0) {
                avgarea = 0;
                for(k = 1; k < data->set[i].y[j].ncpt; k++){
                    avgarea += 0.5*(SPAN(data->set[i].y[j].cpt[k-1])+SPAN(data->set[i].y[j].cpt[k]))*
                                (data->set[i].y[j].t[k]-data->set[i].y[j].t[k-1]);
                }
                avgarea /= data->set[i].y[j].input_t[data->set[i].y[j].npts-1] - data->set[i].y[j].input_t[0];
				fprintf(fptr,"    y[%d] tasle fit with tmin=%f, \n"
						"          yielding height %f with %d control points.\n",
						j,data->set[i].y[j].tmin,avgarea,data->set[i].y[j].ncpt);
			}
			else {
				fprintf(fptr,"    y[%d] range read from file\n",j);
			}
		}
    }
    fprintf(fptr,"\nPARAMETERS:\n");
    fprintf(fptr,"maxboxes = %d\nparms.nboxes = %d\n",parms->maxboxes,parms->nboxes);
    fprintf(fptr,"goal widths:\n");
    for (j=0; j<model->nparms; j++) fprintf(fptr,"  p[%d] : %f\n",j,parms->p_gw[j]);
    totalvol = 0.0;
    nvalid = 0;
    avgdist = 0.0;
    k = 0;
    for (i=0; i<parms->nboxes; i++) {
		if (parms->pbox[i].status > -1) {
			nvalid++;
			vol = 100.0*parms->pbox[i].volume/parms->orig_vol;
			dist = 0.0;
			if ((k = (P_TRUE != NULL))) {
				for (j=0; j<model->nparms; j++) {
					if (parms->pbox[i].p[j][0] > P_TRUE[j] || P_TRUE[j] > parms->pbox[i].p[j][1]) k = 2;
                    dist += pow(0.5*(parms->pbox[i].p[j][0] + parms->pbox[i].p[j][1])-P_TRUE[j],2.0);
				}
				dist = sqrt(dist);
			}
			fprintf(fptr,"    BOX %d, status = %d, vol%% = %10.4e",i,parms->pbox[i].status,vol);
			if (k>0) fprintf(fptr,", dist to true = %10.4e",dist);
			if (k==1) fprintf(fptr,", TRUE PARMS IN THIS BOX:\n");
			else fprintf(fptr,":\n");
			totalvol += vol;
			for (j=0; j<model->nparms; j++) {
				fprintf(fptr,"     p[%d] = %f %f\n",j,parms->pbox[i].p[j][0],parms->pbox[i].p[j][1]);
			}
			avgdist += vol*dist;
		}
    }
    if (totalvol > 0.0) avgdist /= totalvol;
    fprintf(fptr," %d of %d boxes still valid\n",nvalid,parms->nboxes);
    fprintf(fptr," percentage of parameter space still valid = %10.4g\n",totalvol);
    if (k>0) fprintf(fptr," average distance to true parmameter values = %10.4g\n",avgdist);
    if (fclose(fptr) == EOF) {
		fprintf(stderr,"parm_reduction: error closing %s.\n",fname); exit(1);
    }
    if (data->write_vars) {
		for (box=0; box<1; box++) {/* parms->nboxes; box++) { */
			if (parms->pbox[box].status > 0) {
				for (i=0; i<data->nsets; i++) {
					for (j=0; j<model->DEdim; j++) {
						sprintf(datafname,"%s_set%d_var%d_box%d",fname,i,j,box);
						if ((fptr = fopen(datafname,"w")) == NULL) {
							fprintf(stderr,"parm_reduction: error opening %s.\n",datafname); exit(1);
						}
						for (n=0; n<data->set[i].y[j].ncpt; n++) {
							fprintf(fptr,"%16.12g %16.12g %16.12g\n",data->set[i].y[j].t[n],
                                    data->set[i].y[j].cpt[n][0],data->set[i].y[j].cpt[n][1]);
						}
						if (fclose(fptr) == EOF) {
							fprintf(stderr,"parm_reduction: error closing %s.\n",datafname); exit(1);
						}
					}
				}
			}
		}
    }
    for (i=0; i<data->nsets; i++) {
		for (j=0; j<model->DEdim; j++) {
			sprintf(datafname,"%s_set%d_var%d_orig",fname,i,j);
			if ((fptr = fopen(datafname,"w")) == NULL) {
				fprintf(stderr,"parm_reduction: error opening %s.\n",datafname); exit(1);
			}
			for (n=0; n<data->set[i].y[j].ncpt; n++) {
				fprintf(fptr,"%16.12g %16.12g %16.12g\n",data->set[i].y[j].t[n],data->set[i].y[j].cpt[n][0],
                        data->set[i].y[j].cpt[n][1]);
			}
			if (fclose(fptr) == EOF) {
				fprintf(stderr,"parm_reduction: error closing %s.\n",datafname); exit(1);
			}
		}
    }

    /* free model.window */
    free(model->window.split_ratios);
    free(model->window.time_split);
    free(model->window.windep[0]);
    free(model->window.windep[1]);
    free(model->window.windep[2]);
    free(model->window.windep[3]);
    free(model->window.windep[4]);
    free(model->window.windep[5]);
    free(model->window.p_prev_mono);
    for (i=0; i<model->DEdim; i++) free(model->window.windep_msrd[0][i]);
    free(model->window.windep_msrd[0]);
    for (i=0; i<model->DEdim; i++) free(model->window.windep_msrd[1][i]);
    free(model->window.windep_msrd[1]);
    for (i=0; i<model->DEdim; i++) free(model->window.windep_msrd[2][i]);
    free(model->window.windep_msrd[2]);
    free(model->window.p_hull);
    free(model->window.t_ind);
    for (i=0; i<discr_set_list->discr_set[0].ispecs[0]+1; i++) free(model->window.y[i]);
    free(model->window.y);
    for (i=0; i<discr_set_list->discr_set[0].ispecs[0]+1; i++) free(model->window.u[i]);
    free(model->window.u);
    free(model->window.tuse);
    free(model->window.attempt);
    free(model->window.eqn_order);
    free(model->window.partial);
    free(model->window.yext);
    free(model->window.uext);
    for (i=0; i<model->DEdim; i++) free(model->window.yext_sort[i]);
    free(model->window.yext_sort);
    for (i=0; i<model->OFdim; i++) free(model->window.uext_sort[i]);
    free(model->window.uext_sort);
    for(i=0; i<model->DEdim; i++) free(model->window.y_alt[i]);
    for(i=0; i<model->OFdim; i++) free(model->window.u_alt[i]);
    free(model->window.y_alt); 
    free(model->window.y_alt1);
    free(model->window.p_alt); 
    free(model->window.u_alt);
    free(model->window.z);
    free(model->window.p_ord0);
    free(model->window.p_ord1);
    free(model->window.p_ord2);
    free(model->window.p_prev[0]);
    free(model->window.p_prev[1]);
    free(model->window.p_prev[2]);
    
    
    /* free discretizations */
    for (i=0; i<discr_set_list->n_discr_types-1; i++) {
        for (j=0; j<discr_set_list->discr_set[i].n_discr; j++) {
            free(discr_set_list->discr_set[i].discr[j].stencil);
            free(discr_set_list->discr_set[i].discr[j].alpha);
            free(discr_set_list->discr_set[i].discr[j].beta);
        }
        free(discr_set_list->discr_set[i].discr);
        free(discr_set_list->discr_set[i].fspecs);
        free(discr_set_list->discr_set[i].ispecs);
    }
    i = discr_set_list->n_discr_types-1;
    for (j=0; j<discr_set_list->discr_set[i].n_discr; j++) {
        free(discr_set_list->discr_set[i].discr[j].stencil);
        free(discr_set_list->discr_set[i].discr[j].alpha);
        free(discr_set_list->discr_set[i].discr[j].beta);
    }
    free(discr_set_list->discr_set[i].discr);
    free(discr_set_list->discr_set);
    
    /* free pboxes */
    for (box=0; box<parms->nboxes; box++) {
        free(parms->pbox[box].p);
        for (j=0; j<data->nsets; j++) {
            for (i=0; i<model->DEdim; i++) {
                if (data->set[0].unmsrd.y_unmsrd[i]) {
                    free(parms->pbox[box].vp[j][i].cpt);
                    free(parms->pbox[box].vp[j][i].t);
                }
            }
            free(parms->pbox[box].vp[j]);
        }
        free(parms->pbox[box].vp);
        for (j=0; j<model->TEdim; j++) {
            free(parms->pbox[box].inf[j]);
        }
        free(parms->pbox[box].inf);
    }
    free(parms->pbox);
    free(parms->p_gw);

    /* free data */
    for (i=0; i<data->nsets; i++) {
        for (j=0; j<model->DEdim; j++) {
            free(data->set[i].y[j].t);
            free(data->set[i].y[j].cpt);
        }
        for (j=0; j<model->OFdim; j++) {
            free(data->set[i].u[j].t);
            free(data->set[i].u[j].cpt);
        }
        free(data->set[i].y);
        free(data->set[i].u);
        free(data->set[i].change);
        free(data->set[i].unmsrd.y_unmsrd);
        free(data->set[i].unmsrd.y_unmsrd_index);
        free(data->set[i].unmsrd.y_unmsrd_order);
    }
    free(data->y_unmsrd_max);
    free(data->set);
    
    /* free model */
    free(model->p_extreme);
    free(model->y_extreme);
    free(model->ispecs);
    free(model->fspecs);
    free(model->parmdep);
    free(model->vardep);
    free(model->ofdep);
    free(model->parmdep_flag);
    
    /* free other */
    free(P_TRUE);
    return;
}

void reduce_all(struct s_model *model, struct s_discr_set_list *discr_set_list, 
                struct s_data *data, struct s_parms *parms) {
	
    int box,dep_factor,hi,i,j,lo,*sortb,split_direction;
    double inf_max,furthest,slice,temp;
    struct s_pbox *cbox;    
    box = 0; /* number of finished boxes.  A box is finished if its status is good or invalid 
              * or if it is putative_split but the maximum number of boxes has been reached. */
//#if defined(output1) || defined(output4)
    double mu,mu_hull,mu_transformed,p_hull[13][2];
    double p_com[13]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},vol0,vol,volbox,volS,volR; 
    int hull_flag;
    //Range p_temp[7];
//#endif
    /* The sortedboxes array has first the finished boxes, and then
     * the unfinished boxes from largest volume to smallest. */
    if ((sortb = (int*) malloc(parms->maxboxes*sizeof(int))) == NULL) {
		fprintf(stderr,"parm_reduction: error allocating memory for sortb\n"); exit(1);
    }

    sortb[0] = 0;
    while (box < parms->nboxes) {

        cbox = &parms->pbox[sortb[box]]; /* the current box to be used */

		if (cbox->status == putative) {
#ifdef db1
			printf("Calling reduce_box for box=%d...\n",sortb[box]); fflush(NULL);
#endif
			reduce_box(model,discr_set_list,data,cbox,data->minred,parms->p_gw,parms->minred);
            /* If this box is now finished shift the count by one. 
             * Otherwise, place this box in the appropriate place in the sorted array. */
			if (cbox->status == invalid || cbox->status == good) box++;
			else {
                lo = box;
				hi = parms->nboxes;
				while (hi-lo > 1) {
					i = (hi+lo)/2;
					if (parms->pbox[sortb[i]].volume > cbox->volume)
						lo = i;
					else
						hi = i;
				}
				hi = sortb[box];
				for (i=box; i<lo; i++) sortb[i] = sortb[i+1];
				sortb[lo] = hi;
			}
		}
		else if (cbox->status == putative_split) {             
            /* If the current number of boxes is less than the maximum, then split the box */
            if (parms->nboxes < parms->maxboxes) {                
#ifdef db1
                printf("Calling split_box for box=%d...\n",sortb[box]);
#endif                
                /* The parameter range that is the largest relative distance from its goal is identified. */
                lo = 0;
                furthest = -1.0;
                for (i=0; i<model->nparms; i++) {
                    temp = SPAN(cbox->p[i])/parms->p_gw[i];
                    if (temp > furthest) {
                        furthest = temp;
                        lo = i;
                    }
                }
#ifdef db_boxcut
                printf("\tParameter %d is furthest relative to goal.\n",
                       (data->num_y_unmsrd_max > 0 && parms->nboxes < 0.5*parms->maxboxes) ? model->nparms : lo);
                printf("\tChecking all parameters for influence in equations dependent on this parameter.\n");
#endif
                /* set default values for the box split */
                inf_max = 0.0;
                split_direction = lo;
                slice = 0.5;

                /* Now that the parameter furthest from its goal has been identified, search
                 * through each equation which is dependent on this parameter, and determine
                 * which parameter has the largest influence on the value of F_min and F_max.
                 * This is the parameter that should be split since it is the largest impediment
                 * to further progress. */      
#ifdef influence_split
                for (j=0; j<model->nparms; j++) {
#ifdef db_boxcut
                    printf("\t\tParameter %d:\n",j);
#endif
                    for (i=0; i<model->TEdim; i++) {
                        /* If we have used less than half the allowed boxes and there is an 
                         * unmeasured variable, then search among parameters in equations containing
                         * those unmeasured variables. Otherwise, check among parameters in equatons
                         * containing the parameter requiring the most work.*/
                        if (data->num_y_unmsrd_max > 0 && 2*parms->nboxes < parms->maxboxes) {
                            dep_factor = 0;
                            for (hi=0; hi<model->DEdim; hi++) {
                                if (data->y_unmsrd_max[hi] && (model->vardep[i][hi] || i == hi)) {
                                    dep_factor = 1; break;
                                }
                            }
                        }
                        else 
                            dep_factor = model->parmdep[i][lo];
#ifdef db_boxcut
                        if(dep_factor && model->parmdep[i][j]) 
                            printf("\t\t\tin eqn %d, influence is\t%f,%f\n",i,cbox->inf[i][j].pct[0]/
                                   cbox->inf[i][j].ctr[0]++,cbox->inf[i][j].pct[1]/cbox->inf[i][j].ctr[1]++);
#endif
                        temp = 0.5*(cbox->inf[i][j].pct[0]/cbox->inf[i][j].ctr[0] + 
                                    cbox->inf[i][j].pct[1]/cbox->inf[i][j].ctr[1]);
                        if(isnan(temp)) temp = 0.0;
                        if (2*parms->nboxes < parms->maxboxes && dep_factor && 
                            model->parmdep[i][j] && temp > inf_max) {
                            inf_max = temp;
                            split_direction = j;
                            slice = 0.5;
                        }
//                        for (hi=0; hi<2; hi++) {
//                            temp = cbox->inf[i][j].pct[hi]/cbox->inf[i][j].ctr[hi];
//                            if (isnan(temp)) temp = 0.0;
//                            if (2*parms->nboxes < parms->maxboxes && dep_factor && 
//                                model->parmdep[i][j] && temp > inf_max) {
//                                inf_max = temp;
//                                split_direction = j;
//                                temp = (cbox->inf[i][j].pct[hi]/cbox->inf[i][j].ctr[hi]) - 
//                                        (cbox->inf[i][j].pct[!hi]/cbox->inf[i][j].ctr[!hi]);
//                                /* If the difference between the influence on the low edge of a 
//                                 * parameter range and the high edge is large, then make a split
//                                 * closer to 95% of the way across. If the influences are close
//                                 * then make the split in the middle of the range */
//                                temp = 0;
//                                if (isnan(temp)) temp = 0.0;
//                                slice = 0.5 + 0.4*(1-2*hi)*(temp);
//                                slice = 1 - slice;
//                            }
//                        }
                    }
                }
#endif
#ifdef db_boxcut
                printf("\tSplit variable %d. Perform split at fraction %f across range.\n",split_direction,slice);
#endif
                split_box(model,parms,sortb[box],data,split_direction,slice);

                /* If there was a successful split, then enter the two new boxes into
                 * the sorted array.  The current box will have been halved in volume
                 * and a new box of the same size will be in the nboxes-1 location. */
                lo = box;
                hi = parms->nboxes - 1;  /* because the last one is the new box */
                while (hi-lo > 1) {
                    i = (hi+lo)/2;
                    if (parms->pbox[sortb[i]].volume > cbox->volume)
                        lo = i;
                    else
                        hi = i;
                }
                hi = sortb[box];
                for (i=box; i<lo; i++) sortb[i] = sortb[i+1]; 
                sortb[lo] = hi;
                sortb[parms->nboxes-1] = parms->nboxes-1;
                
                lo = box;
                hi = parms->nboxes-1; 
                while (hi-lo > 1) {
                    i = (hi+lo)/2; 
                    if (parms->pbox[sortb[i]].volume > parms->pbox[parms->nboxes-1].volume)
                        lo = i;
                    else
                        hi = i;
                }
                hi = lo + (parms->pbox[parms->nboxes-1].volume <= parms->pbox[sortb[lo]].volume);
                for (i=parms->nboxes-1; i>hi; i--) sortb[i] = sortb[i-1]; 
                sortb[hi] = parms->nboxes-1;
            }
            else box++;
		}
    }    
    
#ifdef unknown_output
    FILE *fptr;
    if ((fptr = fopen("CYTOKINE/tasle_output.txt","w")) == NULL) {
		exit(1);
    }
    Range temp1={0.0,0.0};
    double ave_height = 0.0;
    for (hi = 0; hi < model->DEdim; hi++) {
//    for (hi = 0; hi < 1; hi++) {

        if (data->y_unmsrd_max[hi]) {
//            printf("\n\n");
//            printf("UNKNOWN VARIABLE y[%d] (set 0 only)\n",hi);
//            for (i=0; i<1; i++) {
            for (i=0; i<parms->pbox[0].vp[0][hi].ncpt; i++) {
                lo = 0;
                for (j=0; j<parms->nboxes; j++) {
                    if (parms->pbox[j].status != -1) {
                        if (lo == 0){
                            temp1[0] = parms->pbox[j].vp[0][hi].cpt[i][0];
                            temp1[1] = parms->pbox[j].vp[0][hi].cpt[i][1];
                            lo = 1;
                        }
                        else{
                            temp1[0] = MIN(temp1[0],parms->pbox[j].vp[0][hi].cpt[i][0]);
                            temp1[1] = MAX(temp1[1],parms->pbox[j].vp[0][hi].cpt[i][1]);
//                            temp1[0] += parms->pbox[j].vp[0][hi].cpt[i][0];
//                            temp1[1] += parms->pbox[j].vp[0][hi].cpt[i][1];
                            lo++;
                        }
                    }
                }
                fprintf(fptr,"%f\t%f\t%f\t%f\n",parms->pbox[0].vp[0][hi].t[i],temp1[0],parms->pbox[0].vp[0][hi].t[i],temp1[1]);
//                fprintf(fptr,"%f\t%f\t%f\t%f\n",parms->pbox[0].vp[0][hi].t[i],(temp1[0]+temp1[1])/(2*lo),parms->pbox[0].vp[0][hi].t[i],(temp1[0]+temp1[1])/(2*lo));

                ave_height += SPAN(temp1);
            }
        }
    }
    if (fclose(fptr) == EOF) exit(1);
//    printf("\n\n%f\t",ave_height/parms->pbox[0].vp[0][1].ncpt);
#endif

#ifdef box_summary
    for(i=0; i<parms->nboxes; i++){
//        printf("i=%d, status=%d, ",i,parms->pbox[i].status);
        if(parms->pbox[i].status != -1){
//            printf("p=");
            for(j=0; j<model->nparms; j++){
                printf("%f\t%f\t",parms->pbox[i].p[j][0],parms->pbox[i].p[j][1]);
//                printf("[%f,%f],",parms->pbox[i].p[j][0],parms->pbox[i].p[j][1]);
            }
//            printf("[%f,%f],",parms->pbox[i].vp[0][0].cpt[0][0],parms->pbox[i].vp[0][0].cpt[0][1]);
            printf("\n");
        }
//        printf("\n");
    }
#endif
#ifdef output4
    vol = 0;
    /* This is the calculation for mu in the transformed coordinates */
//    model->nparms = 2; /* SIR only */
//    vol0 = (200-0.01)*(1e3-1e-6); /* SIR */
//    vol0 = 4*499.9*99.9*149.9*49.9*49.9*9.9; /* PHARM */
//    vol0 = 4*499.99*499.99*149.99*49.99*49.99*9.9*90; /* PHARM2 */
//    vol0 = 49.995*59.994*4.9995; /* NLP */
    vol0 = (0.000047-0.0000098)*(0.38-0.047)*(1.99-0.16)*(235-29)*(13-1.4); /* HIV */
    
    for (i=0; i<parms->nboxes; i++) {
        if (parms->pbox[i].status != -1) {
            volbox = 1;
            for (j=0; j<model->nparms; j++) {
                volbox *= SPAN(parms->pbox[i].p[j]);
            }
            vol += volbox;
        }
    }

    mu_transformed = pow(vol/vol0,1/(double)model->nparms);

    /* This is the calculation for mu in the original coordinates */
    vol = 0;

    // make a transformation for each box
    for (i=0; i<parms->nboxes; i++) {
        if (parms->pbox[i].status != -1) {
            // make a transformation of each box (use p_temp as a temporary holder)
//            for (j=0; j<model->nparms; j++) {
//                p_temp[j][0] = parms->pbox[i].p[j][0];
//                p_temp[j][1] = parms->pbox[i].p[j][1];
//            }            
//            parms->pbox[i].p[0][0] = MAX(0.005,1/p_temp[0][1]);
//            parms->pbox[i].p[0][1] = MIN(50.0,1/p_temp[0][0]);
//            
//            parms->pbox[i].p[1][0] = MAX(0.006,p_temp[0][0]/p_temp[2][1]);
//            parms->pbox[i].p[1][1] = MIN(60.0,p_temp[0][1]/p_temp[2][0]);
//            
//            parms->pbox[i].p[2][0] = MAX(0.0005,p_temp[0][0]*p_temp[1][0]/p_temp[2][1]);
//            parms->pbox[i].p[2][1] = MIN(5.0,p_temp[0][1]*p_temp[1][1]/p_temp[2][0]);
//            
//            if (parms->pbox[i].p[0][0] > parms->pbox[i].p[0][1] ||
//                parms->pbox[i].p[1][0] > parms->pbox[i].p[1][1] ||
//                parms->pbox[i].p[2][0] > parms->pbox[i].p[2][1] ) {
//                parms->pbox[i].status = -1;
//            }
//            else {
                volbox = 1;
                for (j=0; j<model->nparms; j++) {
                    volbox *= SPAN(parms->pbox[i].p[j]);
                }
                vol += volbox;
//            }
        }
    }
    mu = MIN(1.0,pow(vol/vol0,1/(double)model->nparms));
    

#endif
#ifdef output4
//    printf("%d\t",parms->maxboxes);
#endif
#ifdef output4
//    printf("%d\t%f\t",data->set[0].unmsrd.base_winsize,data->set[0].step_size[0]);
//    printf("%d\t%f\t%10.4e\t",data->set[0].unmsrd.base_winsize,data->set[0].step_size[0],mu);
    printf("%d\t%f\t%d\t%d\t",data->set[0].unmsrd.base_winsize,data->set[0].step_size[0],model->window.subbox_specs[0],
           model->window.subbox_specs[1]);

#endif
//#ifdef output4
    vol0 = 0; volS = 0; volR = 0;
    for (i=0; i<parms->nboxes; i++) {
        if (parms->pbox[i].status != -1) {
            volbox = 1;
            for (j=0; j<model->nparms; j++) {
                volbox *= SPAN(parms->pbox[i].p[j]);
            }
            vol0 += volbox;
//            volS += SPAN(parms->pbox[i].vp[0][0].cpt[0]);
//            volR += SPAN(parms->pbox[i].vp[0][2].cpt[0]);
        }
    }
    for (i=0; i<parms->nboxes; i++) {
        if (parms->pbox[i].status != -1) {
            volbox = 1;
            for (j=0; j<model->nparms; j++) {
                volbox *= SPAN(parms->pbox[i].p[j]);
            }
            for (j=0; j<model->nparms; j++) {
                p_com[j] += (0.5*(parms->pbox[i].p[j][0]+parms->pbox[i].p[j][1])*volbox)/vol0;
            }
//            p_com[3] += (0.5*(parms->pbox[i].vp[0][0].cpt[0][0]+parms->pbox[i].vp[0][0].cpt[0][1])*
//                         SPAN(parms->pbox[i].vp[0][0].cpt[0]))/volS;
//            p_com[4] += (0.5*(parms->pbox[i].vp[0][2].cpt[0][0]+parms->pbox[i].vp[0][2].cpt[0][1])*
//                         SPAN(parms->pbox[i].vp[0][2].cpt[0]))/volR;
        }
    }
//    printf("%f\t%f\t",p_com[0],p_com[1]);
//    d = 0;
//    for (j=0; j<3; j++) {
//        d += pow(1-p_com[j]/p_true[j],2);
//    }
//    printf("%10.4e\t",sqrt(d));
    hull_flag = 1; /* hull */
    for (i=0; i<parms->nboxes; i++) { 
        if (parms->pbox[i].status != -1) {
            if (hull_flag) {
                for (j=0; j<model->nparms; j++) {
                    p_hull[j][0] = parms->pbox[i].p[j][0];
                    p_hull[j][1] = parms->pbox[i].p[j][1];
                }
//                p_hull[3][0] = parms->pbox[i].vp[0][0].cpt[0][0];
//                p_hull[3][1] = parms->pbox[i].vp[0][0].cpt[0][1];
//                p_hull[4][0] = parms->pbox[i].vp[0][2].cpt[0][0];
//                p_hull[4][1] = parms->pbox[i].vp[0][2].cpt[0][1];
                hull_flag = 0;
            }
            else {
                for (j=0; j<model->nparms; j++) {
                    p_hull[j][0] = MIN(p_hull[j][0],parms->pbox[i].p[j][0]);
                    p_hull[j][1] = MAX(p_hull[j][1],parms->pbox[i].p[j][1]);
                }

//                p_hull[3][0] = MIN(p_hull[3][0],parms->pbox[i].vp[0][0].cpt[0][0]);
//                p_hull[3][1] = MAX(p_hull[3][1],parms->pbox[i].vp[0][0].cpt[0][1]);
//                p_hull[4][0] = MIN(p_hull[4][0],parms->pbox[i].vp[0][2].cpt[0][0]);
//                p_hull[4][1] = MAX(p_hull[4][1],parms->pbox[i].vp[0][2].cpt[0][1]);
            }
        }
    }
//    printf("%f\t%f\t%f\t",p_hull[3][0],p_hull[3][1],p_com[3]);
//    printf("%f\t%f\t%f\t",p_hull[4][0],p_hull[4][1],p_com[4]);
//    printf("%f\t%f\t%f\t%f\t%f\t%f\t",p_hull[0][0],p_hull[0][1],p_com[0],p_hull[1][0],p_hull[1][1],p_com[1]);

//    printf("%f\t%f\t%f\t%f\t%f\t%f\t",p_hull[3][0],p_hull[3][1],p_com[3],0.0,0.0,0.0);
//    printf("%f\t%f\t%f\t%f\t%f\t%f\t",0.0,0.0,0.0,p_hull[4][0],p_hull[4][1],p_com[4]);
//    printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t",p_hull[0][0],p_hull[0][1],p_com[0],p_hull[1][0],p_hull[1][1],p_com[1],p_hull[2][0],p_hull[2][1],p_com[2]);
    for (i=0; i<model->nparms; i++) printf("%f\t%f\t%f\n",p_hull[i][0],p_hull[i][1],p_com[i]);
//#endif
#ifdef output4
//    vol0 = 49.995*59.994*4.9995; /* NLP */
    vol0 = (0.000047-0.0000098)*(0.38-0.047)*(1.99-0.16)*(235-29)*(13-1.4); /* HIV */
//    vol0 = 4*499.99*499.99*149.99*49.99*49.99*9.9*90; /* PHARM2 */
//    vol0 = (200-0.01)*(1e3-1e-6); /* SIR */
    vol = 1;
    for (j=0; j<model->nparms; j++) {
        vol *= SPAN(p_hull[j]);
    }
    mu_hull = pow(vol/vol0,1/(double)model->nparms);
//    printf("%10.4e\t",mu_hull);
//    printf("%10.4e\t%10.4e\t%10.4e\t",mu_transformed,mu,mu_hull);
    printf("%10.4e\t%10.4e\t",mu,mu_hull);
#endif
#ifdef E_test
    discr_error_test(model,discr_set_list,data,parms);
#endif 
    free(sortb);
    
    return;
}

void split_box(struct s_model *model, struct s_parms *parms, int boxnum, 
			  struct s_data *data, int split_direction, double slice) {
	
    /* This function allocates memory in a new box, and splits along the appropriate
     * parameter by the appropriate slice amount. */
    double temp;
    int i,j;
    struct s_pbox *newbox,*oldbox;

    oldbox = &parms->pbox[boxnum];
    newbox = &parms->pbox[parms->nboxes];
    oldbox->status = putative;
    newbox->status = putative;
    if (((newbox->p = (Range *) malloc((model->nparms+data->num_y_unmsrd_max*
                                        data->set[0].unmsrd.base_winsize)*sizeof(Range))) == NULL) ||
        ((newbox->inf = (struct s_inf **) malloc(model->TEdim*sizeof(struct s_inf*))) == NULL)) {
        fprintf(stderr,"parm_reduction: error allocating memory for newbox->p,inf\n"); exit(1);
    }
    
    if ((newbox->vp = (struct s_yu **) malloc(data->nsets*sizeof(struct s_yu *))) == NULL) {
		fprintf(stderr,"parm_reduction: error allocating memory for newbox->vp\n"); exit(1);
    }    
    for (j=0; j<data->nsets; j++) {
        if (((newbox->vp[j] = (struct s_yu *) malloc(model->DEdim*sizeof(struct s_yu))) == NULL)) {
            fprintf(stderr,"parm_reduction: error allocating memory for newbox->vp,vp.cpt,vp.t\n");
            exit(1);
        }
        for (i=0; i<model->DEdim; i++) {
            if (data->set[j].unmsrd.y_unmsrd[i]) {
                if (((newbox->vp[j][i].cpt = (Range *) malloc(oldbox->vp[j][i].ncpt*sizeof(Range))) == NULL) ||
                    ((newbox->vp[j][i].t = (double *) malloc(oldbox->vp[j][i].ncpt*sizeof(double))) == NULL)) {
                    fprintf(stderr,"parm_reduction: error allocating memory for newbox->vp.cpt,t\n"); exit(1);
                }                
            }
        }
    }

    for (i=0; i<model->TEdim; i++) {
        if (((newbox->inf[i] = (struct s_inf*) calloc((model->nparms+1),sizeof(struct s_inf))) == NULL)) {
            fprintf(stderr,"parm_reduction: error allocating memory for inf[i]\n"); exit(1);
        }
    }
    memcpy(newbox->p,oldbox->p,model->nparms*sizeof(Range));    
    temp = (1-slice)*oldbox->p[split_direction][0] + slice*oldbox->p[split_direction][1];
    oldbox->p[split_direction][1] = temp;
    newbox->p[split_direction][0] = temp;
    newbox->volume = slice*oldbox->volume;
    oldbox->volume = (1-slice)*oldbox->volume;
    parms->nboxes++;
    for (j=0; j<data->nsets; j++) {
        for (i=0; i<model->DEdim; i++) {
            if (data->set[j].unmsrd.y_unmsrd[i]) {
                memcpy(newbox->vp[j][i].cpt,oldbox->vp[j][i].cpt,oldbox->vp[j][i].ncpt*sizeof(Range));
                memcpy(newbox->vp[j][i].t,oldbox->vp[j][i].t,oldbox->vp[j][i].ncpt*sizeof(double));
                newbox->vp[j][i].ncpt = oldbox->vp[j][i].ncpt;                
            }
        }
    }
    return;
}

void reduce_box(struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_data *data,
                struct s_pbox *pbox, double y_minred, double *p_gw, double minred) {
	
    /* On input pbox->status is "putative".  On return of this function, 
     * pbox->status will be one of "invalid", "putative_split", or "good". 
     * Continue to reduce variable and parameter ranges as long as substantial
     * progress is being made and the parameter goals have not been met. */
    int i,set_progress,num_noprogress;
	
    if (parm_goals_met(pbox,model->nparms,p_gw)) pbox->status = good;
    num_noprogress = 0;  /* number of consecutive sets with no significant progress */
    i = -1;
    while (pbox->status == putative) {
		if (++i == data->nsets) i=0;
#ifdef db1
        int j;
        for(j=0; j<model->nparms; j++) printf("    p_init[%d] = [%f,%f]\n",j,pbox->p[j][0],pbox->p[j][1]);
#endif
#ifdef db2
		printf("  Calling process_set...\n");
#endif
		/* now we are going to make some progress */
        pbox->vp_set = i;
		set_progress = process_set(model,discr_set_list,&data->set[i],pbox,y_minred,p_gw,minred);

		if (set_progress < 0) { /* inconsistency found */
			pbox->status = invalid;
		}
		else if ((set_progress >= 1) && (parm_goals_met(pbox,model->nparms,p_gw))) {
			/* some progress made and goals met */
			pbox->status = good;
		}
		else if ((set_progress == 2)) {
			/* significant progress for this set */
			num_noprogress = 0;
		}
		else if ((set_progress <= 1) && (++num_noprogress == data->nsets)) { 
			/* no significant progress for all sets */
			pbox->status = putative_split;
		}
    }
    /* recalculate the box volume */
    pbox->volume = 1.0;
    for (i=0; i<model->nparms; i++) pbox->volume *= SPAN(pbox->p[i]);    
    
    return;
}

int parm_goals_met(struct s_pbox *pbox, int nparms, double *p_gw) {
    
    /* Determine if a given parameter has reached its goal width */
	int i;
    for (i=0; i<nparms; i++) { 
        if (SPAN(pbox->p[i]) > p_gw[i]) return 0; 
    }
    return 1;
}

int process_set(struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_set *set,
                struct s_pbox *pbox, double y_minred, double *p_gw, double minred) { 
	
    /* return value:  -1 if there is an inconsistency, 
     *                 0 if no update was made, 
     *                 1 if an insignificant update was made.
     *                 2 if a significant update was made. */
    int i,progress=0,settling=1,win_progress,winsize;
    double durmax,h,max_durmax=0.0,noprog_duration=0.0,overall_duration=0.0;
    double returned_value,settling_duration,settling_durmax,winstart;
    Range *p_prev=model->window.p_prev[0];
	
    /* Keep a copy of the parameter box on input. */
    memcpy(p_prev,pbox->p,model->tparms*sizeof(Range));

    /* Set the overall user inputted window size */
    winsize = discr_set_list->discr_set[0].discr[discr_set_list->discr_set[0].n_discr-1].winsize;
        
#ifdef db2
	printf("      Calling select_window...\n");
#endif
    /* Select an initial window for processing. */
    settling_duration = select_window(set,0.0,3,winsize,&winstart,&h,discr_set_list->loop_specs[0]);
    
    /* Cycle through windows until the maximum has been reached, or an inconsistency is found.
     * First, we cycle through windows for a minimum width of 'settling_durmax' overlapping
     * windows by a user-specified amount. Then we continue until we have made no significant
     * progress for a minimum of 'durmax' amount of window width, up to a max of 'max_durmax' */
    
    durmax = (discr_set_list->loop_specs[2])*SPAN(set->t_range);
    settling_durmax = (discr_set_list->loop_specs[1])*SPAN(set->t_range)/(1-discr_set_list->loop_specs[0]);
    max_durmax = MAX(durmax,(discr_set_list->loop_specs[3])*SPAN(set->t_range));

    while (settling || (noprog_duration < durmax && overall_duration < max_durmax)) {

        if (settling_duration >= settling_durmax) settling = 0;
#ifdef db2
		printf("      selected winstart = %lf,  with h = %lf\n",winstart,h);
        printf("      Calling process_window...\n");
#endif
        /* process chosen window */
        win_progress = process_window(model,discr_set_list,set,pbox,0,y_minred,p_gw,minred,winstart,h);

        if (win_progress < 0) return -1;
        else if (win_progress > 1) noprog_duration = 0.0;
        
        /* choose another window */
        returned_value = select_window(set,winstart+(winsize-1)*h,(settling == 1) ? -2 : win_progress,
                                       winsize,&winstart,&h,discr_set_list->loop_specs[0]);

        if (settling) settling_duration += returned_value;  
        else {
            noprog_duration += returned_value;
            overall_duration += returned_value;
        }
    }
    
    /* If the progress made on a regular parameter is significant, then return 2. If some progress is made
       on a parameter, then return 1. Else, return 0 (default). This information is used in reduce_box to decide
       whether or not to split, continue reducing, or denote as finished. */
    for (i=0; i<model->nparms; i++) {
        if ((SPAN(p_prev[i]) - SPAN(pbox->p[i])) > minred*p_gw[i]) {
            progress = 2; break;
        }
        else if ((SPAN(p_prev[i]) - SPAN(pbox->p[i])) > 0.0) {
            progress = 1;
        }
    }
#ifdef db2
	printf("    process_set returning progress = %d\n",progress);
#endif
    return progress;
}

int process_window(struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_set *set,
                   struct s_pbox *pbox, int eqn, double y_minred, double *p_gw,
                   double minred, double winstart, double h) {
    /* return value:  -1 if there is an inconsistency, 
     *                 0 if no update was made, 
     *                 1 if an insignificant update was made. 
     *                 2 if a significant update was made.  */
        
    struct s_unmsrd *unmsrd=&set->unmsrd;
    int base_winsize = unmsrd->base_winsize,eqn_progress=0,guess=0,i,j,k,m,mono,mono2;
    int progress=0,split_option,unmsrd_progress,*t_ind=model->window.t_ind;
    Range *p_prev=model->window.p_prev[1],**u=model->window.u,diff,**y=model->window.y;
    double prop,*tuse=model->window.tuse;
    intRange *attempt=model->window.attempt;
    
    /* get information required for processing the first time */
    for (i=0; i<base_winsize; i++) tuse[i] = winstart + i*h;
            
    /* this array allows us to stop making attempting to make progress if a particular sub-box
     fails to make progress for a particular parameter */
    for (i=0; i<model->tparms; i++) for(j=0; j<2; j++) attempt[i][j] = 1;

    /* Get the ranges for the y values at the time points within the discretization. */
    for (j=0; j<model->DEdim; j++) {
        for (i=0; i<base_winsize; i++) {
#ifdef db3
            printf("         Calling get_range from process_equation DEdim loop...\n");
#endif
            if(set->unmsrd.y_unmsrd[j]){ /* unmeasured */
                m = unmsrd->y_unmsrd_index[j]; /* location of unmeasured variable in pbox array */
                /* index shift */
                t_ind[j] = get_range(tuse[i],&pbox->vp[pbox->vp_set][j],t_ind[j],y[i][j]);
                pbox->p[m+i][0] = y[i][j][0];
                pbox->p[m+i][1] = y[i][j][1];
            }
            else{
                t_ind[j] = get_range(tuse[i],&set->y[j],t_ind[j],y[i][j]);
            }
            /* save data in another format for use later */
            model->window.y_alt[j][i][0] = y[i][j][0];
            model->window.y_alt[j][i][1] = y[i][j][1];
        }

        /* Determine the maximum/minimum values for all measured variables,
         * for use by the dependency checks later in the program. This code sorts the 
         * time steps for each variable from highest to lowest so that when a subwindow
         * is used, it is easy to find the largest and smallest value within that window. */
        if (!set->unmsrd.y_unmsrd[j]) {
            model->window.yext_sort[j][0][0] = 0;
            model->window.yext_sort[j][0][1] = 0;
            for (i=1; i<base_winsize; i++) {
                /* place the lower values */
                k = i;
                while (k > 0 && y[i][j][0] < y[model->window.yext_sort[j][k-1][0]][j][0]) {
                    model->window.yext_sort[j][k][0] = model->window.yext_sort[j][k-1][0]; k--;
                }
                model->window.yext_sort[j][k][0] = i;
                
                /* place the upper values */
                k = i;
                while (k > 0 && y[i][j][1] < y[model->window.yext_sort[j][k-1][1]][j][1]) {
                    model->window.yext_sort[j][k][1] = model->window.yext_sort[j][k-1][1]; k--;
                }
                model->window.yext_sort[j][k][1] = i;
            }
        }
    }
    /* Get the ranges for the u values at the time points within the discretization. */
    for (j=0; j<model->OFdim; j++) {
        for (i=0; i<base_winsize; i++) {
#ifdef db3
            printf("         Calling get_range from process_equation OFdim loop...\n");
#endif
            t_ind[j] = get_range(tuse[i],&set->u[j],t_ind[j],u[i][j]);
            model->window.u_alt[j][i][0] = u[i][j][0];
            model->window.u_alt[j][i][1] = u[i][j][1];
        }
        
        model->window.uext_sort[j][0][0] = 0;
        model->window.uext_sort[j][0][1] = 0;
        for (i=1; i<base_winsize; i++) {
            /* place the lower values */
            k = i;
            while (k > 0 && u[i][j][0] < u[model->window.uext_sort[j][k-1][0]][j][0]) {
                model->window.uext_sort[j][k][0] = model->window.uext_sort[j][k-1][0]; k--;
            }
            model->window.uext_sort[j][k][0] = i;
            
            /* place the upper values */
            k = i;
            while (k > 0 && u[i][j][1] < u[model->window.uext_sort[j][k-1][1]][j][1]) {
                model->window.uext_sort[j][k][1] = model->window.uext_sort[j][k-1][1]; k--;
            }
            model->window.uext_sort[j][k][1] = i;
        }
    }
    
    /* First allocate some space, and keep a copy of the parameter box on input. */
    memcpy(p_prev,pbox->p,model->tparms*sizeof(Range));
    
    /* Loop through all equations with the full box */
    for (j = 0; j < model->TEdim; j++) {
        
        eqn = model->window.eqn_order[j]; /* use this equation */
        
        /* Process the equation if there is something to reduce */
        if (model->parmdep_flag[eqn]) {
            /* Process as either a ODE or an algebraic equation. The monotonicity check
             * for the algebraic equation happens within the function called. */
            if (eqn < model->DEdim) {                     
#ifdef db2
                printf("         Calling process_equation on eqn %d...\n",eqn);
#endif         
                eqn_progress = process_equation(model,discr_set_list,unmsrd,pbox,eqn,y_minred,tuse,y,u,
                                                h,mono,mono2,attempt);
            }
            else if (eqn >= model->DEdim) {
#ifdef db2
                printf("         Calling process_algebraic_equation on eqn %d...\n",eqn);
#endif         
                eqn_progress = process_algebraic_equation(model,discr_set_list,unmsrd,pbox,eqn,y_minred,
                                                          tuse,y,u,mono,mono2,attempt);
            }
            if (eqn_progress < 0) break;
        }
    }
    
    /* Now, if there are unmeasured variables and the user has indicated that subboxes should
     * be used, process these subboxes */
    if ((split_option = model->window.subbox_specs[2]) > 0 && eqn_progress >= 0) {
        /* For each unmeasured variable, create subboxes based on splitting those
         * unmeasured parameters */
//        for (j = 0; j < set->unmsrd.num_y_unmsrd; j++) {
        j = 0;
#ifdef db2
            printf("       Calling process_subboxes for unmeasured var = %d...\n",unmsrd->y_unmsrd_order[j]);
#endif   
            unmsrd_progress = process_subboxes(model,discr_set_list,set,pbox,eqn,y_minred,tuse,y,u,
                                               h,mono,mono2,attempt,split_option,j);
            if (unmsrd_progress < 0) {
                eqn_progress = unmsrd_progress; 
//                break;
            }
            else eqn_progress = MAX(unmsrd_progress,eqn_progress);
//        }
    }
    
    /* If there was at least one box with no inconsistencies, check significance of progress */
    if (eqn_progress >= 0) {
        for (i=0; i<model->tparms; i++){
            if (SPAN(pbox->p[i]) < SPAN(p_prev[i])) {
                progress = 1;
                if (i < model->nparms)
                    prop = minred*p_gw[i];
                else
                    /* find which unmeasured variable we are dealing with */
                    prop = y_minred*set->y[(i-model->nparms)/(base_winsize)].gw;
                
                if (SPAN(p_prev[i]) - SPAN(pbox->p[i]) > prop) {
#ifdef db2
                    printf("        parm[%d] significant change from [%lf,%lf] to [%lf,%lf]\n",
                           i,p_prev[i][0],p_prev[i][1],pbox->p[i][0],pbox->p[i][1]);
#endif
                    progress = 2; 
//                    break; /* Change this for speed */
                }
            }
        }
        
        /* Put the special parameters selected for this window back into the control points. */
        if (progress > 0) {
            for (i = model->nparms; i < model->tparms; i++) {     
                /* Find which variable we are dealing with */
                k = set->unmsrd.y_unmsrd_order[(i-model->nparms)/(base_winsize)];
                
                guess_binsearch(tuse[0]+h*((i-model->nparms)%base_winsize),
                                pbox->vp[pbox->vp_set][k].ncpt,pbox->vp[pbox->vp_set][k].t,&guess);
                if(guess == 0) guess++;
                
                prop = (pbox->vp[pbox->vp_set][k].t[guess]-tuse[0]-h*((i-model->nparms)%base_winsize))/
                (pbox->vp[pbox->vp_set][k].t[guess]-pbox->vp[pbox->vp_set][k].t[guess-1]);
                
                for (j = 1; j >= 0; j--) {
//                    if (p_prev[i][1] < pbox->vp[pbox->vp_set][k].cpt[guess-j][1]) 
//                        diff[1] = p_prev[i][1] - pbox->p[i][1];
                    if (p_prev[i][1] > pbox->p[i][1])
                        diff[1] = MAX(0,pbox->vp[pbox->vp_set][k].cpt[guess-j][1]-pbox->p[i][1]);
                    else
                        diff[1] = 0.0;
                    
//                    if (p_prev[i][0] > pbox->vp[pbox->vp_set][k].cpt[guess-j][0]) 
//                        diff[0] = pbox->p[i][0] - p_prev[i][0];
                    if (pbox->p[i][0] > p_prev[i][0])
                        diff[0] = MAX(0,pbox->p[i][0]-pbox->vp[pbox->vp_set][k].cpt[guess-j][0]);
                    else
                        diff[0] = 0.0;
                    
                    if (j == 0) prop = 1-prop;
                    pbox->vp[pbox->vp_set][k].cpt[guess-j][0] += prop*diff[0];
                    pbox->vp[pbox->vp_set][k].cpt[guess-j][1] -= prop*diff[1];
                    
                    if (pbox->vp[pbox->vp_set][k].cpt[guess-j][0] > pbox->vp[pbox->vp_set][k].cpt[guess-j][1]) {
                        progress = -1;
                        i = model->tparms; break;
                    }
                }
            }
        }
    }
    else progress = -1;
#ifdef db2
	printf("      Process_window returning progress = %d\n",progress);
#endif
//    if (progress == -1) {
//    printf("current:\n");
//    for (i=0; i<pbox->vp[pbox->vp_set][0].ncpt; i++) {
//        printf("%f\t%f\t%f\t%f\n",pbox->vp[pbox->vp_set][0].t[i],pbox->vp[pbox->vp_set][0].cpt[i][0],pbox->vp[pbox->vp_set][0].t[i],pbox->vp[pbox->vp_set][0].cpt[i][1]);
//    }
//    }    
//    
    return progress;
}

int process_subboxes (struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_set *set,
                      struct s_pbox *pbox, int eqn, double y_minred, double *t, Range **y, Range **u, 
                      double h, int mono, int mono2, intRange *attempt, int split_option, int p_index) {
    
    struct s_unmsrd *unmsrd = &set->unmsrd;
    int base_winsize=unmsrd->base_winsize,box_split,current_split,eqn_progress;
    int hull_obtained=0,i,index,j,max_subboxes,subbox,subbox_use,subbox_flag=0,temp,time_split;
    double *split=model->window.split_ratios;
    Range *p_hull=model->window.p_hull,*p_prev=model->window.p_prev[2];
    
    /* Store pbox-p for reseting after each subbox */
    memcpy(p_prev,pbox->p,model->tparms*sizeof(Range));
    
    /* this array allows us to stop making attempting to make progress if a particular sub-box
     fails to make progress for a particular parameter */
    for (i=0; i<model->tparms; i++) for(j=0; j<2; j++) attempt[i][j] = 1;
    
    /* If there are special unknown-variable parameters, then split the first one into subboxes
     split_option = { 0 = no subboxes (not called here), 
                    { 1 = split the m-th time step into n subboxes (m,n given as specs[0],[1])
                    { 2 = split m time steps into n subboxes each (m,n given as specs[0],[1]) */
    time_split = model->window.subbox_specs[0];
    box_split = model->window.subbox_specs[1];
    split = model->window.split_ratios;
    index = model->nparms + p_index*base_winsize; /* which group of parameters to be split */ 
    /* Check to see if enough progress has been made on the unmeasured that splitting isn't needed */
    j = 0;
    for (i=0; i<base_winsize; i++)
        if (SPAN(p_prev[index+i]) < set->y[unmsrd->y_unmsrd_order[p_index]].gw) j++; 
    if (j == base_winsize) split_option = 0;
    max_subboxes = (split_option == 1) ? 
        pow(box_split,unmsrd->num_y_unmsrd) : pow(pow(box_split,time_split),unmsrd->num_y_unmsrd);
    if (split_option == 1) index += time_split;

    /* Loop through all subboxes */
    for (subbox = 0; subbox < max_subboxes; subbox++) {  
        /* Reset the parameters for the next subbox */
        memcpy(pbox->p,p_prev,model->tparms*sizeof(Range));
        
#ifdef db2
        printf("       Processing subbox = %d/%d...\n",subbox,max_subboxes-1);
#endif
//        printf("       Processing subbox = %d/%d...\n",subbox,max_subboxes-1);

        /* Set up subboxes here according to the option chosen above.*/
        if (split_option == 1) {
            memcpy(pbox->p,p_prev,model->tparms*sizeof(Range));
            subbox_use = subbox;
            for (j=0; j<unmsrd->num_y_unmsrd; j++) {
                temp = subbox_use % (int)pow(box_split,j+1);
                index = model->nparms + j*base_winsize;
                pbox->p[index][0] = p_prev[index][0] + split[temp]*SPAN(p_prev[index]);
                pbox->p[index][1] = p_prev[index][0] + split[temp+1]*SPAN(p_prev[index]);  
                subbox_use = (subbox_use - temp)/pow(box_split,j+1);
//                printf("pbox->p[%d] = [%f,%f]\n",index,pbox->p[index][0],pbox->p[index][1]);
            }
        }
        else if (split_option == 2) {
            /* WAIT! the new splitting has not yet been done on option 2 */
            memcpy(pbox->p,p_prev,model->tparms*sizeof(Range));
            j = 0;
            for (i=0; i<base_winsize; i++) { /* only reset first unmeasured */
                if (model->window.time_split[i]) {
                    current_split = (subbox/(int)(pow(box_split,time_split-j-1))) % box_split;
                    j++;
                    pbox->p[index+i][0] = p_prev[index+i][0] + (split[current_split])*(SPAN(p_prev[index+i]));
                    pbox->p[index+i][1] = p_prev[index+i][0] + (split[current_split+1])*(SPAN(p_prev[index+i]));                    
                }
            }
        }

        /* The equations have been previously ordered by usefulness. Equations that have
         * the unmeasured variables are ranked first, followed by any algebraic equations, 
         * with all other equations following in order. */
        for (j = 0; j < model->TEdim; j++) {
                        
            eqn = model->window.eqn_order[j]; /* use this equation */
            
            /* Process the equation in the subbox only if the equation is dependent on the 
             * unmeasured variable, or it is the first box, and there is a parameter to reduce */
            if (model->parmdep_flag[eqn] == 2) {
                /* Process as either a ODE or an algebraic equation. The monotonicity check
                 * for the algebraic equation happens within the function called. */
                if (eqn < model->DEdim) {                     
#ifdef db2
                    printf("         Calling process_equation on eqn %d...\n",eqn);
#endif         
                    eqn_progress = process_equation(model,discr_set_list,unmsrd,pbox,eqn,y_minred,t,y,u,
                                                    h,mono,mono2,attempt);
                }
                else if (eqn >= model->DEdim) {
#ifdef db2
                    printf("         Calling process_algebraic_equation on eqn %d...\n",eqn);
#endif         
                    eqn_progress = process_algebraic_equation(model,discr_set_list,unmsrd,pbox,eqn,y_minred,
                                                              t,y,u,mono,mono2,attempt);
                }
                if (eqn_progress < 0) break;
            }
        }
    
        /* now that we have looped through each equation, we create an overall hull over all subboxes.
         * Thus, if the first subbox resulted in progress, but the second one didn't, we cannot
         * make overall progress. */
        if (eqn_progress >= 0) {
            if (subbox_flag) {
                for (i=0; i<model->tparms; i++) {
                    p_hull[i][0] = MIN(p_hull[i][0],pbox->p[i][0]);
                    p_hull[i][1] = MAX(p_hull[i][1],pbox->p[i][1]);
                }
            }
            else {
                memcpy(p_hull,pbox->p,model->tparms*sizeof(Range));
                subbox_flag = 1;
            }
            /* If no progress has been made, check to see if the hull for all parameters
             * is already maximal. If so, stop progress on subboxes. */
            for (i=0; i<model->tparms; i++) {
                if (attempt[i][0] && (p_hull[i][0] <= p_prev[i][0])) {
                    attempt[i][0] = 0; hull_obtained++;
                }
                if (attempt[i][1] && (p_hull[i][1] >= p_prev[i][1])) {
                    attempt[i][1] = 0; hull_obtained++;
                }
            }
            if (hull_obtained >= 2*model->tparms) break;
        } 
    }
        
    /* If something was placed into the hull at some point, then copy the hull back into pbox,
     * otherwise, there was an inconsistency in every box, so we can discard the box */
    if (subbox_flag) 
        memcpy(pbox->p,p_hull,model->tparms*sizeof(Range));
    else 
        subbox_flag = -1;
    
    return subbox_flag;
}

int get_range (double t_req, struct s_yu *v, int t_ind, Range val) {
    
    /* Given a requested time value, output the appropriate range of the 
     * dependent variable (y or u).  The return value is the index into the time
     * vector for the nearest control point equal to or before t_req. 
     * t_ind is a suggested index into the time array that is close to t_req.  */
    int hi,ind,lo;
    double req;
	
    if (t_ind >= v->ncpt) ind = v->ncpt-1;
    else if (t_ind < 0) ind = 0;
    else ind = t_ind;
    req = t_req*DBL_EPSILON;
    if (t_req - v->t[0] < -req || t_req - v->t[v->ncpt-1] > req) {
		fprintf(stderr,"parm_reduction: error in get_range.  Requested t is %lf but vector"
				" limits are [%lf,%lf].\n",t_req,v->t[0],v->t[v->ncpt-1]); exit(1);
    }
    else { /* This deals with round off error */
		req = MAX(MIN(v->t[v->ncpt-1],t_req),v->t[0]);
    }
#ifdef db3
	printf("           calling guess_binsearch...\n");
#endif
    if ((hi = guess_binsearch(req,v->ncpt,v->t,&ind)) == 0) lo = hi++;
    else lo = hi - 1;
    
    req = (req - v->t[lo])/(v->t[hi] - v->t[lo]);
    val[0] = v->cpt[lo][0] + (v->cpt[hi][0] - v->cpt[lo][0])*req;
    val[1] = v->cpt[lo][1] + (v->cpt[hi][1] - v->cpt[lo][1])*req;

#ifdef db3
	printf("           get_range returning ind = %d\n",ind);
#endif
    return ind;
}

int guess_binsearch (double req, int n, double *v, int *guess) {
    /* Binary search to find the smallest index, ind, into the nondecreasing vector v 
     * (of length n) which gives req <= v[ind].  The following assumptions are made:
     * 1) v is nondecreasing
     * 2) n > 1
     * 3) 0 <= guess <= n-1
     * 4) v[0] <= req <= v[n-1]
     * On input, guess is an initial guess for ind; on output it is ind. */
	
    int hi,i,lo;
    if (req > v[*guess]) {
		if (req <= v[*guess+1]) return ++(*guess);
		lo = *guess;
		hi = n-1;
    }
    else {
		if (*guess == 0 || req > v[*guess-1]) return *guess;
		lo = 0;
		hi = *guess;
    }
    while (hi - lo > 1) {
		if (v[(i = (lo + hi)/2)] < req) lo = i;
		else hi = i;
    }
    return (*guess = hi);
}

void check_monotonicity (struct s_model *model, struct s_discr *discr, struct s_unmsrd *unmsrd,
                         struct s_pbox *pbox, int eqn, double *t, Range **y, Range **u, double h, 
                         int base_winsize, int *windep[3], int *mono, int *mono2) {

    /* Function inserts monotonicity properties into windep, and into the values mono and mono2 */
    Range partial,text,*p_prev = model->window.p_prev_mono;
    int cor,flag,i,j,k;
    double alpha_flag,temp;
    static int map[2][3] = {
		{ 0, 0, 1},
		{ 1, 0, 0}};
    
    /* We first determine extreme values over the given interval if needed for the equation. */
    for (j=0; j<model->DEdim; j++) {
        /* Update any unmeasured variables first */
        if (unmsrd->y_unmsrd[j] && (model->vardep[eqn][j] || eqn==j)) {
                k = unmsrd->y_unmsrd_index[j]; /* location of unmeasured variables in pbox array */
                for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) {
                    y[i][j][0] = pbox->p[k+i][0];
                    y[i][j][1] = pbox->p[k+i][1];
                    model->window.y_alt[j][i][0] = pbox->p[k+i][0];
                    model->window.y_alt[j][i][1] = pbox->p[k+i][1];
                }
        }
        /* If the variable appears in any partial derivative, then determine its max/min values
         * If the variable is measured, then you can use the max/min ordering performed earlier.
         * If the variable is unmeasured, then some more work needs to be done */
        if (model->varpartialdep[eqn][j]) {
            if (unmsrd->y_unmsrd[j]) {
                model->window.yext[j][0] = y[discr->sub_window[0]][j][0];
                model->window.yext[j][1] = y[discr->sub_window[0]][j][1];

                for (i=discr->sub_window[0]+1; i<discr->sub_window[1]; i++) {   
                    model->window.yext[j][0] = MIN(model->window.yext[j][0],y[i][j][0]);
                    model->window.yext[j][1] = MAX(model->window.yext[j][1],y[i][j][1]);
                }
            }
            else {
                flag=0;
                while (model->window.yext_sort[j][flag][0] < discr->sub_window[0] || 
                       model->window.yext_sort[j][flag][0] >= discr->sub_window[1] ) flag++;
                model->window.yext[j][0] = y[model->window.yext_sort[j][flag][0]][j][0];

                flag=base_winsize-1;
                while (model->window.yext_sort[j][flag][1] < discr->sub_window[0] || 
                       model->window.yext_sort[j][flag][1] >= discr->sub_window[1] ) flag--;
                model->window.yext[j][1] = y[model->window.yext_sort[j][flag][1]][j][1];
            }
        }
    }
    for (j=0; j<model->OFdim; j++) {
        if (model->ofpartialdep[eqn][j]) {
            flag=0;
            while (model->window.uext_sort[j][flag][0] < discr->sub_window[0] || 
                   model->window.uext_sort[j][flag][0] >= discr->sub_window[1] ) flag++;
            model->window.uext[j][0] = u[model->window.uext_sort[j][flag][0]][j][0];
            
            flag=base_winsize-1;
            while (model->window.uext_sort[j][flag][1] < discr->sub_window[0] || 
                   model->window.uext_sort[j][flag][1] >= discr->sub_window[1] ) flag--;
            model->window.uext[j][1] = u[model->window.uext_sort[j][flag][0]][j][1];
        }
    }
    
    /* Determine the max/min time values if they are needed for monotonicity purposes */
    if (model->varpartialdep[eqn][model->DEdim]) {
        text[0] = t[discr->sub_window[0]];
        text[1] = t[discr->sub_window[1]-1];        
    }
    
    /* Determine the monotonicity of each parameter */
    for (j=0; j<model->nparms; j++) {
        if (model->parmdep[eqn][j]) (*model->determine_dependence)(eqn,pbox->p,text,model->window.yext,
                                                                   model->window.uext,model->partial,
                                                                   model,windep[2],discr,j);
        else windep[2][j] = 0;
    }
                
    /* If any of the windeps are returned as state -3, then there may be further checking involved */
    *mono = 1;
    for (i=0; i<model->nparms; i++) {
        if (windep[2][i] == -3) {
            partial[0] = 0.0; partial[1] = 0.0;
            
            for (j=discr->sub_window[0]; j<discr->sub_window[1]; j++) {
                if (model->varpartialdep[eqn][model->DEdim]) {
                    text[0] = t[j]; text[1] = t[j];
                }
                (*model->determine_dependence)(eqn,pbox->p,text,y[j],u[j],model->partial,model,windep[2],discr,i);
                partial[0] += h*discr->beta[j]*model->partial[0];
                partial[1] += h*discr->beta[j]*model->partial[1];
            }
            
            /* reset the windep if appropriate */
            if (ISZERO_INTERIOR(partial)) windep[2][i] = -2;
            else if ((temp = partial[0] + partial[1]) > 0.0) windep[2][i] = 1;
            else if (temp < 0.0) windep[2][i] = -1;
            else windep[2][i] = 0;
        }
        if(windep[2][i] == -2) *mono = 0;
    }
        
    /* Now we check monotonicity of unmeasured and (if necessary) measured variables */
    for (cor = 2; cor >= 0; cor--) {
                
        if (cor < 2) {
            for (i=0; i<model->nparms; i++) {
                p_prev[i][0] = pbox->p[i][0];
                p_prev[i][1] = pbox->p[i][1];
                if (windep[2][i]*windep[2][i] == 1) {
                    pbox->p[i][map[cor][windep[2][i]+1]] = pbox->p[i][!map[cor][windep[2][i]+1]];
                }
                windep[cor][i] = windep[2][i];
            }
        }

        /* Now check for monotonicity for unmeasured variables */
        for (j=0; j<unmsrd->num_y_unmsrd; j++) { /* for each unmeasured variable */
            k = unmsrd->y_unmsrd_order[j];
            flag = unmsrd->y_unmsrd_index[k];
            if (model->vardep[eqn][k] || eqn == k) { /* if the unmeasured is in the equation */
                for (i = discr->sub_window[0]; i < discr->sub_window[1]; i++) {
                    if (cor == 2 || windep[2][flag+i] == -2) {
                        if (discr->stencil[i] >= 2) {
                            if (model->varpartialdep[eqn][model->DEdim]) {
                                text[0] = t[i]; text[1] = t[i];
                            }
                            model->var = k;
                            (*model->determine_dependence)(eqn,pbox->p,text,y[i],u[i],model->partial,
                                                           model,windep[cor],discr,model->nparms);
                            multiply(h*discr->beta[i],model->partial);
                            alpha_flag = (eqn == k) ? discr->alpha[i] : 0;
                            
                            if (IS_INTERIOR(alpha_flag,model->partial)) windep[cor][flag+i] = -2;
                            else {
                                if (model->partial[1] > alpha_flag) windep[cor][flag+i] = 1;
                                else if (model->partial[0] < alpha_flag) windep[cor][flag+i] = -1;
                                else  windep[cor][flag+i] = 0;
                            }
                        }
                        else if (eqn == k){
                            if (discr->alpha[i] > 0) windep[cor][flag+i] = -1;
                            else if (discr->alpha[i] < 0) windep[cor][flag+i] = 1;
                            else windep[cor][flag+i] = 0;
                        }
                        else windep[cor][flag+i] = 0;              
                    }
                    else {
                        windep[cor][flag+i] = windep[2][flag+i];
                    }
                }
            }
            else {
                for (i = discr->sub_window[0]; i < discr->sub_window[1]; i++) windep[cor][flag+i] = 0;
            }
        }  
#ifdef measured_variable_monotonicity
        /* if monotonicity fails in any way, we will have to check monotonicity of the
         * meaured variables for the interval vector field calculation */
        if (!(*mono)) {
            for (j = 0; j < model->DEdim; j++) {
                if (!unmsrd->y_unmsrd[j]) { /* if variable is not unmeasured */
                    for (i = discr->sub_window[0]; i < discr->sub_window[1]; i++) {
                        if (cor == 2 || model->window.windep_msrd[2][j][i] == -2) {
                            if (discr->stencil[i] >= 2) {
                                if (model->varpartialdep[eqn][model->DEdim]) {
                                    text[0] = t[i]; text[1] = t[i];
                                }
                                model->var = j;
                                (*model->determine_dependence)(eqn,pbox->p,text,y[i],u[i],model->partial,
                                                               model,windep[cor],discr,model->nparms);
                                multiply(h*discr->beta[i],model->partial);
                                alpha_flag = (eqn == j) ? discr->alpha[i] : 0;
                                
                                if (IS_INTERIOR(alpha_flag,model->partial)) {
                                    /* put y into y_alt as is */
                                    model->window.windep_msrd[cor][j][i] = -2;
                                }
                                else {
                                    if (model->partial[1] > alpha_flag) {
                                        model->window.windep_msrd[cor][j][i] = 1;
                                    }
                                    else if (model->partial[0] < alpha_flag) {
                                        model->window.windep_msrd[cor][j][i] = -1;
                                    }
                                    else {
                                        model->window.windep_msrd[cor][j][i] = 0;
                                    }
                                }
                            }
                            else if (eqn == j){
                                if (discr->alpha[i] > 0) model->window.windep_msrd[cor][j][i] = -1;
                                else if (discr->alpha[i] < 0) model->window.windep_msrd[cor][j][i] = 1;
                                else model->window.windep_msrd[cor][j][i] = 0;
                            }
                            else model->window.windep_msrd[cor][j][i] = 0;
                        }
                        else {
                            model->window.windep_msrd[cor][j][i] = model->window.windep_msrd[2][j][i];
                        }                        
                    }
                }
            }
        }
#endif           
        if (cor < 2) { /* copy back the parameter values */
            for (i=0; i<model->nparms; i++) {
                pbox->p[i][0] = p_prev[i][0];
                pbox->p[i][1] = p_prev[i][1];
            }
        }    
        else {
            for (i=0; i<model->tparms; i++) {
                windep[0][i] = windep[2][i];
                windep[1][i] = windep[2][i];
            }
        }
    }
    
    /* If monotonicity with respect to any unmeasured variables is corner dependent, then
     * flag this information for future use */
    model->window.windep_cor_flag = 0;
    for (j=0; j<unmsrd->num_y_unmsrd; j++) {
        k = unmsrd->y_unmsrd_order[j];
        flag = unmsrd->y_unmsrd_index[k];
        for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) {
            if (windep[0][flag+i] != windep[1][flag+i]) {
                model->window.windep_cor_flag = 1; break;
            }
        }
    }
#ifdef measured_variable_monotonicity
    if(!model->window.windep_cor_flag) {
        for (j=0; j<model->DEdim; j++) {
            if (!unmsrd->y_unmsrd[j]) {
                for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) {
                    if (model->window.windep_msrd[0][j][i] != model->window.windep_msrd[1][j][i]) {
                        model->window.windep_cor_flag = 1; break;
                    }
                }
            }
        }
    }
#endif
//    if(model->window.windep_cor_flag) printf("monotonicity change occured\n");
    
//    for (i=0; i<model->tparms; i++) {
//        if (i!=3){
//            windep[0][i] = windep[2][i];
//            windep[1][i] = windep[2][i];            
//        }
//    }

//    if (model->window.windep_cor_flag == 1){
//        printf("windep[][3] = %d,%d,%d\n",windep[0][3],windep[1][3],windep[2][3]);
//        if (windep[0][3] != windep[2][3] || windep[1][3] != windep[2][3])
//        exit(1);
//    }
//    if (model->window.windep_cor_flag == 1) {
//        printf("WINDEP_COR_FLAG == 1\n");
//        
//        for (i=0; i<model->nparms; i++) {
//            printf("windep(%d) = %d\t%d\t%d\n",i,windep[0][i],windep[1][i],windep[2][i]);
//        }
//
//        for (j=0; j<unmsrd->num_y_unmsrd; j++) {
//            k = unmsrd->y_unmsrd_order[j];
//            flag = unmsrd->y_unmsrd_index[k];
//            for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) {
//                printf("windep(%d) = %d\t%d\t%d\n",flag+i,windep[0][flag+i],windep[1][flag+i],windep[2][flag+i]);
//            }
//        }
//
//        for (j=0; j<model->DEdim; j++) {
//            if (!unmsrd->y_unmsrd[j]) {
//                for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) {
//                    printf("windep(%d,%d) = %d\t%d\t%d\n",j,i,model->window.windep_msrd[0][j][i],model->window.windep_msrd[1][j][i],model->window.windep_msrd[2][j][i]);
//                }
//            }
//        }
//        exit(1);
//    }
    return;
}

double select_window (struct s_set *set, double winend_prev, int prev_progress, int winsize, 
                      double *winstart, double *h, double overlap) {
    
    /* Returns the length of the window chosen. */
//    int loc,steepness;
    double r1,r2,temp,winlength;
    
    /* Depending on the usage of the window... */
    switch (prev_progress) {
        case -2:
            /* we are in the 'settling' phase of window selection */
            winlength = *h*(winsize-1);
            *winstart = winend_prev - overlap*winlength; 
            if(*winstart + winlength > set->t_range[1]) *winstart = set->t_range[0];
            break;
        case 3:
            /* this is a first window selection */
            winlength = set->t_range[0]*(winsize-1);
            r1 = (double)rand()/(double)RAND_MAX;
            *winstart = MIN(set->t_range[1]-winlength,r1*SPAN(set->t_range)+set->t_range[0]);
            break;
        default:
            /* this is a random selection after the settling phase */
            r1 = (double)rand()/(double)RAND_MAX;
            r2 = (double)rand()/(double)RAND_MAX;
            temp = sqrt(-2*log(r1))*cos(6.28318531*r2);
            
            winlength = *h*(winsize-1);
            *winstart = winend_prev + (*h)*temp;
            if (*winstart + winlength > set->t_range[1]) {
                if (*winstart + 0.75*winlength > set->t_range[1]) *winstart = fabs(temp)*(*h);
                else *winstart = set->t_range[1]-winlength;
            } 
            break;
    }
    *winstart = MAX(set->t_range[0],*winstart);
    
//    steepness = 0;
//    loc = (int) ((*winstart-set->t_range[0])/(set->step_size[0])); 
//    while ((temp = set->t_range[0]+loc*set->step_size[0]) < set->t_range[1] && 
//           temp < *winstart+winsize*set->step_size[0]) { /* may want to use set->step_size[1] here */
//        steepness = MAX(steepness,set->change[loc]); loc++;
//    }

    /* Select a random step size between the user-defined values */
//    switch (steepness) {
//        case 0: /* flat */
//            *h = set->step_size[1];
//            break;
//        case 1: /* moderately flat */
//            *h = set->step_size[0] + ((double)rand()/(double)RAND_MAX)*SPAN(set->step_size);
//            *h = 0.5*(set->step_size[0]+set->step_size[1]) + ((double)rand()/(double)RAND_MAX)*0.5*SPAN(set->step_size);
//            *h = set->step_size[1];
//            break;
//        case 2: /* moderately steep */
//            *h = 0.5*(set->step_size[0]+set->step_size[1]) + ((double)rand()/(double)RAND_MAX)*0.5*SPAN(set->step_size);
//            *h = set->step_size[0] + ((double)rand()/(double)RAND_MAX)*0.5*SPAN(set->step_size);
//            *h = set->step_size[0] + ((double)rand()/(double)RAND_MAX)*SPAN(set->step_size);
//            break;
//        default: /* steep */
//            *h = set->step_size[1];
//            *h = set->step_size[0] + ((double)rand()/(double)RAND_MAX)*SPAN(set->step_size);
//            *h = set->step_size[0] + ((double)rand()/(double)RAND_MAX)*0.5*SPAN(set->step_size);
//            break;
//    }
    *h = set->step_size[0] + ((double)rand()/(double)RAND_MAX)*SPAN(set->step_size);
//    *h = set->step_size[0];
//    *h = set->step_size[1];
    
    *winstart = MIN(*winstart,set->t_range[1]-*h*(winsize-1));
    
#ifdef db3
	printf("        select_window returning winlength = %lf\n",winlength);
#endif
    return winlength;
}

int process_equation (struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_unmsrd *unmsrd,
                      struct s_pbox *pbox, int eqn, double y_minred, double *t, Range **y, Range **u, 
                      double h, int mono, int mono2, intRange *attempt) {
    /* return value:  -1 if there is an inconsistency, 
     *                 0 if no update was made, 
     *                 1 if an insignificant update was made. 
     *                 2 if a significant update was made. */
        
    int cons_progress,hull_progress,i,progress=0,*windep[3],winsize;
    for (i=0; i<3; i++) windep[i] = model->window.windep[i];
    struct s_discr *discr;
    double *p_ord[2];
    p_ord[0] = model->window.p_ord0;
    p_ord[1] = model->window.p_ord1;

    /* Perform 'alpha'-hull consistency */
#ifdef db2
    printf("           Calling alpha-hull_consistency\n");
#endif
    hull_progress = hull_consistency(model,discr_set_list,unmsrd,pbox,eqn,t,y,u,h,attempt,p_ord,0);

    if (hull_progress < 0) {
#ifdef db2
        printf("         Process_equation returning progress = %d\n",hull_progress);
#endif
        return -1;
    }
    
    /* Perform 'beta'-hull consistency */ 
#ifdef db2
    printf("           Calling beta-hull_consistency\n");
#endif
    hull_progress = hull_consistency(model,discr_set_list,unmsrd,pbox,eqn,t,y,u,h,attempt,p_ord,1);
    
    if (hull_progress < 0) {
#ifdef db2
        printf("         Process_equation returning progress = %d\n",hull_progress);
#endif
        return -1;
    }
    
    /* Declare the base discretization to be used for box_consistency,
     * check monotonicty, and get the sets of corner parameters */
    discr = &discr_set_list->discr_set[0].discr[discr_set_list->discr_set[0].n_discr-1];
    winsize = discr->winsize;
    discr->sub_window[0] = 0;
    discr->sub_window[1] = discr->winsize;
    check_monotonicity(model,discr,unmsrd,pbox,eqn,t,y,u,h,winsize,windep,&mono,&mono2);
    get_put_ordered_parms(0,model->tparms-1,pbox,windep[0],p_ord);
    

    /* Perform the box consistency test on each corner. */
#ifdef db2
    printf("           Calling parm_box_consistency for corner = 0\n");
#endif
//    if (mono && mono2) {
    cons_progress = parm_box_consistency(model,discr,unmsrd,pbox,0,p_ord[0],p_ord[1],eqn,t,y,u,
                                         windep[0],h,attempt,mono); /* mono replaced */

    if (cons_progress < 0) {
#ifdef db2
        printf("         Process_equation returning progress = %d\n",progress);
#endif
        return -1;
    }
    else if (cons_progress == 1) progress = 1;
    
    if (model->window.windep_cor_flag) {
        get_put_ordered_parms(1,model->tparms-1,pbox,windep[0],p_ord);
        get_put_ordered_parms(0,model->tparms-1,pbox,windep[1],p_ord);
    }
    
#ifdef db2
    printf("           Calling parm_box_consistency for corner = 1\n");
#endif
    cons_progress = parm_box_consistency(model,discr,unmsrd,pbox,1,p_ord[1],p_ord[0],eqn,t,y,u,
                                         windep[1],h,attempt,mono); /* mono replaced */
    

    if (cons_progress < 0) {
#ifdef db2
        printf("         Process_equation returning progress = %d\n",cons_progress);
#endif
        return -1;
    }
    else if (cons_progress == 1) progress = 1;
        
    /* Re-insert the ordered parameters back into the ranges. */
    get_put_ordered_parms(1,model->nparms+winsize*unmsrd->num_y_unmsrd-1,pbox,windep[1],p_ord);

//    }
//    else {
//        cons_progress = 0;
//    }
   
    
#ifdef db2
	printf("         Process_equation returning progress = %d\n",cons_progress);
#endif
    return progress;
}

int process_algebraic_equation (struct s_model *model, struct s_discr_set_list *discr_set_list,
                                struct s_unmsrd *unmsrd, struct s_pbox *pbox, int eqn, double y_minred,  
                                double *t, Range **y, Range **u, int mono, int mono2, 
                                intRange *attempt) {
        
    int cons_progress,i,j,k,h=1.0,progress=0,*windep[3];
    for (i=0; i<3; i++) windep[i] = model->window.windep[i+3];
    double discr_error_temp = model->discr_error,*p_ord[2];
    struct s_discr *discr;
    p_ord[0] = model->window.p_ord0;
    p_ord[1] = model->window.p_ord1;
    model->discr_error = 0.0;
        
    /* use a set-up fake stencil so that we can utilize pre-existing code for box consistency */
    discr = &discr_set_list->discr_set[3].discr[0];
    h = 1.0;
    for (i=0; i<unmsrd->num_y_unmsrd*discr->winsize; i++) {
        for (j=0; j<3; j++) windep[j][model->nparms+i] = 0;
    }   

    /* Perform a regular parameter reduction for each time step on the algebraic equation. */
    for (i=0; i<discr->winsize; i++) {
        discr->sub_window[0] = i;
        discr->sub_window[1] = i+1;
        
        /* Check monotonicity on this time step */
        check_monotonicity(model,discr,unmsrd,pbox,eqn,t,y,u,h,discr->winsize,windep,&mono,&mono2);
        get_put_ordered_parms(0,model->nparms+discr->winsize*unmsrd->num_y_unmsrd-1,pbox,windep[0],p_ord);
        
        /* Reset windeps for time steps not being used */
        for (j=0; j<unmsrd->num_y_unmsrd; j++)
            if (i > 0) for(k=0; k<3; k++) windep[k][unmsrd->y_unmsrd_index[unmsrd->y_unmsrd[j]]+i-1] = 0;
        /* TEST THIS */
        
#ifdef db2
        printf("           Calling parm_box_consistency for corner=%d.\n",0);
#endif
        cons_progress = parm_box_consistency(model,discr,unmsrd,pbox,0,p_ord[0],p_ord[1],eqn,t,y,u,
                                             windep[0],h,attempt,mono); /* mono replaced */
        if (cons_progress < 0) { progress = -1; break; }
        else progress = cons_progress;
        
        if (model->window.windep_cor_flag) {
            get_put_ordered_parms(1,model->tparms-1,pbox,windep[0],p_ord);
            get_put_ordered_parms(0,model->tparms-1,pbox,windep[1],p_ord);
        }

#ifdef db2
        printf("           Calling parm_box_consistency for corner=%d.\n",1);
#endif
        cons_progress = parm_box_consistency(model,discr,unmsrd,pbox,1,p_ord[1],p_ord[0],eqn,t,y,u,
                                             windep[1],h,attempt,mono); /* mono replaced */
        if (cons_progress < 0) { progress = -1; break; }
        else progress = MAX(progress,cons_progress);
        
        /* put back the updated parameters */
        get_put_ordered_parms(1,model->nparms+discr->winsize*unmsrd->num_y_unmsrd-1,pbox,windep[1],p_ord);
    }
    
    model->discr_error = discr_error_temp;
    return progress;
}

void get_put_ordered_parms (int flag, int tparms, struct s_pbox *pbox, int *windep, double *p_ord[2]) {
    
    /* windep has an entry of 0 if the eqn on this window does not depend on 
     * the parameter, 1 means the eqn is monotonically non-decreasing, and 
     * -1 means the eqn is monotonically nonincreasing. p_ord[0] will be the parameters
     * at the low corner and p_ord[1] at the upper corner. */
    
    int i;
    static int map[2][3] = {
		{ 1, 0, 0},
		{ 0, 0, 1}};  /* only entries in the first and third columns are relevant*/
    	
    if (flag == 0) { /* get */
        for (i=tparms; i>=0; i--) {
            if (windep[i] == -2) {
                p_ord[0][i] = pbox->p[i][0];
                p_ord[1][i] = pbox->p[i][1];
            }
            else if (windep[i] != 0) {
                p_ord[0][i] = pbox->p[i][map[0][windep[i]+1]];
                p_ord[1][i] = pbox->p[i][map[1][windep[i]+1]];  
            }
		}
    }
    else { /* flag == 1, put */
        for (i=tparms; i>=0; i--) {
            if (windep[i] == -2) {
                pbox->p[i][0] = p_ord[0][i];
                pbox->p[i][1] = p_ord[1][i];
            }
            else if (windep[i] != 0) {
                pbox->p[i][map[0][windep[i]+1]] = p_ord[0][i];
                pbox->p[i][map[1][windep[i]+1]] = p_ord[1][i];  
            }
		}
    }
    return;
}

int hull_consistency (struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_unmsrd *unmsrd,
                      struct s_pbox *pbox, int eqn, double *t, Range **y, Range **u,
                      double h, intRange *attempt, double *p_ord[2], int flag) {
    
    /* The value of flag indicates the purpose of the function call 
     * (flag == 0): alpha-hull consistency, (flag == 1): beta-hull consistency */
    int base_winsize,get_put_flag,hull_flag,i,index,j,j_max,j_min,k,mono,mono2,progress=0;
    int use_discr,*windep[3];
    for (i=0; i<3; i++) windep[i] = model->window.windep[3+i];
    struct s_discr *discr;
    double divisor,F_temp;
    Range F,*y_alt=model->window.y_alt1;
    intRange *subwindow;
    static int map[2][4] = {
		{ 0, 1, 0, 1},
		{ 0, 0, 1, 1}};

    /* Depending on which hull consistency being performed, different
     * subwindow parameters are called from the input file.*/
    if (flag == 0) subwindow = &discr_set_list->sub_win_var;
    else subwindow = &discr_set_list->sub_win_parm;
    i = discr_set_list->discr_set[0].discr[discr_set_list->discr_set[0].n_discr-1].winsize;
    base_winsize = ((*subwindow)[0] != -1)*i;
    j_min = (*subwindow)[0] - 2;
    j_max = (*subwindow)[1] - 1;
    
    /* Loop through all subwindows */
    for (i=1; i<base_winsize; i++) {
        for (j=MAX(0,i-j_max); j<i-j_min; j++) { 
            
            /* Choose the appropriately size discretization, depending on the purpose */
            if (flag == 0) use_discr = 1;
            else use_discr = 2;
            discr = &discr_set_list->discr_set[use_discr].discr[i-j-1];
            
            /* If we are performing an 'alpha'-hull consistency, then only proceed if
             * the end time step is variable dependent. If we are performing a 'beta'-hull
             * consistency, only proceed if the end time step is vector field-dependent. */
            if (map[flag][discr->stencil[i-j]]) {
                
                discr->sub_window[0] = j;
                discr->sub_window[1] = i+1;
                
                /* Shift the discretization values into their correct places */
                for (k=i; k>=j; k--) {
                    discr->alpha[k] = discr->alpha[k-j];
                    discr->beta[k] = discr->beta[k-j];
                    discr->stencil[k] = discr->stencil[k-j];
                }
                if (flag == 0) {
                    divisor = discr->alpha[i];
                    discr->alpha[i] = 0.0;
                    discr->stencil[i] -= 1;
                }
                else {
                    divisor = -h*discr->beta[i];
                    discr->beta[i] = 0.0;
                    discr->stencil[i] -= 2;
                }

                /* Check monotonicity given the subwindow and stencil */
                check_monotonicity(model,discr,unmsrd,pbox,eqn,t,y,u,h,base_winsize,windep,&mono,&mono2);
                
                if (unmsrd->num_y_unmsrd > 0)
                    get_put_flag = model->nparms+(unmsrd->num_y_unmsrd-1)*base_winsize+i;
                else
                    get_put_flag = model->nparms-1;
                
                get_put_ordered_parms(0,get_put_flag,pbox,windep[0],p_ord);
                
                /* Determine the hull of the remaining discretization components, and 
                 * divide to create the hull consistency compare interval. */
                F[0] = constraint_eval(model,discr,unmsrd,0,-model->nparms,p_ord[0],p_ord[1],eqn,
                                       t,y,u,h,mono,windep[0]); /* mono replaced */
                
                if (model->window.windep_cor_flag) {
                    get_put_ordered_parms(1,get_put_flag,pbox,windep[0],p_ord);
                    get_put_ordered_parms(0,get_put_flag,pbox,windep[1],p_ord);
                }
                
                F[1] = constraint_eval(model,discr,unmsrd,1,-model->nparms,p_ord[1],p_ord[0],eqn,
                                       t,y,u,h,mono,windep[1]); /* mono replaced */
                k = unmsrd->y_unmsrd_index[eqn];

                if (divisor > 0) {
                    F[0] = (F[0] - model->discr_error)/divisor;
                    F[1] = (F[1] + model->discr_error)/divisor;
                }
                else {
                    F_temp = F[0];
                    F[0] = (F[1] + model->discr_error)/divisor;
                    F[1] = (F_temp - model->discr_error)/divisor;
                }

                if (flag == 0) {
                    /* If we are performing an 'alpha'-hull consistency and an unmeasured
                     * variable is in the derivative, then make an intersection and update
                     * as needed. If we are in a different equation, then check to see if
                     * we can discard the box with an inconsistency. */
                    if (unmsrd->y_unmsrd[eqn]) { 
                        /* Determine if the variable corresponding to the derivative in 
                         * this equation is unmeasured. If so, determine the index in which 
                         * these variable parameter will be located in p_ord, pbox */
                        k = unmsrd->y_unmsrd_index[eqn];
                        windep[1][k+i] = 1;
                        p_ord[0][k+i] = pbox->p[k+i][0];
                        p_ord[1][k+i] = pbox->p[k+i][1];
#ifdef db3
                        printf("              on sub_window {%d,%d} intersection of F=[%f,%f], p=[%f,%f]\n",j,i,F[0],F[1],p_ord[0][k+i],p_ord[1][k+i]);
#endif
                        if (F[0] > p_ord[1][k+i] || F[1] < p_ord[0][k+i]) {
                            progress = -1;
                        }
                        else {
                            p_ord[0][k+i] = MAX(F[0],p_ord[0][k+i]);
                            p_ord[1][k+i] = MIN(F[1],p_ord[1][k+i]);      
                        }
                    }
                    else {
#ifdef db3
                        printf("              on sub_window {%d,%d} intersection of F=[%f,%f], y=[%f,%f]\n",j,i,F[0],F[1],y[i][eqn][0],y[i][eqn][1]);
#endif
                        if (F[0] > y[i][eqn][1] || F[1] < y[i][eqn][0]) {
                            progress = -1;
                        }
                        else {
                            y[i][eqn][0] = MAX(y[i][eqn][0],F[0]);
                            y[i][eqn][1] = MIN(y[i][eqn][1],F[1]);
                            model->window.y_alt[eqn][i][0] = y[i][eqn][0];
                            model->window.y_alt[eqn][i][1] = y[i][eqn][1];
                        }
                    }
                }
                else {
                    /* If we are conducting a 'beta'-hull consistency, then update any
                     * unmeasured variables while transferring the data into the 
                     * y_alt array (sorted differently for each of interval manipulation) */
                    if (unmsrd->num_y_unmsrd > 0) {
                        for (k = 0; k < model->DEdim; k++){
                            if (model->vardep[eqn][k] || eqn == k) {
                                if(unmsrd->y_unmsrd[k]){ /* if variable is unmeasured */
                                    index = unmsrd->y_unmsrd_index[k]; /* TEST this */
                                    y_alt[k][0] = MIN(p_ord[0][index+i],p_ord[1][index+i]);
                                    y_alt[k][1] = MAX(p_ord[0][index+i],p_ord[1][index+i]);
                                }
                                else{
                                    y_alt[k][0] = y[i][k][0];
                                    y_alt[k][1] = y[i][k][1];
                                }
                            }
                        }
                    }
                    
                    /* For each parameter, if it is included in the vector field, then attempt to 
                     * invert. If the reuturn of hull_flag is 0, then there was nothing coded by
                     * the user for inverting the vector field. If the return is 1, then some progress
                     * was made, if the return is -1, there was an inconsistency */                    
                    for (k = 0; k < model->nparms; k++) {
                        if (model->parmdep[eqn][k]) {
                            if(unmsrd->num_y_unmsrd > 0){
                                hull_flag = (*model->invert_vec_field)(eqn,k,F,model,pbox->p,t[i],y_alt,u[i]);
                            }
                            else{
                                hull_flag = (*model->invert_vec_field)(eqn,k,F,model,pbox->p,t[i],y[i],u[i]);
                            }
                            if (hull_flag == -1) {  progress = -1; break; }
                        }
                    }
                    /* If the model file includes an inversion for the unmeasured variable */
                    if (unmsrd->num_y_unmsrd > 0 && progress >= 0) {
                        for (k = 0; k < unmsrd->num_y_unmsrd; k++) {
                            if (model->vardep[eqn][unmsrd->y_unmsrd_order[k]]) {
                                index = model->nparms+k*base_winsize;
                                model->var = unmsrd->y_unmsrd_order[k];
                                hull_flag = (*model->invert_vec_field)(eqn,index+i,F,model,pbox->p,t[i],y_alt,u[i]);
                                if (hull_flag == -1) {  progress = -1; break; }
                            }   
                        }
                    }
                }
                
                /* Reset the discretization. */
                if (flag == 0) {
                    discr->alpha[i] = divisor;
                    discr->stencil[i] += 1;                    
                }
                else {
                    discr->beta[i] = divisor/(-h);
                    discr->stencil[i] += 2;
                }
                for (k=0; k<=i-j; k++) {
                    discr->alpha[k] = discr->alpha[k+j];
                    discr->beta[k] = discr->beta[k+j];
                    discr->stencil[k] = discr->stencil[k+j];
                }
                
                if (progress < 0){  i = base_winsize; break; }     
                
                get_put_ordered_parms(1,get_put_flag,pbox,windep[1],p_ord);
            }
        }
    }
#ifdef db2
	printf("           hull_consistency returning progress = %d\n",progress);
#endif
    return progress;
}

int parm_box_consistency (struct s_model *model, struct s_discr *discr, struct s_unmsrd *unmsrd,
                          struct s_pbox *pbox, int cor, double *p_cor, double *p_opp, int eqn,  
                          double *t, Range **y, Range **u, int *windep, double h, intRange *attempt, int mono) {
    /* return value:  -1 if there is an inconsistency 
     *	               0 if no update was made, 
     *	               1 if an update was made. */
            
    static double pm1[2] = {1.0, -1.0};
    static int map[2][3] = {
		{ 0, 0, 1},
		{ 1, 0, 0}};  /* only entries in the first and third columns are relevant*/
    int base_winsize=unmsrd->base_winsize,bisect,ctr,i,index,j;
    int opp_cor = (cor == 0) ? 1 : 0,progress=0,time_step,variable;
    double Fest[2],Fext,*p_adj=model->window.p_ord2,pbracket_old[2],pbracket[2],pest[2];
    double target,temp,temp1,temp2;

    /* update any unmeasured variables if necessary. */
    for (j=0; j<unmsrd->num_y_unmsrd; j++) {
        variable = unmsrd->y_unmsrd_order[j];
        index = unmsrd->y_unmsrd_index[variable];
        if (model->vardep[eqn][unmsrd->y_unmsrd_order[j]] || eqn==variable) {
            if (mono) {
                for (i = 0; i < discr->winsize; i++) {
                    y[i][variable][0] = p_cor[index + i];
                    y[i][variable][1] = p_cor[index + i];
                }
            }
            else {
                for (i = 0; i < discr->winsize; i++) {
                    if (windep[index+i] == -2) {
                        model->window.y_alt[variable][i][cor] = p_cor[index + i];
                        model->window.y_alt[variable][i][opp_cor] = p_opp[index + i];
                    }
                    else if (windep[index+i] != 0) {
                        model->window.y_alt[variable][i][0] = p_cor[index + i];
                        model->window.y_alt[variable][i][1] = p_cor[index + i];
                    }
                }
            }
        }
    }
    
    /* Calculate the constraint function at the extreme corner. */
    Fext = constraint_eval(model,discr,unmsrd,cor,-model->nparms,p_cor,p_opp,eqn,t,y,u,h,mono,windep); 

#ifdef db3
	printf("              Fext = %g\n",Fext);
#endif

    /* Check for inconsistency. */
    if (Fext*pm1[cor] > 0.0) {
#ifdef db2
		printf("                INCONSISTENCY parm_box Fext = %g, t = %g\n",Fext,t[0]);
#endif        
//        for (i=0; i<model->tparms; i++) {
//            printf("p[%d] = %f,%f\n",i,p_cor[i],p_opp[i]);
//        }
//        printf("Fext = %f, eqn %d, t = %f\n",Fext,eqn,t[0]);
		return -1;
    }

    memcpy(p_adj,p_cor,model->tparms*sizeof(double));

    /* Loop through the parameters doing box consistency checks. */   
    for (i=0; i<model->tparms; i++) {
                
        /* If the discretization is non-monotonic with respect the parameter */
        if (windep[i] == -2) {
            /* to have made it this far in the code, we must have been consistent, and we
             * know that the extreme value occured at some interior parameter value.
             * Now we must check both end points to see if we cross 0 at the extremes. 
             * Start on the LEFT (lower) end of the parameter interval. */
            if (attempt[i][cor]) {
                
                temp1 = p_cor[i];
                temp2 = p_opp[i];
                p_adj[i] = p_cor[i];
                p_opp[i] = p_cor[i];
                
                /* Calculate the edge value */
                Fest[1] = constraint_eval(model,discr,unmsrd,cor,i-model->nparms,p_adj,p_opp,eqn,
                                          t,y,u,h,mono,windep);
                
#ifdef influence_split
                if (fabs(Fest[1]-Fext) < fabs(Fext)) { /* update influence */
                    temp = MIN(fabs(Fest[1]-Fext)/fabs(Fext),1.0);
                    pbox->inf[eqn][MIN(i,model->nparms)].pct[opp_cor] += temp;
                    pbox->inf[eqn][MIN(i,model->nparms)].ctr[opp_cor]++;
                }
#endif
                target = MAX(model->discr_error*0.75,0.001)*pm1[cor];
                p_opp[i] = temp2;
                
                /* If we pass the discretization error target */
                if (Fest[1]*pm1[cor] > 1.25*fabs(target)) {
//                    if (i > 6 && cor == 1) printf("%f\t%f\n",Fext,Fest[1]);

                    pbracket[0] = p_cor[i];
                    pbracket_old[1] = p_opp[i];
                    pbracket[1] = 0.5*(p_cor[i]+p_opp[i]);
                    
                    p_adj[i] = p_cor[i];
                    p_opp[i] = pbracket[1];
                    
                    Fest[0] = constraint_eval(model,discr,unmsrd,cor,i-model->nparms,p_adj,p_opp,eqn,
                                              t,y,u,h,mono,windep);
                    ctr = 0;

                    while (fabs(Fest[0]-target) > 0.25*fabs(target) && ++ctr <= 20) {
                        /* Adjust the bracket until the entire interval result of the 
                         * constraint eval lies entirely past the target */
                        if (IS_INTERIOR(target,Fest)) {
                            pbracket_old[1] = pbracket[1];
                            pbracket[1] = 0.5*(pbracket[0]+pbracket[1]);
                        }
                        else {
                            if (Fest[0]*pm1[cor] > fabs(target)) pbracket[0] = pbracket[1];
                            pbracket[1] = 0.5*(pbracket[1]+pbracket_old[1]);
                        }
                        p_opp[i] = pbracket[1];
                        Fest[0] = constraint_eval(model,discr,unmsrd,cor,i-model->nparms,p_adj,p_opp,eqn,
                                                  t,y,u,h,mono,windep);
                    }
                    progress = 1;
                    p_cor[i] = p_opp[i];
                }
                else{
                    p_cor[i] = temp1;
                }
                p_opp[i] = temp2;
#ifdef db3
                printf("                  windep = %d, cor = %d, Updating p[%d] = %f %f, to %f %f\n",windep[i],cor,i,temp1,temp2,p_cor[i],p_opp[i]);
#endif
            }   
            if (attempt[i][opp_cor]) {
                /* Now, proceed on the RIGHT (upper) end of the parameter interval. */
                temp1 = p_cor[i];
                temp2 = p_opp[i];
                
                p_adj[i] = p_opp[i];
                Fest[1] = constraint_eval(model,discr,unmsrd,cor,i-model->nparms,p_adj,p_opp,eqn,
                                          t,y,u,h,mono,windep);
#ifdef influence_split               
                if (fabs(Fest[1]-Fext) < fabs(Fext)) { /* update influence */
                    temp = MIN(fabs(Fest[1]-Fext)/fabs(Fext),1.0);
                    pbox->inf[eqn][MIN(i,model->nparms)].pct[cor] += temp;
                    pbox->inf[eqn][MIN(i,model->nparms)].ctr[cor]++;
                }
#endif         
                target = MAX(model->discr_error*0.75,0.001)*pm1[cor];
                
                p_opp[i] = temp2;
                
                /* If we pass the discretization error target */
                if (Fest[1]*pm1[cor] > 1.25*fabs(target)) {            
                    pbracket_old[0] = p_cor[i];
                    pbracket[0] = 0.5*(p_cor[i]+p_opp[i]);
                    pbracket[1] = p_opp[i];
                    
                    p_adj[i] = pbracket[0];
                    Fest[0] = constraint_eval(model,discr,unmsrd,cor,i-model->nparms,p_adj,p_opp,eqn,
                                              t,y,u,h,mono,windep);
                    ctr = 0;
                    while (fabs(Fest[0]-target) > 0.25*fabs(target) && ++ctr <= 20) {
                        /* Adjust the bracket until the entire interval result of the 
                         * constraint eval lies entirely past the target */
                        if (IS_INTERIOR(target,Fest)) {
                            pbracket_old[0] = pbracket[0];
                            pbracket[0] = 0.5*(pbracket[0]+pbracket[1]);
                        }
                        else {
                            if (Fest[0]*pm1[cor] > fabs(target)) pbracket[1] = pbracket[0];
                            pbracket[0] = 0.5*(pbracket[0]+pbracket_old[0]);
                        }
                        p_adj[i] = pbracket[0];
                        Fest[0] = constraint_eval(model,discr,unmsrd,cor,i-model->nparms,p_adj,p_opp,eqn,
                                                  t,y,u,h,mono,windep);
                    }

                    progress = 1;
                    p_opp[i] = p_adj[i];
                }
                else{
                    p_opp[i] = temp2;
                }
                p_cor[i] = temp1;    
                p_adj[i] = p_cor[i];
#ifdef db3
                printf("                  windep = %d cor = %d, Updating p[%d] = %f %f, to %f %f\n",windep[i],cor,i,temp1,temp2,p_cor[i],p_opp[i]);
#endif
            }
        }      
		else if (windep[i] != 0 && attempt[i][map[cor][windep[i]+1]]) {     
            /* Otherwise, the discretization is monotonic with respect to this parameter,
             * so we know that the extreme occurs at the appropriate end point */
			Fest[0] = Fext;
            
            /* Calculate the constraint function at the adjacent corner. */
			p_adj[i] = p_opp[i]; 
			Fest[1] = constraint_eval(model,discr,unmsrd,cor,i-model->nparms,p_adj,p_opp,eqn,
                                      t,y,u,h,mono,windep);
#ifdef db3
			printf("                Fadj for parm %d = %f\n",i,Fest[1]);
#endif
#ifdef influence_split
            if (fabs(Fest[1]-Fext) < fabs(Fext)) { /* update influence */
                temp = MIN(fabs(Fest[1]-Fext)/fabs(Fext),1.0);
                pbox->inf[eqn][MIN(i,model->nparms)].pct[map[opp_cor][windep[i]+1]] += temp;
                pbox->inf[eqn][MIN(i,model->nparms)].ctr[map[opp_cor][windep[i]+1]]++;                    
            }
#endif
			/* If the function crosses +-[0.0,(1.0+AIM)*Fderr], then a 
			 * reduction in the range is possible. Set the target value */
            target = MAX(model->discr_error/0.75,0.001)*pm1[cor];
            
            /* If we pass the target value */
			if (Fest[1]*pm1[cor] > model->discr_error) {

				pbracket[0] = p_cor[i];
				pbracket[1] = p_adj[i];
                pest[0] = pbracket[0];
                pest[1] = pbracket[1];
				ctr = 0;
				bisect = 0;
				/* Use a secant/bisection method to find the parameter value at
				 * which the F value is in target*[0.5,1.5].  Use secant
				 * unless the guess puts you outside the current bracket around the
				 * target value, in which case do a bisection of the bracket.  Only 
				 * allow at most 20 iterations. */
				while (fabs(Fest[1] - target) > 0.25*fabs(target) && ++ctr <= 20) {
					p_adj[i] = pest[1] + (target - Fest[1])*SPAN(pest)/SPAN(Fest);
                    
					if (!IS_INTERIOR(p_adj[i],pbracket)) { /* Do a bisection. */
						p_adj[i] = 0.5*(pbracket[0] + pbracket[1]);
						bisect = 1;
					}
					temp = constraint_eval(model,discr,unmsrd,cor,i-model->nparms,p_adj,p_opp,eqn,
                                           t,y,u,h,mono,windep);                    
#ifdef db3
					printf("                  pbracket = [%g %g]\n",pbracket[0],pbracket[1]);
					printf("                    new p_adj = %g,  new Fest[1] = %g\n",p_adj[i],temp);
#endif
					/* j will be one if the p_adj[i] value is still beyond
					 * (from the extreme corner) to where we are aiming and zero
					 * otherwise. j is the pbracket end we will update. */
					j = temp*pm1[cor] > target*pm1[cor];
					if (bisect) {
						bisect = 0;
						/* The new estimates for re-starting the secant method will
						 * be (p_adj[i],temp) and one of (pest[0],Fest[0]) or
						 * (pest[1],Fest[1]).   But (p_adj[i],temp) will be put into
						 * the (pest[1],Fest[1]) spot so that the while test clause
						 * will be correct.  At least one of pest[0] or pest[1] 
						 * must be equal to either (old) pbracket[0] or pbracket[1].
						 * If only one is, choose it, otherwise choose the one 
						 * that will also be equal to the endpoint of the new 
						 * pbracket that will not be updated. */
						if (pbracket[!j] == pest[1] || (pbracket[!j] != pest[0] && pbracket[j] == pest[1])) {
                            pest[0] = pest[1];
                            Fest[0] = Fest[1];
                        }
						pest[1] = p_adj[i];
						Fest[1] = temp;
					}
					else {
						pest[0] = pest[1];
						pest[1] = p_adj[i];
						Fest[0] = Fest[1];
						Fest[1] = temp;
					}
					/* Update the bracket. */
					pbracket[j] = p_adj[i];
				}
				if (ctr == 21) { /* Above did not converge, use best bracket. */
					p_adj[i] = pbracket[1];
#ifdef db3
					printf("Hit the limit.\n");
#endif
				}
				/* Update. */	
#ifdef db3
				printf("                  bc ctr=%d cor = %d, Updating p[%d] = %f %f, to %f %f\n",ctr,cor,i,
                       p_cor[i],p_opp[i],p_cor[i],p_adj[i]);
#endif
				p_opp[i] = p_adj[i];
				progress = 1;
			}
			/* Reset p_adj[i] to the corner value so that the next adjacent corner
			 * will have the correct p_adj values. */
			p_adj[i] = p_cor[i];            
		}
        
        /* Reset y,y_alt if this is an unmeasured parameter */
        if (i >= model->nparms) {
            variable = unmsrd->y_unmsrd_order[(i-model->nparms)/base_winsize];  
            if (model->vardep[eqn][variable]) {
                time_step = (i-model->nparms)%base_winsize;
                if (mono) {
                    y[time_step][variable][0] = p_cor[i];
                    y[time_step][variable][1] = p_cor[i];
                }
                else {
                    if (windep[i] == -2) {
                        model->window.y_alt[variable][time_step][cor] = p_cor[i];
                        model->window.y_alt[variable][time_step][opp_cor] = p_opp[i];                                    
                    }
                    else if (windep[i] != 0) {
                        model->window.y_alt[variable][time_step][0] = p_cor[i];
                        model->window.y_alt[variable][time_step][1] = p_cor[i];                                    
                    }
                }
            }
        }
    }
    
#ifdef db2
	printf("           parm_box_consistency returning progress = %d\n",progress);
#endif
    return progress;
}

double constraint_eval (struct s_model *model, struct s_discr *discr, struct s_unmsrd *unmsrd, int cor, 
                        int um_parm, double *p, double *p_opp, int eqn, double *t, Range **y, Range **u, 
                        double h, int mono, int *windep) {
	    
    int i,time_step,variable;
    double F=0.0,betadt;
    static int map[2][2] = {
		{ 0, 1},
		{ 1, 0}};
#ifdef measured_variable_monotonicity
    int j;
    static int map2[2][3] = {
		{ 1, 0, 0},
		{ 0, 0, 1}};  /* only entries in the first and third columns are relevant*/
#endif
    
    if (um_parm >= 0) { /* we are processing a unmeasured variable */
        variable = unmsrd->y_unmsrd_order[um_parm/unmsrd->base_winsize];    
        time_step = um_parm%unmsrd->base_winsize;
    }
    
    if (mono) { 
        /* IF the discretization is fully monotonic, the LMS constraint function is of the form 
         *    h*SUM(beta[i]*f(y^i,t_i)) -  SUM(alpha[i]*y^i)
         * so if alpha[i] > 0 and beta[i] == 0, the function is decreasing with y^i.
         * Then map[alpha[i] > 0.0][cor] will pick out the appropriate range end for y^i.
         * The vec_field function evaluates f(y^i,t_i) + a*y^i, so if beta[i] != 0, then
         * the contribution to the constraint is
         *    h*beta[i]*vec_field(...,a,...) , where a = -alpha[i]/(beta[i]*h).  */

        /* Update unmeasured variables if necessary */
        if (um_parm >= 0) {
            /* Which unmeasured variable are we currently processing as a parameter */
            if (model->vardep[eqn][variable]) {
                y[time_step][variable][0] = p[model->nparms+um_parm];
                y[time_step][variable][1] = p[model->nparms+um_parm];
            }
        }
        
        /* Calculate vector field */        
        for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) {
            if (discr->stencil[i] == 1) {
                F -= discr->alpha[i]*y[i][eqn][map[discr->alpha[i] > 0.0][cor]];
            }
            else if (discr->stencil[i] > 1) {
                betadt = discr->beta[i]*h;
                F += betadt*(*model->vec_field)(cor,eqn,p,-discr->alpha[i]/betadt,t[i],y[i],u[i],model);
            }
        }
    }
    else { 
        /* if F fails monotonicity, then calculate using interval arithmetic 
         * the interval vector field accesses information in a different way, 
         * so we need to update model->window.y_alt */
        if (um_parm >= 0 && model->vardep[eqn][variable]) {
            if (windep[model->nparms+um_parm] == -2) {
                /* Extra monotonicity check here too */
                model->window.y_alt[variable][time_step][cor] = p[model->nparms+um_parm];
                model->window.y_alt[variable][time_step][!cor] = p_opp[model->nparms+um_parm];
            }
            else if (windep[model->nparms+um_parm] != 0) {
                model->window.y_alt[variable][time_step][0] = p[model->nparms+um_parm];
                model->window.y_alt[variable][time_step][1] = p[model->nparms+um_parm];  
            }
        }
        
        /* sp is used to reassemble p and p_opp into a usable array of ranges we cannot simply 
         * use pbox->p since we are using these arrays as part of a box consistency context */
        for (i=0; i<model->nparms; i++) {
            if (model->parmdep[eqn][i]) {
                if (windep[i] == -2) {
                    if (cor == 0) {
                        model->window.p_alt[i][0] = p[i];
                        model->window.p_alt[i][1] = p_opp[i];
                    }
                    else {
                        model->window.p_alt[i][0] = p_opp[i];
                        model->window.p_alt[i][1] = p[i];
                    }
                }
                else if (windep[i] != 0) {
                    model->window.p_alt[i][0] = p[i];
                    model->window.p_alt[i][1] = p[i];
                }
            }
        }

#ifdef measured_variable_monotonicity
        /* update measured variables using monotonicity information */
        for (j=0; j<model->DEdim; j++) {
            if (!unmsrd->y_unmsrd[j] && (model->vardep[eqn][j] || eqn == j)) {
                for (i = discr->sub_window[0]; i < discr->sub_window[1]; i++) {
                    if (model->window.windep_msrd[cor][j][i] == -2) {
                        model->window.y_alt[j][i][0] = y[i][j][0];
                        model->window.y_alt[j][i][1] = y[i][j][1];
                    }
                    else if (model->window.windep_msrd[cor][j][i] != 0) {
                        model->window.y_alt[j][i][0] = y[i][j][map2[cor][model->window.windep_msrd[cor][j][i]+1]];
                        model->window.y_alt[j][i][1] = model->window.y_alt[j][i][0];                        
                    }
                }
            }
        }
#endif
        
        /* compute the vector field, first by using the alphas, and then by adding
         * the interval result of the vector field itself. */
        for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) {
            if(discr->stencil[i]%2){                
                F -= discr->alpha[i]*model->window.y_alt[eqn][i][map[discr->alpha[i] > 0][cor]];
            }
        }    
        F += h*(*model->interval_vec_field)(eqn,cor,model->window.p_alt,t,model->window.y_alt,
                                            model->window.u_alt,model,discr);
    }
    return F;
}
