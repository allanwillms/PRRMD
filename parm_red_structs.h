#ifndef TypedefRange
#define TypedefRange
typedef double Range[2];
typedef int intRange[2];
typedef int intRange3[3];
#endif

#define SPAN(a) ((a)[1] - (a)[0])
#define ISZERO_INTERIOR(a) (((a)[0] * (a)[1]) < 0.0 ? 1 : 0 )
#define IS_INTERIOR(a,b) ((((b)[0] < (a) && (a) < (b)[1]) || ((b)[1] < (a) && (a) < (b)[0])) ? 1 : 0 )

struct s_yu {
    char filename[100];  /* file name which sourced this data */
    int npts;            /* number of time points in data file */
    double *input_t;     /* input time values */
    double *input;       /* array of input points, input[npts] */
    double gw;     /* goal width */
    double tmin;   /* minimum time for tasle, negative indicates data not fit by tasle */
    int ncpt;      /* number of control points */
    double *t;     /* array of t control points */
    Range *cpt;   /* array of control points */
};

struct s_unmsrd {
    int num_y_unmsrd; /* number of unmeasured variables */
    int *y_unmsrd; /* array indicating if a variable is measured (0) or unmeasured (1) */
    int *y_unmsrd_order; /* array numbering the unmeasured variables in the set */
    int *y_unmsrd_index; /* array indicating where in pbox->p each unmeasured variable parameter begins */
    int base_winsize; /* number of time steps in a default window (used for index management) */
};

struct s_set {
    Range t_range;   /* range of input t values (first and last t values) */
    Range step_size;  /* maximum discretization step size to use for this data set */
    struct s_yu *y;  /* array of y info y[DEdim] */
    struct s_yu *u;  /* array of u info u[OFdim] */
    struct s_unmsrd unmsrd;  /* contains information regarding unmeasured variables in this set */
    enum e_change *change;  /* contains information regarding relative steepness of measured variables 
                             * size of number of control points of smallest h value */
};

enum e_change { flat = 0, moderate_flat, moderate_steep, steep};
    /* determining the relative steepness of our measured variables at a specific time step */

struct s_data {
    int nsets; /* number of data sets */
    int num_y_unmsrd_max; /* maximum number of unmeasured variables among all sets */
    int *y_unmsrd_max; /* summary over all sets of whether a variable is ever unmeasured */
    int write_vars; /* flag for writing variables at output */
    double minred; /* mininum reduction fraction (of goal width) for success */
    struct s_set *set;  /* array of sweeps */
};

enum e_box_status { invalid = -1, putative, putative_split, good};
    /*
     * invalid	    : parameter ranges are invalid 
     * putative	    : parameter ranges not shown to be invalid
     * putative_split  : needs to be split for further progress 
     * good	    : parameter goals have been met */

struct s_inf {
    Range pct; /* Tracks the percentage influence of each variable/equation pair */
    intRange ctr; /* Counts how many percentages have been tracked for averaging purposes */
};

struct s_pbox {
    enum e_box_status status;  
    double volume;  /* Current volume of parameter space */
    Range *p;  /* Current set of p ranges, p[nparms][2] */
    int vp_set; /* Indicating the set of data we are processing */
    struct s_yu **vp;  /* Current set of the overall hull of special p ranges for the entire time range 
                        * This is indexed by vp[nsets][DEdim] */
    struct s_inf **inf; /* Influence on F of each variable/equation pair 
                         * indexed by [DEdim+AEdim][nparms+1 (for unmeasured variables)] */
};

struct s_parms {
    int nboxes; /* number of boxes */
    int maxboxes; /* maximum number of boxes */
    double *p_gw;  /* array of p goal widths, p_gw[nparms] */
    double minred; /* mininum reduction fraction (of goal width) for success */
    double orig_vol; /* volume of original parameter space box */
    struct s_pbox *pbox;  /* array of parameter boxes, pbox[maxboxes] */
};

enum e_discrtype {parm_bc, var_bc, var_hc};

enum e_stencil {not_used=0, var=1, vec=2, var_vec=3 };
    /* Indicator of how a point in a discretization window is used in the discretization 
     * not_used (alpha and beta both zero)
     * var_only (alpha nonzero, beta zero)
     * vecfield_only (alpha zero, beta nonzero)
     * var_vecfield (alpha and beta nonzero) */

struct s_discr {
    enum e_discrtype discrtype; /* Discretization type. */
    int winsize;  /* size of (# of points in) discretization window, counting
                   both end points (so a one step formula would have size=2.*/
    enum e_stencil *stencil; /* Stencil for the discretization (arry of
                              size winsize) */
    double *alpha; /* coefficients of the variables, (array of size winsize) */
    double *beta;  /* coefficients of the vector field (array of size winsize) */
    int sub_window[2]; /* start and stop indices of the subwindow */
};

struct s_window { /* s_window contains information used by the program in various
                   * lower level functions - they are malloced and freed in acquire_input
                   * and write_output to avoid unnecessary time */
    int subbox_specs[3];
    double *split_ratios;
    int *time_split;
    int *eqn_order;
    int *windep[6];
    int windep_cor_flag;
    Range *p_prev_mono;
    int **windep_msrd[3];
    Range *partial;
    intRange *attempt;
    Range *p_hull;
    double *tuse;
    Range **y;
    Range **u;
    int *t_ind;
    Range *yext;
    Range *uext;
    intRange **yext_sort;
    intRange **uext_sort;
    Range **y_alt;
    Range *y_alt1;
    Range **u_alt;
    Range *p_alt;
    Range *z;
    double *p_ord0;
    double *p_ord1;
    double *p_ord2;
    Range *p_prev[3];
};

struct s_model {
    char *name;  /* model name */
    int DEdim;  /* number of dimensions for the differential equation */
    int OFdim;  /* number of dimensions for the external functions */
    int AEdim;  /* number of dimensions for the algebraic state equations */
    int TEdim;  /* number of total dimensions for the system (DEdim + AEdim) */
    int nparms; /* number of parameters in the differential equation */
    int tparms;  /* number of total parameters (including those of unmeasured variables) */
    Range *p_extreme; /* a priori bounds on the parameters */
    Range *y_extreme; /* a priori bounds on the variables */
    int n_ispecs; /* number of integer specifications */
    int *ispecs; /* model integer specifications */
    int n_fspecs; /* number of floating point specifications */
    double *fspecs; /* model floating point specifications */
    void (*determine_dependence)(int eqn, Range *p, Range t, Range *y, Range *u, Range partial,
                                struct s_model *model, int *dep, struct s_discr *discr, int parm);
    int (*determine_special_dependence)(int eqn, struct s_model *model, Range *p, 
                                        int var_vparms, double alpha, double beta, double h, 
                                        Range t, Range *y, Range *u, int* dep, int i);
    double (*vec_field)(int cor,int eqn, double *p, double a, double t, Range *y, Range *u, 
                        struct s_model *model); 
    double (*interval_vec_field)(int eqn, int cor, Range *p, double *t, Range **y, Range **u, 
                                 struct s_model *model, struct s_discr *discr);
    int (*invert_vec_field)(int eqn, int i, Range F, struct s_model *model, Range *p, double t, 
                            Range *y, Range *u);
    int **parmdep;  /* matrix indicating dependence of equation i of DE on 
                     * parameter j, parmdep[i][j], is (DEdim+AEdim) x nparms */
    int *parmdep_flag; /* array indicating if the given equation has any dependence on any parameter
                             * or unmeasured variable (summary over all reducable quantities) */
    int **vardep;  /* matrix indicating dependence of equation i of DE on 
                    * variable j, vardep[i][j], is (DEdim+AEdim) x DEdim */
    int **ofdep;  /* matrix indicating dependence of equation i of DE on 
                   * other function j, ofdep[i][j], is (DEdim+AEdim) x OFdim */
    int **varpartialdep;  /* matrix indicating dependence of partial derivatives of equation i of DE on 
                    * variable j, vardep[i][j], is (DEdim+AEdim) x DEdim */
    int **ofpartialdep;  /* matrix indicating dependence of partial derivatives of equation i of DE on 
                   * other function j, ofdep[i][j], is (DEdim+AEdim) x OFdim */
    double discr_error; /* used to make cuts at the appropriate place */
    Range partial; /* partial derivative for monotonicity check */
    struct s_window window;  /* struct containing storeage for use in the program */
    int var; /* temp variable for passing information to the invert function */
};

struct s_discr_set {
    char *name;  /* name of discretization set */
    int n_ispecs; /* number of integer specifications */
    int *ispecs; /* discretization integer specifications */
    int n_fspecs; /* number of floating point specifications */
    double *fspecs; /* discretization floating point specifications */
    int n_discr;   /* number of discretizations */
    intRange flag;  /* a flag to decide whether or not to read in specs, or accept them from default */
    struct s_discr *discr;  /* array of discretizations */
    int input_size;
};

struct s_discr_set_list {
    int n_discr_types; /* number of discretization types (i.e. A1OUT, B1OUT, etc.) */
    struct s_discr_set *discr_set; /* an array of discretization types */
    intRange sub_win_var; /* min/max subwindow lengths for unmeasured variable hull consistency */
    intRange sub_win_parm; /* min/max subwindow lengths for parameter hull consistency */
    double loop_specs[4]; /* user input: number of loops allowed without significant progress */
};
