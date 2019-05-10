void aquire_input(struct s_model *model, struct s_discr_set_list *discr_set_list,
                  struct s_data *data, struct s_parms *parms, char *fname, char *argv[]);

int parse_colspec(char *expr, int *file, double *val);

double parse_gw(char *expr, double mmm[3]);

void write_output(struct s_model *model, struct s_discr_set_list *discr_set_list,
                  struct s_data *data, struct s_parms *parms, char *fname);

void reduce_all(struct s_model *model, struct s_discr_set_list *discr_set_list, 
                struct s_data *data, struct s_parms *parms);

void split_box(struct s_model *model, struct s_parms *parms, int boxnum, 
              struct s_data *data, int split_direction, double slice);

void reduce_box(struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_data *data,
                struct s_pbox *pbox, double y_minred, double *p_gw, double minred);

int parm_goals_met(struct s_pbox *pbox, int nparms, double *p_gw);

int process_set(struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_set *set,
                struct s_pbox *pbox, double y_minred, double *p_gw, double minred);

int process_window(struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_set *set,
                   struct s_pbox *pbox, int eqn, double y_minred, double *p_gw,
                   double minred, double winstart, double h);

int process_subboxes(struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_set *set,
                     struct s_pbox *pbox, int eqn, double y_minred, double *t, Range **y, Range **u, 
                     double h, int mono, int mono2, intRange *attempt, int split_option, 
                     int p_index);

int get_range(double t_req, struct s_yu *v, int t_ind, Range val);

int guess_binsearch(double req, int n, double *v, int *guess);

void check_monotonicity(struct s_model *model, struct s_discr *discr, struct s_unmsrd *unmsrd,
                        struct s_pbox *pbox, int eqn, double *t, Range **y, Range **u, double h, 
                        int base_winsize, int *windep[3], int *mono, int *mono2);

double select_window(struct s_set *set, double winend_prev, int prev_progress, 
                      int winsize, double *winstart, double *h, double overlap);

int process_equation(struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_unmsrd *unmsrd,
                     struct s_pbox *pbox, int eqn, double y_minred, double *t, Range **y, Range **u, 
                     double h, int mono, int mono2, intRange *attempt);

int process_algebraic_equation(struct s_model *model, struct s_discr_set_list *discr_set_list,
                               struct s_unmsrd *unmsrd, struct s_pbox *pbox, int eqn, double y_minred,  
                               double *t, Range **y, Range **u, int mono, int mono2, 
                               intRange *attempt);

void get_put_ordered_parms(int flag, int tparms, struct s_pbox *pbox, int *windep, double *p_ord[2]);

int hull_consistency (struct s_model *model, struct s_discr_set_list *discr_set_list, struct s_unmsrd *unmsrd,
                      struct s_pbox *pbox, int eqn, double *t, Range **y, Range **u,
                      double h, intRange *attempt, double *p_ord[2], int flag);

int parm_box_consistency(struct s_model *model, struct s_discr *discr, struct s_unmsrd *unmsrd,
                         struct s_pbox *pbox, int cor, double *p_cor, double *p_opp, int eqn,  
                         double *t, Range **y, Range **u, int *windep, double h, intRange *attempt, int mono);

double constraint_eval(struct s_model *model, struct s_discr *discr, struct s_unmsrd *unmsrd, int cor, 
                       int um_parm, double *p, double *p_opp, int eqn, double *t, Range **y, Range **u, 
                       double h, int mono, int *windep);

void int_sine(Range x, Range result);

void int_cosine(Range x, Range result);

int intersect(Range x, Range y);

void int_vfield_sum(Range *x, Range result, struct s_discr *discr);

void multiply(double x, Range result);

void int_divide(Range x, Range y, Range result);

void int_multiply(Range x, Range y, Range result);

void int_subtract(Range x, Range y, Range result);

void int_dependent_subtract(Range x, Range y, Range result);

void int_add(Range x, Range y, Range result);

void int_pow(Range x, Range y, Range result);

void int_log(Range x, Range result);

void int_put(Range x, Range result);
