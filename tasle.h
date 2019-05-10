/* May 25/12 - current for non-constant band */

#ifndef TypedefRange
#define TypedefRange
typedef double Range[2];
typedef int intRange[2];
typedef int intRange3[3];
#endif

struct s_band {
    intRange *pivots;
    int *interval_list;
    Range *tjoinpoints;
    Range *yjoinpoints;
    intRange *status;
};

int set_status(intRange *status, int interval, int step, int status_number, int info);

int slope_compare(int a, int b, int c, int d, double *t, double *y, Range **store, intRange **check, int flag);

double slope(int a, int b, double *t, double *y, Range **store, intRange **check, int flag);

int distance_check(int a, int b, double *t, Range **store, intRange **check, double user_tmin);

int distance_compare(int a, int b, int c, int d, double *t, Range **store, intRange **check);

double distance(int a, int b, double *t, Range **store, intRange **check);

int check_validity(double *t, double *y, intRange *pivots, Range **store, intRange **check, int p, double user_tmin, int size, int top);

int convex(double *t, double *y, int start, int size, intRange *pivots, int initial, int top, int insert, Range **store, intRange **check);

int tasle(int N, double *t, double *y, double user_tmin, double min_height, double **Tout, Range **Yout);

int process_status(double *t, double *y, intRange *status, int size, double user_tmin, intRange *pivots, Range **store, intRange **check, int *interval_list, intRange *encroach_pivot, int **potential_pivot, int p, int k, int iteration, int top, int counter);

int process_interval(double *t, double *y, intRange *status, int size, double user_tmin, intRange *pivots, int interval_start, int interval_finish, int interval_next, int interval_before, int p, int k, int iteration, Range **store, intRange **check, int *interval_list, intRange *encroach_pivot, int **potential_pivot, int top);

int find_next_pivot(double *t, double *y, intRange *pivots, int p, double user_tmin, int interval_start, int interval_finish, int interval_next, int interval_before, Range **store, intRange **check, intRange *encroach_pivot, int **potential_pivot, int flag_left, int flag_right, int e, int top);
    
int phase_three(double *t, double *y, intRange *pivots, int *interval_list, int p, Range **store, intRange **check, double user_tmin, intRange *status, int iteration, int top);

int find_inflection_points(double *t, double *y, double user_tmin, int *interval_list, Range **store, intRange **check, int top, int k, int up_down, int *left_encroach, int *right_encroach, int flag);

int concave(double *t, double *y, int start, int finish, int before, int next, int start_alt, int finish_alt, intRange *pivots, Range **store, intRange **check, int *alternate, int alternate_counter, double input_area, double user_tmin, int initial, int top, int insert);