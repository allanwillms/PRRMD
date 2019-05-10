#ifndef TypedefRange
#define TypedefRange
typedef double Range[2];
typedef int intRange[2];
typedef int intRange3[3];
#endif

void int_sine(Range x, Range result);

void int_cosine(Range x, Range result);

void int_vfield_sum(Range *x, Range result, struct s_discr *discr);

int intersect(Range x, Range y);

void multiply(double x, Range result);

void int_divide(Range x, Range y, Range result);

void int_multiply(Range x, Range y, Range result);

void int_subtract(Range x, Range y, Range result);

void int_add(Range x, Range y, Range result);

void int_square(Range x, Range result);

void int_pow(Range x, Range y, Range result);

void int_log(Range x, Range result);

void int_log10(Range x, Range result);

void int_put(Range x, Range result);