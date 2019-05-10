#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "parm_red_structs.h"
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define ISZERO_INTERIOR(a) (((a)[0] * (a)[1]) < 0.0 ? 1 : 0 )

void int_vfield_sum(Range *x, Range result, struct s_discr *discr) {

    int i;
    result[0] = 0.0; result[1] = 0.0;
    
    for (i=discr->sub_window[0]; i<discr->sub_window[1]; i++) {
        result[0] += discr->beta[i]*x[i][0];
        result[1] += discr->beta[i]*x[i][1];
    }
    return;
}

int intersect(Range x, Range y) {
    /* intersect the the interval x with y, and update x */
	/* return value is -1 if intersection is empty, 1 if x is contained in y, 2 if x is updated. */
    int ok = 1;
    
    if (x[0] > y[1] || x[1] < y[0]) {
         ok = -1;
    }
    else {
        if (x[0] < y[0]) {
            x[0] = y[0];
            ok = 2;
        }
		if (x[1] > y[1]) {
            x[1] = y[1];
            ok = 2;
        }
    }
    return ok;
}

void int_multiply(Range x, Range y, Range result) {
    /* Interval multiplication function. */
        
    double temp[4];
    temp[0] = x[0]*y[0];
    temp[1] = x[1]*y[1];
    
    if (x[0] >= 0 && y[0] >= 0) {
        result[0] = temp[0];
        result[1] = temp[1];
    }
    else {
        temp[2] = x[1]*y[0];
        temp[3] = x[0]*y[1];
        
        if (temp[0] < temp[1]) {
            result[0] = temp[0];
            result[1] = temp[1];
        }
        else {
            result[0] = temp[1];
            result[1] = temp[0];
        }
        if (temp[2] > result[1]) result[1] = temp[2];
        else if (temp[2] < result[0]) result[0] = temp[2];
        
        if (temp[3] > result[1]) result[1] = temp[3];
        else if (temp[4] < result[0]) result[0] = temp[4];
    }
    return;
}

void int_divide(Range x, Range y, Range result) {
    /* Interval divide function. */
    if(y[0] <= 0.0 && y[1] >= 0.0){
        result[0] = -DBL_MAX;
        result[1] = DBL_MAX;
    }
    else{
        double temp1,temp2,temp3 = x[1]/y[0],temp4 = x[0]/y[1];
        
        if ((temp1=(x[0]/y[0])) < (temp2=(x[1]/y[1]))) {
            result[0] = temp1;
            result[1] = temp2;
        }
        else {
            result[0] = temp2;
            result[1] = temp1;
        }
        
        if (temp3 > result[1]) result[1] = temp3;
        else if (temp3 < result[0]) result[0] = temp3;
        
        if (temp4 > result[1]) result[1] = temp4;
        else if (temp4 < result[0]) result[0] = temp4;
    }
    return;
}

void int_subtract(Range x, Range y, Range result) {
    /* Interval subtract function. */
    result[0] = x[0] - y[1];
    result[1] = x[1] - y[0];
    return;
}

void int_add(Range x, Range y, Range result) {
    /* Interval add function. */
    result[0] = x[0] + y[0];
    result[1] = x[1] + y[1];
    return;
}

void multiply(double x, Range result) {
    /* Interval multiply by scalar function. */
    if (x >= 0) {
        result[0] = x*result[0];
        result[1] = x*result[1];
    }
    else{
        double temp = result[0];
        result[0] = x*result[1];
        result[1] = x*temp;
    }
    return;
}


void int_square(Range x, Range result){
    /* Interval square function. */
	result[0] = x[0]*x[0];
	result[1] = x[1]*x[1];
	if (result[0] > result[1]) {
		double temp = result[0];
		result[0] = result[1];
		result[1] = temp;
	}
    if (ISZERO_INTERIOR(x)) result[0] = 0.0;
    return;
}

void int_pow(Range x, Range y, Range result){
    
    double a[4];
    int i;
    
    if(x[0] < 0.0 || (x[0] == 0.0 && y[0] <= 0.0)) {
        printf("ERROR in int_pow: base must be strictly greater than 0 or,\n   if x is zero, y must be strictly greater than zero.\n"); exit(1);
    }
    a[0] = pow(x[0],y[0]);
    a[1] = pow(x[0],y[1]);
    a[2] = pow(x[1],y[0]);
    a[3] = pow(x[1],y[1]);
    
    result[0] = a[0];
    result[1] = a[0];
    for (i=1; i<4; i++) {
        if (a[i] < result[0]) result[0] = a[i];
		else if (a[i] > result[1]) result[1] = a[i];
    }
    return;
}

void int_log(Range x, Range result) {
    /* Interval natural logarithm function. */
    
	for (int i=0; i<2; i++) {
		if (x[i] <= 0.0) result[i] = -DBL_MAX;
		else  result[i] = log(x[i]);
	}
    return;
}

void int_log10(Range x, Range result) {
    /* Interval natural logarithm base 10 function. */
    
	for (int i=0; i<2; i++) {
		if (x[i] <= 0.0) result[i] = -DBL_MAX;
		else  result[i] = log10(x[i]);
	}
    return;
}


void int_put(Range x, Range result){
    
    result[0] = x[0];
    result[1] = x[1];
    return;
}

void int_cosine(Range x, Range result) {
    /* Interval cosine function. */
    
    static double pi = 3.14159265358979323846;
    Range z;
    double temp;
    
    /* Scale x by 1/pi and shift so that the low end is in range [0,2) */
    z[0] = fmod(x[0]/pi,2.0);
    if (z[0] < 0.0) z[0] += 2.0;
    /* Shift x[1]/pi the same amount. */
    z[1] = (x[1] - x[0])/pi + z[0];  /* x[1]/pi - (x[0]/pi - z[0]) */
    
    if (z[1] >= 3.0) {
        result[0] = -1.0;
        result[1] = 1.0;
    }
    else if (z[1] >= 2.0) {
        result[1] = 1.0;
        if (z[0] <= 1.0) 
            result[0] = -1.0;
        else {
            result[0] = cos(x[0]);
            temp = cos(x[1]);
            if (temp < result[0]) result[0] = temp;
        }
    }
    else if (z[1] >= 1.0) {
        if (z[0] <= 1.0) {
            result[0] = -1.0;
            result[1] = cos(x[0]);
            temp = cos(x[1]);
            if (temp > result[1]) result[1] = temp;
        }
        else {
            result[0] = cos(x[0]);
            result[1] = cos(x[1]);
        }
    }
    else {
        result[0] = cos(x[1]);
        result[1] = cos(x[0]);
    }
    return;
}

void int_sine(Range x, Range result) {
    /* Interval sine function. */
    
    static double pi = 3.14159265358979323846;
    Range z;
    double temp;
    
    /* Scale x by 1/pi and shift so that the low end is in range [0,2) */
    z[0] = fmod(x[0]/pi,2.0);
    if (z[0] < 0.0) z[0] += 2.0;
    /* Shift x[1]/pi the same amount. */
    z[1] = (x[1] - x[0])/pi + z[0];  /* x[1]/pi - (x[0]/pi - z[0]) */
    if (z[1] >= 3.5) {
        result[0] = -1.0;
        result[1] = 1.0;
    }
    else if (z[1] >= 2.5) {
        result[1] = 1.0;
        if (z[0] <= 1.5) 
            result[0] = -1.0;
        else {
            result[0] = sin(x[0]);
            temp = sin(x[1]);
            if (temp < result[0]) result[0] = temp;
        }
    }
    else if (z[1] >= 1.5) {
        if (z[0] <= 0.5) {
            result[0] = -1.0;
            result[1] = 1.0;
        }
        else if (z[0] <= 1.5) {
            result[0] = -1.0;
            result[1] = sin(x[0]);
            temp = sin(x[1]);
            if (temp > result[1]) result[1] = temp;
        }
        else {
            result[0] = sin(x[0]);
            result[1] = sin(x[1]);
        }
    }
    else if (z[0] <= 0.5) {
        if (z[1] >= 0.5) {
            result[1] = 1.0;
            result[0] = sin(x[0]);
            temp = sin(x[1]);
            if (temp < result[0]) result[0] = temp;
        }
        else {
            result[0] = sin(x[0]);
            result[1] = sin(x[1]);
        }
    }
    else {
        result[0] = sin(x[1]);
        result[1] = sin(x[0]);
    }
    return;
}
