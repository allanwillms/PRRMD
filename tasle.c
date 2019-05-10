#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "tasle.h"
#include <time.h>

// valgrind --leak-check=yes --dsymutil=yes ./tasle-shell 32
// valgrind --tool=cachegrind --dsymutil=yes ./tasle-shell 50
// cg_annotate cachegrind.out.90871 --auto=yes

//  ### means something has to be changed still
//  ### going with alternate definition of concave hull - checking right slope error (TBD)

//#define db1 /* db1 includes interval and procedural information at every iteration */
//#define db2 
//#define db3
//#define db4  /* db4 includes pivots */
//#define db5 /* db5 includes all phase 3 information */
//#define db6 /* db6 includes flop count information */
//#define db7 /* db7 includes timing information */
//#define db8a /* db8 includes band expansions at the ends */
#define db8b
//#define db9 /* db9 includes testing information, check validity and slope call index errors */ 

#define abs_value(x) (x >= 0 ? x : -x) 
#define MIN(a,b) (a < b ? a : b)
#define MIN3(a,b,c) (MIN(a, MIN(b,c)))
#define YSPAN(a) ((a)[1] - (a)[0])

/* flop_counter used for (+,-,*,/) operations only, flop_counter_2 used for (double comparisons) */
#ifdef db6
    long flop_counter = 0, flop_counter_2 = 0;
#endif

/* setting the status of a given interval */
int set_status(intRange *status, int interval, int step, int status_number, int info){
    status[2*interval][step%2] = status_number;
    status[2*interval+1][step%2] = info;
    return 1;
}

int slope_compare(int a, int b, int c, int d, double *t, double *y, Range **store, intRange **check, int flag){
    
#ifdef db6
    if(check[a][b-a-1][0] == 0)  flop_counter += 3;
    if(check[c][d-c-1][0] == 0)  flop_counter += 3; 
    flop_counter_2 ++;
#endif
#ifdef db9
    if(a == b || a > b || a < 0 || b < 0){ printf("slope_compare a = %d,b = %d\n",a,b); exit(1); }
    if(c == d || c > d || c < 0 || d < 0){ printf("slope_compare c = %d,d = %d\n",c,d); exit(1); }
#endif
    if(check[a][b-a-1][0] == 0){
        store[a][b-a-1][0] = flag == 0 ? (y[b]-y[a])/(t[b]-t[a]) : (y[a]-y[b])/(t[b]-t[a]);
        check[a][b-a-1][0] = 1;
    }
    if(check[c][d-c-1][0] == 0){
        store[c][d-c-1][0] = flag == 0 ? (y[d]-y[c])/(t[d]-t[c]) : (y[c]-y[d])/(t[d]-t[c]);
        check[c][d-c-1][0] = 1;
    }
    
    return flag == 0 ? (store[a][b-a-1][0] <= store[c][d-c-1][0]) : (store[a][b-a-1][0] >= store[c][d-c-1][0]);
}

double slope(int a, int b, double *t, double *y, Range **store, intRange **check, int flag){
    
#ifdef db6
    if(check[a][b-a-1][0] == 0)  flop_counter += 3;
#endif
#ifdef db9
    if(a == b || a > b || a < 0 || b < 0){ printf("slope a = %d,b = %d\n",a,b); exit(1); }
#endif
    if(check[a][b-a-1][0] == 0){
        store[a][b-a-1][0] = flag == 0 ? (y[b]-y[a])/(t[b]-t[a]) : (y[a]-y[b])/(t[b]-t[a]);
        check[a][b-a-1][0] = 1;
    }
    return store[a][b-a-1][0];
}

/* calculates the horizontal distance between two data points */
int distance_check(int a, int b, double *t, Range **store, intRange **check, double comparison){
    
#ifdef db6
    if(check[a][b-a-1][0] == 0)  flop_counter++;
    flop_counter_2 ++;
#endif
#ifdef db9
    if(a == b || a > b || a < 0 || b < 0){ printf("distance_check a = %d,b = %d\n",a,b); exit(1); }
#endif
    if(check[a][b-a-1][1] == 0){
        store[a][b-a-1][1] = t[b]-t[a];
        check[a][b-a-1][1] = 1;
    }
    
    return comparison <= store[a][b-a-1][1];
}

int distance_compare(int a, int b, int c, int d, double *t, Range **store, intRange **check){
    
#ifdef db6
    if(check[a][b-a-1][0] == 0)  flop_counter++;
    if(check[c][d-c-1][0] == 0)  flop_counter++;
    flop_counter_2 ++;
#endif
#ifdef db9
    if(a == b || a > b || a < 0 || b < 0){ printf("distance_compare a = %d,b = %d\n",a,b); exit(1); }
    if(c == d || c > d || c < 0 || d < 0){ printf("distance_compare c = %d,d = %d\n",c,d); exit(1); }
#endif
    if(check[a][b-a-1][1] == 0){
        store[a][b-a-1][1] = t[b]-t[a];
        check[a][b-a-1][1] = 1;
    }
    if(check[c][d-c-1][1] == 0){
        store[c][d-c-1][1] = t[d]-t[c];
        check[c][d-c-1][1] = 1;
    }

    return store[a][b-a-1][1] < store[c][d-c-1][1];
}

double distance(int a, int b, double *t, Range **store, intRange **check){
    
#ifdef db6
    if(check[a][b-a-1][1] == 0)  flop_counter++;
#endif
#ifdef db9
    if(a == b || a > b || a < 0 || b < 0){ printf("distance a = %d,b = %d\n",a,b); exit(1); }
#endif
    
    if(check[a][b-a-1][1] == 0){
        store[a][b-a-1][1] = t[b]-t[a];
        check[a][b-a-1][1] = 1;
    }
    
    return store[a][b-a-1][1];
}

int check_validity(double *t, double *y, intRange *pivots, Range **store, intRange **check, int p, double user_tmin, int size, int top){
    
    double slope1, slope2, slope3;
    int q;
    
    for(q = 2; q <= p-2; q++){
        if(t[pivots[q][top]]-t[pivots[q-1][top]] < user_tmin){
            if(pivots[q][top] != pivots[q-1][top] && pivots[q+1][top] != pivots[q][top] && pivots[q-1][top] != pivots[q-2][top]){
                slope1 = (y[pivots[q-2][top]]-y[pivots[q-1][top]])/(t[pivots[q-2][top]]-t[pivots[q-1][top]]);
                slope2 = (y[pivots[q-1][top]]-y[pivots[q][top]])/(t[pivots[q-1][top]]-t[pivots[q][top]]);
                slope3 = (y[pivots[q][top]]-y[pivots[q+1][top]])/(t[pivots[q][top]]-t[pivots[q+1][top]]);
                
                if(!((slope3 >= slope2 && slope2 >= slope1) || (slope1 >= slope2 && slope2 >= slope3))){
                    printf("FAILURE q = %d\n",q);
                    printf("   %f %f\n",t[pivots[q-1][top]],t[pivots[q][top]]);
                    printf("   (%d) %d %d (%d)\n",pivots[q-2][top],pivots[q-1][top],pivots[q][top],pivots[q+1][top]);
                    printf("   distance from (%d,%d) = %f\n",pivots[q-1][top],pivots[q][top],t[pivots[q][top]]-t[pivots[q-1][top]]);
                    printf("   slopes = (%d,%d) %f, (%d,%d) %f, (%d,%d) %f\n",pivots[q-2][top],pivots[q-1][top],slope1,pivots[q-1][top],pivots[q][top],slope2,pivots[q][top],pivots[q+1][top],slope3);
                    exit(1);
                }
            }
            
        }
    }  
//    int pivot_shift = 1;
//    for(q = 2; q <= size-1; q++){
//        if(t[q] > t[pivots[pivot_shift+1][top]]){
//            pivot_shift++;
//        }
//        if(q != pivots[pivot_shift+1][top]){
//            if(pivots[pivot_shift][top] != pivots[pivot_shift+1][top]){
//                if(0.975*y[q] > y[pivots[pivot_shift][top]] + (y[pivots[pivot_shift][top]]-y[pivots[pivot_shift+1][top]])/(t[pivots[pivot_shift][top]]-t[pivots[pivot_shift+1][top]])*(t[q]-t[pivots[pivot_shift][top]]) ){
//                    printf("FAILURE\n");
//                    printf("   %d %f  (y-value check %f <= %f)\n",q,t[q],y[q],y[pivots[pivot_shift][top]] + (y[pivots[pivot_shift][top]]-y[pivots[pivot_shift+1][top]])/(t[pivots[pivot_shift][top]]-t[pivots[pivot_shift+1][top]])*(t[q]-t[pivots[pivot_shift][top]]));
//                    exit(1);
//                }
//            }
//        }
//    }
    
    return 1;
}

int convex(double *t, double *y, int start, int size, intRange *pivots, int initial, int top, int insert, Range **store, intRange **check){
    
    int *hull; hull = (int*) malloc((size-start+1)*sizeof(int));
    int i, q, k = 0;
    
    /* add/remove points from the hull as appropriate */
    for (i = start; i < size; ++i) {
        while (k >= 2 && slope_compare(hull[k-2],i,hull[k-2],hull[k-1],t,y,store,check,top)) --k;
        hull[k++] = i;
    }    
    for (i = size-2, q = k+1; i >= start; --i) {
        while (k >= q && slope_compare(i,hull[k-2],hull[k-1],hull[k-2],t,y,store,check,top)) --k;
        hull[k++] = i;
    }
    
    /* storing points in the list of pivots as appropriate - all bookkeeping below */
    if(initial == 1){
        for(q = 0; q < k; q++){
            pivots[q][1] = hull[q];
            if(hull[q] == size-1){   pivots[q+1][1] = 0; q = k;  }
        }
        for(q = k-1; q > 0; q--){
            pivots[k-1-q][0] = hull[q];
            if(hull[q] == size-1){   pivots[k-q][0] = 0; q = 0;  }
        }
    }
    else if(initial == 0){
        for(q = k-1; q > 0; q--){
            if(hull[q] == size-1){   i = k-1-q-1; q = 0;  }
            else   pivots[insert+k-1-q][top] = hull[q];
        }
    }
    else{
        for(q = k-1; q > 0; q--){
            pivots[insert+k-1-q][top] = hull[q];
            if(hull[q] == size-1){   i = k-1-q; q = 0;  }
        }
    }    
    free(hull);
    return i;
}

/* tasle function starts here - input/output described above */
int tasle(int size, double *t, double *y, double user_tmin, double min_height, double **Tout, Range **Yout){
    
	int i, q, p, k, interval_counter, iteration, completion_counter, completion, top, p_upper, p_lower; 
	double mindt = user_tmin+1, final_height = 0;
    
    struct s_band sol;
    sol.pivots = (intRange*) malloc((size+1)*sizeof(intRange)); 
	sol.interval_list = (int*) malloc((size+1)*sizeof(int));
    sol.status = (intRange*) calloc(2*size,sizeof(intRange));   

    Range **store; store = (Range**) malloc((size)*sizeof(Range*));
    for(i = 0; i < size; i++)    store[i] = (Range*) malloc((size-i)*sizeof(Range));

    intRange **check; check = (intRange**) calloc((size),sizeof(intRange*));
    for(i = 0; i < size; i++)    check[i] = (intRange*) calloc((size-i),sizeof(intRange));
    
    intRange *encroach_pivot; encroach_pivot = (intRange*) calloc(size,sizeof(intRange));

    int **potential_pivot; potential_pivot = (int**) malloc(size*sizeof(int*));
    for(i = 0; i < size; i++)    potential_pivot[i] = (int*) malloc(5*sizeof(int));
    
	/* START - CONSIDER THE SIMPLEST CASE FIRST - we start by considering the simplest case - find the minimum distance between two t-values - if this is less than the user specified minimum, then the best solution is simply to connect all data points - final sol will have a height of zero */
	for(i = 0; i <= size-2; i++){
        if(distance_check(i,i+1,t,store,check,mindt) == 0)   mindt = distance(i,i+1,t,store,check);
    }
            
	if(user_tmin <= mindt){	
		*Tout = (double*) malloc((size)*sizeof(double));
		*Yout = (Range*) malloc((size)*sizeof(Range));
		for(i = 0; i <= size-1; i++){
			(*Tout)[i] = t[i]; //(*Tout)[i][1] = t[i];
			(*Yout)[i][0] = y[i] - 0.5*min_height; 
            (*Yout)[i][1] = y[i] + 0.5*min_height;
		}        
        printf("%f\n",final_height);
		return size; 
	}
    
    /* TIMING begins */
#ifdef db7
    clock_t start, end;
    double cpu_time_used;
    start = clock();
#endif
    
    /* convex hull to initialize - second and third-last arguments dummy arguments for this usage */
    convex(t, y, 0, size, sol.pivots, 1, 0, 0, store, check);
    
    /* for the top band and then the reflected lower band */
    for (top = 0; top <= 1; top++) {
        /* if we are banding the bottom, then flip all the data */
        if(top == 1){   for(i = 0; i <= size-1; i++)   y[i] *= -1;  }
        
        /* bookkeeping */
        q = 1;  while(sol.pivots[q][top] != 0)  q++;   p = q; 
        if(top == 1){   for(i = 0; i <= 2*p; i++)  sol.status[i][0] = 0;  }
        iteration = 0; completion = 0;
        
        /* now we loop until we can make no further progress in each interval */
        while(completion == 0){
            
            iteration++; completion_counter = 0;
#ifdef db1
            printf("Iteration %d\n",iteration);
#endif
            /* place all pivots from the last iteration into interval_list for the next interation */
            for(q = 0; q <= p-1; q++)   sol.interval_list[q] = sol.pivots[q][top];
            interval_counter = p-1; p = 0;
            
            /* for all intervals in the interval_list from the previous iteration */
            for(k=0; k <= interval_counter-1; k++){
                
                /* just a check to see if we are done completely and can exit the largest loop */
                completion_counter += sol.status[2*k][(iteration-1)%2] == 1 ? 1 : 0;
                /* place the start of the interval in the pivot list */
                sol.pivots[p][top] = sol.interval_list[k]; p++;
                /* process each interval differently depending on the current status value */
                p = process_status(t, y, sol.status, size, user_tmin, sol.pivots, store, check, sol.interval_list,  encroach_pivot, potential_pivot, p, k, iteration, top, interval_counter-1);
            }
            /* put the last interval end point in the list of pivots */
            sol.pivots[p][top] = sol.interval_list[k]; p++;
            
            /* if we cannot make any more progress, then set the flag to 1 so that we can exit this loop */
            completion = completion_counter == interval_counter ? 1 : 0;
        }
   
#ifdef db4
        printf("PIVOTS:\n");
        for(q = 0; q <= p; q++)  printf("%d ",sol.pivots[q][top]);
        printf("\n");
#endif
        /* now begin the smoothing procedure around convex areas of the band */
        p = phase_three(t, y, sol.pivots, sol.interval_list, p, store, check, user_tmin, sol.status, iteration, top);
        
#ifdef db9
        /* perform a quick test here to check to see if the final answer is correct - remove this for the final version */
        check_validity(t, y, sol.pivots, store, check, p, user_tmin, size, top);
#endif
        if(top == 0)  p_upper = p;
        else  p_lower = p;
    }
    
    for (i = 0; i <= size-1; i++)  y[i] *= -1; /* flip data back for final consideration */
    
    /* calculation of final list of pivots - based off entire list of upper and lower pivots */
    sol.tjoinpoints = (Range*) malloc((p_upper+p_lower+1)*sizeof(Range));
    sol.yjoinpoints = (Range*) malloc((p_upper+p_lower+1)*sizeof(Range));
    sol.tjoinpoints[0][0] = t[0];
    sol.yjoinpoints[0][0] = y[0];
    sol.yjoinpoints[0][1] = y[0];
    p = sol.pivots[0][1] == sol.pivots[1][1] ? 2 : 1;
    q = sol.pivots[0][0] == sol.pivots[1][0] ? 2 : 1;
    i = 1;
    while (p < p_lower && q < p_upper) {
        if (sol.pivots[p][1] == sol.pivots[q][0]) {  
            sol.tjoinpoints[i][0] = t[sol.pivots[q][0]];
            sol.yjoinpoints[i][0] = y[sol.pivots[q][0]];
            sol.yjoinpoints[i][1] = y[sol.pivots[q][0]];
            p++; q++; 
        }
        else if (sol.pivots[p][1] < sol.pivots[q][0]) {
            sol.tjoinpoints[i][0] = t[sol.pivots[p][1]];
            sol.yjoinpoints[i][0] = y[sol.pivots[p][1]];
#ifdef db6
            flop_counter += 2; 
#endif
            sol.yjoinpoints[i][1] = y[sol.pivots[q-1][0]] + slope(sol.pivots[q-1][0],sol.pivots[q][0],t,y,store,check,0)*distance(sol.pivots[q-1][0],sol.pivots[p][1],t,store,check);
            p++;
        }
        else{
            sol.tjoinpoints[i][0] = t[sol.pivots[q][0]];
#ifdef db6
            flop_counter += 2;
#endif
            sol.yjoinpoints[i][0]= y[sol.pivots[p-1][1]] + slope(sol.pivots[p-1][1],sol.pivots[p][1],t,y,store,check,0)*distance(sol.pivots[p-1][1],sol.pivots[q][0],t,store,check);
            sol.yjoinpoints[i][1] = y[sol.pivots[q][0]];
            q++;
        }
#ifdef db6
        if(i > 1) {  flop_counter += 12; flop_counter_2 += 2; }
#endif
        /* eliminating any collinear sets of three points */
        if((i <= 1 || (sol.yjoinpoints[i][0]-sol.yjoinpoints[i-1][0])/(sol.tjoinpoints[i][0]-sol.tjoinpoints[i-1][0]) != (sol.yjoinpoints[i-1][0]-sol.yjoinpoints[i-2][0])/(sol.tjoinpoints[i-1][0]-sol.tjoinpoints[i-2][0]) || (sol.yjoinpoints[i][1]-sol.yjoinpoints[i-1][1])/(sol.tjoinpoints[i][0]-sol.tjoinpoints[i-1][0]) != (sol.yjoinpoints[i-1][1]-sol.yjoinpoints[i-2][1])/(sol.tjoinpoints[i-1][0]-sol.tjoinpoints[i-2][0]))){
            i++;
        }
        else{
            sol.tjoinpoints[i-1][0] = sol.tjoinpoints[i][0];
            sol.yjoinpoints[i-1][0] = sol.yjoinpoints[i][0];
            sol.yjoinpoints[i-1][1] = sol.yjoinpoints[i][1];
        }
        
    }
    sol.tjoinpoints[i-1][0] = t[size-1];
    
#ifdef db8a
    /* expand bands at the start of the data set - first is the first non-zero vertical distance between bands considering pivot points from left to right - second is the maximum amount we can expand without violating the slope condition on the second segment on the top and then bottom - finally, we add the minimum equally to the top and bottom */
    k = 1; double a[3]; intRange b;
    while((a[0] = 0.5*(sol.yjoinpoints[k][1]-sol.yjoinpoints[k][0])) == 0)  k++;
    
    for(k = 0; k <= 1; k++){
        if(!distance_check(sol.pivots[1][k],sol.pivots[2][k],t,store,check,user_tmin) && slope_compare(sol.pivots[0][k],sol.pivots[1][k],sol.pivots[1][k],sol.pivots[2][k],t,y,store,check,top) == 1-k && slope_compare(sol.pivots[1][k],sol.pivots[2][k],sol.pivots[2][k],sol.pivots[3][k],t,y,store,check,top) == 1-k){
            a[k+1] = y[sol.pivots[1][k]] - slope(sol.pivots[1][k],sol.pivots[2][k],t,y,store,check,top)*(distance(sol.pivots[0][k],sol.pivots[1][k],t,store,check)) - y[sol.pivots[0][k]];
#ifdef db6
            flop_counter += 3;
#endif;
            a[k+1] = abs_value(a[k+1]);
        }
        else a[k+1] = a[0];
    }
    
    
#ifdef db6
    flop_counter += 2;
    flop_counter_2 += 2;
#endif
    double min3expand = MIN3(a[0],a[1],a[2]);
    min3expand = 1.5*a[0];
    sol.yjoinpoints[0][0] -= min3expand;
    sol.yjoinpoints[0][1] += min3expand;

    /* go along the top */
    p=0;
    while(sol.pivots[p][0] == 0) p++;
    k = 1;
    while(sol.tjoinpoints[k][0] < t[sol.pivots[p][0]]){
        sol.yjoinpoints[k][1] = y[sol.pivots[p][0]] + (y[sol.pivots[p][0]]-sol.yjoinpoints[0][1])/(t[sol.pivots[p][0]]-sol.tjoinpoints[0][0])*(sol.tjoinpoints[k][0]-t[sol.pivots[p][0]]);
        k++;
#ifdef db6
        flop_counter += 5;
#endif
    }
    /* go along the top */
    p = 0;
    while(sol.pivots[p][1] == 0) p++;
    k = 1;    
    while(sol.tjoinpoints[k][0] < t[sol.pivots[p][1]]){
        sol.yjoinpoints[k][0] = y[sol.pivots[p][1]] + (y[sol.pivots[p][1]]-sol.yjoinpoints[0][0])/(t[sol.pivots[p][1]]-sol.tjoinpoints[0][0])*(sol.tjoinpoints[k][0]-t[sol.pivots[p][1]]);
        k++;
#ifdef db6
        flop_counter += 5;
#endif
    }
    
    /* expand bands at the right hand side of the data set */
    k = i-2; b[0] = p_upper; b[1] = p_lower;
    while((a[0] = 0.5*(sol.yjoinpoints[k][1]-sol.yjoinpoints[k][0])) == 0)  k--;
    
    for(k = 0; k <= 1; k++){
        if(!distance_check(sol.pivots[b[k]-3][k],sol.pivots[b[k]-2][k],t,store,check,user_tmin) && slope_compare(sol.pivots[b[k]-4][k],sol.pivots[b[k]-3][k],sol.pivots[b[k]-3][k],sol.pivots[b[k]-2][k],t,y,store,check,top) == 1-k && slope_compare(sol.pivots[b[k]-3][k],sol.pivots[b[k]-2][k],sol.pivots[b[k]-2][k],sol.pivots[b[k]-1][k],t,y,store,check,top) == 1-k){
            a[k+1] = (y[sol.pivots[b[k]-2][k]] - slope(sol.pivots[b[k]-3][k],sol.pivots[b[k]-2][k],t,y,store,check,top)*distance(sol.pivots[b[k]-2][k],sol.pivots[b[k]-1][k],t,store,check) - y[sol.pivots[b[k]-1][k]]);
            a[k+1] = abs_value(a[k+1]);    
#ifdef db6
            flop_counter += 3;
#endif;
        }
        else a[k+1] = a[0];
    }
    
#ifdef db6
    flop_counter += 2;
    flop_counter_2 += 2;
#endif
    min3expand = MIN3(a[0],a[1],a[2]);
    min3expand = a[0];
    sol.yjoinpoints[i-1][0] -= min3expand;
    sol.yjoinpoints[i-1][1] += min3expand;
    sol.tjoinpoints[i-1][0] = t[size-1];
    //sol.tjoinpoints[i-1][1] = t[size-1];
    
    if(sol.pivots[p_upper-2][0] < sol.pivots[p_lower-2][1]){ /* if last pivot on top occurs before bottom, adjust top */
        k = 2;
        while(sol.pivots[p_upper-2][0] < sol.pivots[p_lower-k][1]){ /* FLOPS */
            sol.yjoinpoints[i-k][1] = y[sol.pivots[p_upper-2][0]] + (y[sol.pivots[p_upper-2][0]]-sol.yjoinpoints[i-1][1])/(t[sol.pivots[p_upper-2][0]]-sol.tjoinpoints[i-1][0])*distance(sol.pivots[p_upper-2][0],sol.pivots[p_lower-k][1],t,store,check);
            k++;
#ifdef db6
            flop_counter += 5;
#endif
        }
    }
    else if(sol.pivots[p_lower-2][1] < sol.pivots[p_upper-2][0]){ /* or the other way around */
        k = 2;
        while(sol.pivots[p_lower-2][1] < sol.pivots[p_upper-k][0]){ /* FLOPS */
            sol.yjoinpoints[i-k][0] = y[sol.pivots[p_lower-2][1]] + (y[sol.pivots[p_lower-2][1]]-sol.yjoinpoints[i-1][0])/(t[sol.pivots[p_lower-2][1]]-sol.tjoinpoints[i-1][0])*distance(sol.pivots[p_lower-2][1],sol.pivots[p_upper-k][0],t,store,check);
            k++;
#ifdef db6
            flop_counter += 5;
#endif
        }
    }
#endif
    /* end of band expansion */
    
#ifdef db8b
    /* band expansion at any point */
    double height_added;
    for (p=0; p<=i-1; p++) {
        
        if(YSPAN(sol.yjoinpoints[p]) < min_height){
            //printf("at t=%f, SPAN fails (height %f)\n",sol.tjoinpoints[p][0],YSPAN(sol.yjoinpoints[p]));
            height_added = 0.5*(min_height - YSPAN(sol.yjoinpoints[p]));
            /* head in both directions and add height to the minimum until you reach a wide segment */
            
            /* First, head left */
            q = p-1;
            sol.yjoinpoints[q+1][0] -= height_added;
            sol.yjoinpoints[q+1][1] += height_added;
            while (q >= 0 && sol.tjoinpoints[q+1][0]-sol.tjoinpoints[q][0] < user_tmin){
                //printf("\tdistance = %f->%f -> add height at %d\n",sol.tjoinpoints[q][0],sol.tjoinpoints[q+1][0],q+1);
                q--;
                sol.yjoinpoints[q+1][0] -= height_added;
                sol.yjoinpoints[q+1][1] += height_added;
            }
            
            /* Second, head right */
            q = p;
            while (q+1 <= i-1 && sol.tjoinpoints[q+1][0]-sol.tjoinpoints[q][0] < user_tmin){
                //printf("\tdistance = %f->%f -> add height at %d\n",sol.tjoinpoints[q][0],sol.tjoinpoints[q+1][0],q);
                q++;
                sol.yjoinpoints[q][0] -= height_added;
                sol.yjoinpoints[q][1] += height_added;
            }
            //exit(1);
        }
        
        
    }
    
#endif
    
    /* TIMING stop */
#ifdef db7
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%f\t",cpu_time_used);
#endif  
        
    /* final storeage of pivot list */
    *Tout = (double *) malloc(i*sizeof(double));
	*Yout = (Range *) malloc(i*sizeof(Range));
    for(p = 0; p <= i-1; p++){
//      (*Tout)[p][0] = sol.tjoinpoints[p][0];
		(*Yout)[p][0] = sol.yjoinpoints[p][0];
        (*Tout)[p] = sol.tjoinpoints[p][0];
//      (*Tout)[p][1] = sol.tjoinpoints[p][0];
		(*Yout)[p][1] = sol.yjoinpoints[p][1];
        final_height += p > 0 ? 0.5*(sol.yjoinpoints[p][1]-sol.yjoinpoints[p][0]+sol.yjoinpoints[p-1][1]-sol.yjoinpoints[p-1][0])*(sol.tjoinpoints[p][0]-sol.tjoinpoints[p-1][0]) : 0;
	}
    
//    printf("%f\n",final_height); /* total area */
    final_height /= t[size-1]-t[0];
//    printf("%f\n",final_height); /* average height */
#ifdef db6
    printf("%ld\t%ld\t",flop_counter,flop_counter_2);
#endif
    
    free(sol.pivots); free(sol.interval_list); free(sol.status);
    for(p = 0; p < size; p++) free(store[p]); 
	free(store);
    for(p = 0; p < size; p++) free(check[p]); 
	free(check);
    free(sol.tjoinpoints); free(sol.yjoinpoints);
    free(encroach_pivot);
    for(p = 0; p < size; p++) free(potential_pivot[p]); 
	free(potential_pivot);
    
  	return i;   /* return the number of data join points */
}

int process_status(double *t, double *y, intRange *status, int size, double user_tmin, intRange *pivots, Range **store, intRange **check,  int *interval_list, intRange *encroach_pivot, int **potential_pivot, int p, int k, int iteration, int top, int counter){
    
    int interval_before, interval_next, temp_left, temp_right, temp_p, i, adjacent_flag, adjacent_left, adjacent_right;
    
#ifdef db1
    printf("   k=%d(%d.0), {p=%d}, [%d,%d] = [%f,%f] - Status %d (%d)\n",k,1,p,interval_list[k],interval_list[k+1],t[interval_list[k]],t[interval_list[k+1]],status[2*k][(iteration-1)%2],status[2*k+1][(iteration-1)%2]);
#endif
    
    /* if the interval has been banned from future splitting, then just transfer the points (no status required) */
    if(status[2*k][(iteration-1)%2] == 1){
        set_status(status,p-1,iteration,status[2*k][(iteration-1)%2],status[2*k+1][(iteration-1)%2]);
    }
    /* if the interval is a potential pivot shift */
    else if(status[2*k][(iteration-1)%2] >= 4){
                        
        if((k == counter || status[2*(k+1)][(iteration-1)%2] > 0) && (k == 0 || status[2*(p-2)][(iteration-1)%2] > 0)){
            
            /* extract information to see where we might be pivoting and which additional points should be considered */
            temp_right = status[2*k+1][(iteration-1)%2] < 0 ? interval_list[k]+(abs_value(status[2*k+1][(iteration-1)%2]))%(interval_list[k+1]-interval_list[k]+1) : (status[2*k][(iteration-1)%2] == 5 ? status[2*k+1][(iteration-1)%2] : interval_list[k+1]);
            temp_left = status[2*k+1][(iteration-1)%2] < 0 ? interval_list[k]+((abs_value(status[2*k+1][(iteration-1)%2]))-(temp_right-interval_list[k]))/(interval_list[k+1]-interval_list[k]+1) : (status[2*k][(iteration-1)%2] == 5 ? interval_list[k] : status[2*k+1][(iteration-1)%2]);
            
            /* first check to see if such a pivot shift is still valid given other changes to the algorithm */
            if((status[2*k][(iteration-1)%2] == 4 && slope_compare(pivots[p-2][top],pivots[p-1][top],pivots[p-1][top],temp_left,t,y,store,check,top) && (p <= 2 || distance_check(pivots[p-3][top],pivots[p-2][top],t,store,check,user_tmin) || slope_compare(pivots[p-2][top],temp_left,pivots[p-3][top],pivots[p-2][top],t,y,store,check,top)) && distance_check(temp_left,temp_right,t,store,check,2*user_tmin) && (status[2*k+1][(iteration-1)%2] >= 0 || k >= counter || slope_compare(interval_list[k+1],interval_list[k+2],temp_right,interval_list[k+1],t,y,store,check,top))) || (status[2*k][(iteration-1)%2] == 5 && slope_compare(temp_right,interval_list[k+1],interval_list[k+1],interval_list[k+2],t,y,store,check,top) && distance_check(temp_left,temp_right,t,store,check,2*user_tmin) && (k+1 == counter || distance_check(interval_list[k+2],interval_list[k+3],t,store,check,user_tmin) || slope_compare(interval_list[k+2],interval_list[k+3],temp_right,interval_list[k+2],t,y,store,check,top))  && (status[2*k+1][(iteration-1)%2] >= 0 || p <= 2 || slope_compare(interval_list[k],temp_left,pivots[p-2][top],interval_list[k],t,y,store,check,top)))){

                /* check to see if this is the interval we should be processing at this time: first, we travel left to see if we encounter first a unfinished interval, a finished interval or a pivot shift interval with a wider distance, then we travel right to look for the same thing */
                adjacent_flag = 0; i = 2;
                while(adjacent_flag == 0){
                    
                    if(status[2*(p-i)][iteration%2] == 1 || i > p)  adjacent_flag = 2;
                    else if(status[2*(p-i)][iteration%2] <= 0)  adjacent_flag = 1;
                    else{ 
                        adjacent_right = status[2*(p-i)+1][(iteration)%2] < 0 ? pivots[p-i][top]+(abs_value(status[2*(p-i)+1][(iteration)%2]))%(pivots[p-i+1][top]-pivots[p-i][top]+1) : (status[2*(p-i)][(iteration)%2] == 5 ? status[2*(p-i)+1][(iteration)%2] : pivots[p-i+1][top]);
                        adjacent_left = status[2*(p-i)+1][(iteration)%2] < 0 ? pivots[p-i][top]+((abs_value(status[2*(p-i)+1][(iteration)%2]))-(adjacent_right-pivots[p-i][top]))/(pivots[p-i+1][top]-pivots[p-i][top]+1) : (status[2*(p-i)][(iteration)%2] == 5 ? pivots[p-i][top] : status[2*(p-i)+1][(iteration)%2]);
                        
                        if(distance_compare(temp_left,temp_right,adjacent_left,adjacent_right,t,store,check))   adjacent_flag = 1;
                    }
                    i++;
                }
                if(adjacent_flag == 2){
                    adjacent_flag = 0; i = 1;
                    while(adjacent_flag == 0){
                        if(status[2*(k+i)][(iteration-1)%2] == 1 || k+i+1 > counter)  adjacent_flag = 2;
                        else if(status[2*(k+i)][(iteration-1)%2] <= 0)  adjacent_flag = 1;
                        else{
                            adjacent_right = status[2*(k+i)+1][(iteration-1)%2] < 0 ? interval_list[(k+i)]+(abs_value(status[2*(k+i)+1][(iteration-1)%2]))%(interval_list[(k+i)+1]-interval_list[(k+i)]+1) : (status[2*(k+i)][(iteration-1)%2] == 5 ? status[2*(k+i)+1][(iteration-1)%2] : interval_list[(k+i)+1]);
                            adjacent_left = status[2*(k+i)+1][(iteration-1)%2] < 0 ? interval_list[(k+i)]+((abs_value(status[2*(k+i)+1][(iteration-1)%2]))-(adjacent_right-interval_list[(k+i)]))/(interval_list[(k+i)+1]-interval_list[(k+i)]+1) : (status[2*(k+i)][(iteration-1)%2] == 5 ? interval_list[(k+i)] : status[2*(k+i)+1][(iteration-1)%2]);
                            if(distance_compare(temp_left,temp_right,adjacent_left,adjacent_right,t,store,check)) adjacent_flag = 1;
                        }
                        i++;
                    }
                    
                    if(adjacent_flag == 2){
                        if(status[2*k][(iteration-1)%2] == 5){  
                            
                            if(temp_left != interval_list[k]){
                                temp_p = convex(t,y,interval_list[k],temp_left+1,pivots,2,top,p-1,store,check);
                                for(i = p-1; i <= p+temp_p; i++)  set_status(status,i,iteration,1,0);
                                p += temp_p;  
                            }
#ifdef db1
                            printf("      Successful Shift - Status 0\n");
#endif
                            set_status(status,k+1,iteration-1,1,0);
                            set_status(status,p-1,iteration,-5,interval_list[k+1]);
                            interval_list[k+1] = temp_right;
                        }
                        else if(status[2*k][(iteration-1)%2] == 4){
#ifdef db1
                            printf("      Successful Shift - Status 0\n");
#endif
                            set_status(status,p-2,iteration,1,0);  
                            set_status(status,p-1,iteration,-4,pivots[p-1][top]);     
                            pivots[p-1][top] = temp_left;
                            
                            if(temp_right != interval_list[k+1]){
                                temp_p = convex(t,y,temp_right,interval_list[k+1]+1,pivots,2,top,p,store,check);
                                for(i = p; i <= p+temp_p-1; i++)   set_status(status,i,iteration,1,0);
                                p += temp_p;
                            } 
                        }                         
                    }
                    else{
                        set_status(status,p-1,iteration,status[2*k][(iteration-1)%2],status[2*k+1][(iteration-1)%2]);
                    }
                }
                else  set_status(status,p-1,iteration,status[2*k][(iteration-1)%2],status[2*k+1][(iteration-1)%2]);
            }
            else{
                set_status(status,p-1,iteration,0,status[2*k+1][(iteration-1)%2]);
                adjacent_flag = 1;
            }
        }
        else{
            set_status(status,p-1,iteration,status[2*k][(iteration-1)%2],status[2*k+1][(iteration-1)%2]);
        }
    }
    /* else if the interval has only two points or is already less than 2*user_tmin wide */
    else if(k != counter && k != 0 && (interval_list[k] + 1 == interval_list[k+1] || !distance_check(interval_list[k],interval_list[k+1],t,store,check,2*user_tmin))){
        set_status(status,p-1,iteration,1,0);
    }
    /* if none of the above conditions apply then process the interval */
    else{        
        interval_before = interval_list[k] == 0 ? -1 : pivots[p-2][top];
        interval_next = interval_list[k+1] == size-1 ? -1 : interval_list[k+2];
        
        temp_p = p;
        p = process_interval(t, y, status, size, user_tmin, pivots, interval_list[k], interval_list[k+1], interval_next, interval_before, p, k, iteration, store, check, interval_list, encroach_pivot, potential_pivot, top);
        
        /* if progres was made, and the current pair had been shifted, then the adjacent interval is opened */
        if(temp_p != p){
            if(status[2*k][(iteration-1)%2] == -5)  set_status(status,k+1,iteration-1,0,0);
            if(status[2*k][(iteration-1)%2] == -4)  set_status(status,temp_p-2,iteration,0,0);
        }
    }
    
    return p;
}

int process_interval(double *t, double *y, intRange *status, int size, double user_tmin, intRange *pivots, int interval_start, int interval_finish, int interval_next, int interval_before, int p, int k, int iteration, Range **store, intRange **check, int *interval_list, intRange *encroach_pivot, int **potential_pivot, int top){
    
    int q, e_shift, e = 1, e_left, e_right, temp_e_left, temp_e_right, j, i, pivot_count = 0, temp_interval_next, temp_interval_before, encroach_change = 0, flag_left, flag_right, original_start, original_finish, temp_info, temp_pivot_count, slope_change, slope_before, slope_after;   
    double  temp_difference;
        
    /* find pivots given the original values of interval_start and interval_finish */
    flag_left = interval_before == -1 || distance_check(interval_before,interval_start,t,store,check,user_tmin);
    flag_right = interval_next == -1 || distance_check(interval_finish,interval_next,t,store,check,user_tmin);
#ifdef db1
    printf("            calling find_next_pivot on [%d,%d]\n",interval_start,interval_finish);
#endif
    temp_pivot_count = find_next_pivot(t, y, pivots, p+pivot_count, user_tmin, interval_start, interval_finish,interval_next, interval_before, store, check, encroach_pivot, potential_pivot, flag_left, flag_right, e, top);
    
    /* if everything is fine and we obtained some sort of sight_pivot satisfying user_tmin condition, then increase pivot count and finish this function - set the status of each interval to 0 since they are now new candidates for splitting in future iterations */
    if(temp_pivot_count >= 0){
        pivot_count += temp_pivot_count;
        for(q = p-1; q <= p+temp_pivot_count-1; q++)   set_status(status,q,iteration,0,0);
    }
    /* if the temp_pivot_count is negative, it means we found a potential pivot, but it was too close to the start of the interval to be of any good - but may be useful again to make the interval narrower to give a better region of viability for a future search */
    else{
        /* make a list of pivots for considering slope information
         - the first index represents how far along into the interval we are (note that the first entry is not actually in the interval itself, and the second entry is the start/finish of the interval itself
         - the second index denotes 0 for left, 1 for right */
        encroach_pivot[0][0] = pivots[p+pivot_count-1][top] == 0 ? -1 : pivots[p+pivot_count-2][top];
        encroach_pivot[1][0] = pivots[p+pivot_count-1][top];
        encroach_pivot[0][1] = interval_next; 
        encroach_pivot[1][1] = interval_finish;
        original_start = encroach_pivot[1][0];
        original_finish = encroach_pivot[1][1];
        e = 2;
        
        /* while we are still making valid progress into the interval */
        while(encroach_change == 0){
            
            /* extract the new points we might be interested in using as the start and end of the new interval - the find next pivot function coded the new interval start and interval finish in the output - note that these are all integer calculations */
            encroach_pivot[e][1] = encroach_pivot[e-1][0]+(abs_value(temp_pivot_count))%(encroach_pivot[e-1][1]-encroach_pivot[e-1][0]+1);
            encroach_pivot[e][0] = encroach_pivot[e-1][0]+((abs_value(temp_pivot_count))-(encroach_pivot[e][1]-encroach_pivot[e-1][0]))/(encroach_pivot[e-1][1]-encroach_pivot[e-1][0]+1);
                        
            /* if all we have is tight_pivots - we only want to pick one to readjust the interval - we pick the one which is closest vertically to the line connecting the interval_start and interval_finish - if one of them has been reset has part of the process (or never chosen in the first place - then we have an easy decision, just pick that one. Otherwise, we must select between them and reset the other */
            if(encroach_pivot[e][0] != encroach_pivot[e-1][0] && encroach_pivot[e][1] != encroach_pivot[e-1][1]){
                
                e_shift = 0;
                while(encroach_pivot[e-2-e_shift][0]==encroach_pivot[e-1][0])   e_shift++;
                if(e_shift == e-2){
                    e_shift = encroach_pivot[1][1];
                    temp_e_left = encroach_pivot[1][0];
                }
                else{
                    temp_e_left = encroach_pivot[e-2-e_shift][0];
                    e_shift = encroach_pivot[e-1][0];
                }
//                printf("LEFT  %d,%d - %d,%d\n",temp_e_left,e_shift,encroach_pivot[e-1][0],encroach_pivot[e][0]);
                temp_difference = slope(temp_e_left,e_shift,t,y,store,check,top) - slope(encroach_pivot[e-1][0],encroach_pivot[e][0],t,y,store,check,top);
                
                e_shift = 0;
                while(encroach_pivot[e-2-e_shift][1]==encroach_pivot[e-1][1])   e_shift++;
                temp_e_right = encroach_pivot[e-2-e_shift][1];
                if(e_shift == e-2){
                    e_shift = encroach_pivot[1][0];
                    temp_e_right = encroach_pivot[1][1];
                }
                else{
                    temp_e_right = encroach_pivot[e-2-e_shift][1];
                    e_shift = encroach_pivot[e-1][1];
                }
//                printf("RIGHT %d,%d - %d,%d\n",encroach_pivot[e][1],encroach_pivot[e-1][1],e_shift,temp_e_right);
                temp_difference -= slope(encroach_pivot[e][1],encroach_pivot[e-1][1],t,y,store,check,top) - slope(e_shift,temp_e_right,t,y,store,check,top);
                
//                temp_difference = (slope(encroach_pivot[1][0],encroach_pivot[1][1],t,y,store,check,top) - slope(encroach_pivot[1][0],encroach_pivot[e][0],t,y,store,check,top)) - (slope(encroach_pivot[e][1],encroach_pivot[1][1],t,y,store,check,top) - slope(encroach_pivot[1][0],encroach_pivot[1][1],t,y,store,check,top));
//                temp_difference = slope(encroach_pivot[1][0],encroach_pivot[1][1],t,y,store,check,top)*distance(encroach_pivot[e][0],encroach_pivot[e][1],t,store,check) - y[encroach_pivot[e][0]] + y[encroach_pivot[e][1]];
#ifdef db6
                flop_counter += 3;
#endif
                
                if(temp_difference < 0)  encroach_pivot[e][1] = encroach_pivot[e-1][1];
                else  encroach_pivot[e][0] = encroach_pivot[e-1][0];
            }
            
            /* find the next point along the progression that is different so we can use it for slope purposes */
            temp_e_right = encroach_pivot[e-2][1];
            if(encroach_pivot[e-2][1] == encroach_pivot[e-1][1]){
                e_shift = 0;
                while(encroach_pivot[e-2-e_shift][1]==encroach_pivot[e-1][1])   e_shift++;
                temp_e_right = encroach_pivot[e-2-e_shift][1];
            }
            temp_e_left = encroach_pivot[e-2][0];
            if(encroach_pivot[e-2][0] == encroach_pivot[e-1][0]){
                e_shift = 0;
                while(encroach_pivot[e-2-e_shift][0]==encroach_pivot[e-1][0])   e_shift++;
                temp_e_left = encroach_pivot[e-2-e_shift][0];
            }
            
            
            /* calculate the slope condition on either side of the interval */
            slope_before = encroach_pivot[e][0] != encroach_pivot[e-1][0] && temp_e_left != -1 ? slope_compare(encroach_pivot[e-1][0],encroach_pivot[e][0],temp_e_left,encroach_pivot[e-1][0],t,y,store,check,top) : 1;
            slope_after = encroach_pivot[e][1] != encroach_pivot[e-1][1] && temp_e_right != -1 ? slope_compare(encroach_pivot[e-1][1],temp_e_right,encroach_pivot[e][1],encroach_pivot[e-1][1],t,y,store,check,top) : 1;
            
#ifdef db1
            printf("          slope condition before (using %d,%d,%d): %d or %d\n",temp_e_left,encroach_pivot[e-1][0],encroach_pivot[e][0],slope_before == 1,encroach_pivot[e-1][0]==encroach_pivot[e][0]);
            printf("          slope condition after (using %d,%d,%d): %d or %d\n",encroach_pivot[e][1],encroach_pivot[e-1][1],temp_e_right,slope_after == 1,encroach_pivot[e-1][1]==encroach_pivot[e][1]); 
#endif
            
            /* if we fail the slope condition at the end of the interval and it is NOT the first time we have moved in on the right, then we can just simply bypass the problem by connecting over top of the slope error back to the last time the slope worked fine! - the function immediately below describes the case on the other side */
            slope_change = 0;
            if(slope_after == 0 && encroach_pivot[e-1][1] != encroach_pivot[1][1]){
                j = 2;
                while(encroach_pivot[e-j+1][1] != encroach_pivot[1][1] && slope_compare(encroach_pivot[e][1],encroach_pivot[e-j+1][1],encroach_pivot[e][1],encroach_pivot[e-j][1],t,y,store,check,top))   j++;
                
                for(i = e-j+2; i <= e-1; i++)   encroach_pivot[i][1] = encroach_pivot[e][1];
                
                if(encroach_pivot[e-j+1][1] != encroach_pivot[1][1] || (encroach_pivot[0][1] == -1 || slope_compare(encroach_pivot[1][1],encroach_pivot[0][1],encroach_pivot[e][1],encroach_pivot[1][1],t,y,store,check,top)))   slope_after = 1;
                else   slope_change = 1;
            }
            if(slope_before == 0 && encroach_pivot[e-1][0] != encroach_pivot[1][0]){
                
                j = 2;
                while(encroach_pivot[e-j+1][0] != encroach_pivot[1][0] && slope_compare(encroach_pivot[e-j][0],encroach_pivot[e][0],encroach_pivot[e-j+1][0],encroach_pivot[e][0],t,y,store,check,top))   j++;
                
                for(i = e-j+2; i <= e-1; i++)   encroach_pivot[i][0] = encroach_pivot[e][0];
                
                if(encroach_pivot[e-j+1][0] != encroach_pivot[1][0] || (encroach_pivot[0][0] == -1 || slope_compare(encroach_pivot[1][0],encroach_pivot[e][0],encroach_pivot[0][0],encroach_pivot[1][0],t,y,store,check,top)))   slope_before = 1;
                else   slope_change = 1;
            }
                        
            /* in order to proceed here, we must have a valid reason to continue on both sides of the interval - this is obtained in one of three ways
             - if the slope condition before is satisfied OR
             - the start or end of the interval hasn't changed since the last processing OR
             - we are dealing with the first or last interval in the data set */
            if(encroach_change == 0 && slope_change == 0 && (slope_before == 1 || encroach_pivot[e-1][0] == encroach_pivot[e][0] || temp_e_left == -1) && (slope_after == 1 || encroach_pivot[e-1][1] == encroach_pivot[e][1] || temp_e_right == -1)){
                /* set the next interval for processing */
                interval_start = encroach_pivot[e][0];
                interval_finish = encroach_pivot[e][1];
#ifdef db1
                printf("       (%d.%d), {p=%d}, [%d,%d] = [%f,%f] - ",1,e-1,p+pivot_count,interval_start,interval_finish,t[interval_start],t[interval_finish]);
#endif
                /* if the interval is too narrow, then exit all loops */
                if(!distance_check(interval_start,interval_finish,t,store,check,user_tmin) || interval_start+1 == interval_finish){
                    set_status(status,p-1,iteration,1,0);
                    encroach_change = 1;
#ifdef db1
                    printf("too narrow - Status %d (with info %d)\n",status[2*(p-1)][iteration%2],status[2*(p-1)+1][iteration%2]);
#endif                 
                }
                /* if the interval is wide enough to accomodate a potential sight pivot point, then continue */
                else{
#ifdef db1
                    printf("process\n");
#endif
                    /* find the most recent different pivot on the left and the right for slope checking purposes - the variable flag_left/flag_right indicates whether or not the slope condition must be checked on the respective sides - if the adjacent interval to the left is wide, then we do not, for example */
                    j = e-1;
                    while(encroach_pivot[j][1] == encroach_pivot[e][1])   j--;
                    temp_interval_next = encroach_pivot[j][1];
                    flag_right = j == 0 && (encroach_pivot[0][1] == -1 || distance_check(encroach_pivot[1][1],encroach_pivot[0][1],t,store,check,user_tmin)) ? 1 : 0;
                    
                    j = e-1;
                    while(encroach_pivot[j][0] == encroach_pivot[e][0])   j--;
                    temp_interval_before = encroach_pivot[j][0];
                    flag_left = j == 0 && (encroach_pivot[0][0] == -1 || distance_check(encroach_pivot[0][0],encroach_pivot[1][0],t,store,check,user_tmin)) ? 1 : 0;
                    
                    /* process the interval using the new smaller interval and the just calculated values of the pivots before and after the current interval */
#ifdef db1
                    printf("            calling find_next_pivot on [%d,%d]\n",interval_start,interval_finish);
#endif
                    temp_pivot_count = find_next_pivot(t, y, pivots, p+pivot_count, user_tmin, interval_start, interval_finish, temp_interval_next, temp_interval_before, store, check, encroach_pivot, potential_pivot, flag_left, flag_right, e, top);
                    
                    encroach_change = temp_pivot_count > 0 ? 1 : 0;
                }
            }
            /* if instead we fail one of the conditions indicated above, but are still allowed to proceed */
            else if(encroach_change == 0){
#ifdef db1
                printf("       (%d.%d), {p=%d}, [%d,%d] = [%f,%f] (info NULL) - FAILED - Status ",1,e-1,p+pivot_count,encroach_pivot[e][0],encroach_pivot[e][1],t[encroach_pivot[e][0]],t[encroach_pivot[e][1]]);
                printf("\n");
#endif
                /* if the slope fails on the left side of the interval */
                if(slope_before == 0){
                    encroach_change = 1;
                    /* if the point on which it fails has been used before, then stop the potential looping */
                    if(status[2*(p-1)+1][(iteration-1)%2] == encroach_pivot[e][0]){
                        set_status(status,p-1,iteration,1,encroach_pivot[e][0]);
                        encroach_change = 1;
                    }
                    /* otherwise, allow a conditional shift that will only be allowed if good pivots are eventually found - note that we must be very careful here - since we have conditionally included pivots on the right side of the interval, we must ensure that such a shift would not excise those points from the band, if so, they must also be included for a potential shift */
                    else{
                        
                        /* start by putting the slope-failing point as info - and check below if this should change */
                        set_status(status,p-1,iteration,4,encroach_pivot[e][0]);
                        
                        /* for each pivot that has been conditionally included on the right */
                        for(i = e-1; i >= 2; i--){
                            e_shift = 1;
                            while(encroach_pivot[i-e_shift][1] == encroach_pivot[i][1])   e_shift++;
                            temp_e_right = encroach_pivot[i-e_shift][1];
                            
                            /* if one of the points added on the right would be excised from the band if a shift was made onto the slope failing point on the left */
                            if(temp_e_right == -1 || !slope_compare(encroach_pivot[e][0],encroach_pivot[i][1],encroach_pivot[i][1],temp_e_right,t,y,store,check,top)){
                                if(distance_check(encroach_pivot[e][0],encroach_pivot[i][1],t,store,check,user_tmin)){
                                    /* this point causes a potential problem with any future shift and must be accounted for when performing the shift */
                                    temp_info = -((encroach_pivot[1][1]-encroach_pivot[1][0]+1)*(encroach_pivot[e][0]-encroach_pivot[1][0])+(encroach_pivot[i][1]-encroach_pivot[1][0]));
                                    set_status(status,p-1,iteration,4,temp_info);
                                    i = 0;
                                }
                                else{
                                    /* no shift will ever be allowed here */
                                    set_status(status,p-1,iteration,1,0);
                                    i = 0;
                                }
                            }
                        }
                    }  
                }
                else if(slope_after == 0){
                    
                    encroach_change = 1;
                    /* if the point on which it fails has been used before, then stop the potential looping */
                    if(status[2*(p-1)+1][(iteration-1)%2] == encroach_pivot[e][1]){
                        set_status(status,p-1,iteration,1,encroach_pivot[e][1]);
                        encroach_change = 1;
                    }
                    else{
                        /* start by putting the slope-failing point as info - and check below if this should change */
                        set_status(status,p-1,iteration,5,encroach_pivot[e][1]);
                        /* for each pivot that has been conditionally included on the left */
                        for(i = e-1; i >= 2; i--){
                            e_shift = 1;
                            while(encroach_pivot[i-e_shift][0] == encroach_pivot[i][0])   e_shift++;
                            temp_e_left = encroach_pivot[i-e_shift][0];
                            
                            /* if one of the points added on the left would be excised from the band if a shift was made onto the slope failing point on the right */
                            if(temp_e_left == -1 || !slope_compare(temp_e_left,encroach_pivot[i][0],encroach_pivot[i][0],encroach_pivot[e][1],t,y,store,check,top)){
                                if(distance_check(encroach_pivot[i][0],encroach_pivot[e][1],t,store,check,user_tmin)){
                                    /* this point causes a potential problem with any future shift and must be accounted for when performing the shift */
                                    temp_info = -((encroach_pivot[1][1]-encroach_pivot[1][0]+1)*(encroach_pivot[i][0]-encroach_pivot[1][0])+(encroach_pivot[e][1]-encroach_pivot[1][0]));
                                    set_status(status,p-1,iteration,5,temp_info);
                                    i = 0;
                                }
                                else{
                                    /* no shift will ever be allowed here */
                                    set_status(status,p-1,iteration,1,0);
                                    i = 0;
                                }
                            }
                        }
                    }  
                }
                else{
#ifdef db1
                    printf("1 AHHHH - line 691\n"); exit(1);
#endif
                    set_status(status,p-1,iteration,1,0);
                    encroach_change = 1;
                } /* why would there be an else? */
            }
            e++;
        }
        
        /* bookkeeping below */
        e_left = 0; e_right = 0;
        if(temp_pivot_count > 0){     
            pivots[p-1][top] = encroach_pivot[1][0];
            for(q = p-1; q <= p-1+(e-1)+temp_pivot_count+1-2; q++)   set_status(status,q,iteration,0,0);
            for(q = 2; q <= e-1; q++){
                if(encroach_pivot[q][0] != encroach_pivot[q-1][0] && encroach_pivot[q][0] != pivots[p+pivot_count-1][top]){
                    /* shift the pivots that have been entered above */
                    for(e_shift = p+pivot_count+temp_pivot_count; e_shift >= p+pivot_count; e_shift--){
                        pivots[e_shift+e_left+1][top] = pivots[e_shift+e_left][top];
                    }
                    pivots[p+pivot_count+(e_left++)][top] = encroach_pivot[q][0];
                }
            }
            for(q = e-1; q >= 2; q--){
                if(encroach_pivot[q][1] != encroach_pivot[q-1][1]){
                    pivots[p+pivot_count+temp_pivot_count+e_left+(e_right++)][top] = encroach_pivot[q][1];
                }
            }
        }
        else   temp_pivot_count = 0;
        
        pivot_count += temp_pivot_count+e_left+e_right;
    }
    
    if(pivot_count > 0){
        /* if interval has two points in it, set that status to 1, there is nothing more we can do to this interval */
        for(q = p; q < p+pivot_count; q++){
            if(pivots[q][top] == 1 + pivots[q-1][top])   set_status(status,q-1,iteration,1,0);
        }
        
        if(interval_finish == 1 + pivots[p+pivot_count-1][top])   set_status(status,p+pivot_count-1,iteration,1,0);
    }
    
    /* this information below is to stop cycling */
    if(pivot_count == 0 && status[2*k][(iteration-1)%2] == 0 && status[2*k+1][(iteration-1)%2] == status[2*(p-1)+1][(iteration)%2]){
        set_status(status,p-1,iteration,1,status[2*k+1][(iteration-1)%2]);
    }
    
    return p+pivot_count;
}

int find_next_pivot(double *t, double *y, intRange *pivots, int p, double user_tmin, int interval_start, int interval_finish, int interval_next, int interval_before, Range **store, intRange **check, intRange *encroach_pivot, int **potential_pivot, int flag_left, int flag_right, int e, int top){
    
    int i, j = 0, k = 0, q, pivot_counter[4] = {0,0,0,0}, back_track = 0, back_counter = 0, slope_condition_left, slope_condition_right, distance_condition_left, distance_condition_right, left_right_condition, max_back_index, test_slope_point, temp_a, temp_b, test_finish;
    double interval_width;
    int *back_encroach; back_encroach = (int*) calloc(e,sizeof(int));
    
    /* potential pivots are to be stored in the array 'potential_pivot', the first index refers to the placement of the potential pivots in the list. The second index refers to the type of potential pivots (i.e. 0 is Region I pivots, 1 is Region IIa/IIIa pivots, and 2 is Region IIb/IIIb pivots */ 
    
    /* if we have found Region III points on both sides of the original interval, then we may have to do some back-checking on previous combinations of interval starts and finishes */
    if(encroach_pivot[e][0] != encroach_pivot[1][0] && encroach_pivot[e][1] != encroach_pivot[1][1]){
        /* if the interval start has changed, we will have to test the new interval start with all previous points used as interval finishes as part of the same loop process, else if the interval finish has changed, we will have to test will previous points used as interval starts */
        if(encroach_pivot[e][0] != encroach_pivot[e-1][0]){
            back_encroach[0] = encroach_pivot[e][1]; back_track = 2;
            for(j = e-1; j >= 0; j--){
                if(encroach_pivot[j][1] != encroach_pivot[e][1] && encroach_pivot[j][1] != back_encroach[back_counter]){
                    back_encroach[++back_counter] = encroach_pivot[j][1];
                }
            }
        }
        else{
            for(j = 0; j <= e-1; j++){
                if(encroach_pivot[j][0] != encroach_pivot[e][0] && (back_counter == 0 || encroach_pivot[j][0] != back_encroach[back_counter-1])){
                    back_encroach[back_counter++] = encroach_pivot[j][0];
                }
            }
            back_encroach[back_counter] = encroach_pivot[e][0]; back_track = 1;
        }
    }
    
    /* test_slope contains the currently largest slope being considered - any future potential pivots must have a slope greater than this, or the point used to calculate test_slope would be blocking the test point from the start */
    test_slope_point = interval_start + 1;
    
    /* test every point in the interval */
    for(q = interval_start+1; q < interval_finish; q++){
        
        /* first, we test all existing potential pivots to see if the current point would cause a blockage to the end of the interval - if this is the case, then remove the pivot from the given list and shift all other pivots accordingly */
        test_finish = 2 + (back_track == 1);
        for(i = 0; i <= test_finish; i++){ /* for each type of pivot */
            for(j = 0; j <= pivot_counter[i]-1; j++){ /* for each pivot in that list */
                /* if there is a point above the line connecting our potential pivot to the finish, then our point is not a potential pivot, so should be removed */
                if(!slope_compare(potential_pivot[j][i],interval_finish,q,interval_finish,t,y,store,check,top)){
                    for(k = j; k <= pivot_counter[i]-1; k++)   potential_pivot[k][i] = potential_pivot[k+1][i];
                    pivot_counter[i]--; j--;
                }
            }
        }
        
        /* we need to be checking to see if these points are valid - this is the same as the previous check, but this time, we are checking to see if line of sight fails to the point stored as the corresponding end point */
        if(back_track == 2){  /* i is just fixed to be 3 here */  
            for(j = 0; j <= pivot_counter[3]-1; j++){ /* for each pivot in that list */
                if(!slope_compare(potential_pivot[j][3],back_encroach[potential_pivot[j][4]],q,back_encroach[potential_pivot[j][4]],t,y,store,check,top)){
                    for(k = j; k <= pivot_counter[3]-1; k++){   
                        potential_pivot[k][3] = potential_pivot[k+1][3];  
                        potential_pivot[k][4] = potential_pivot[k+1][4];
                    }
                    pivot_counter[3]--; j--;
                }
            }
        }
        
        /* next, if the point has a clear line of sight to the start, then it will be added to one of the three lists */
        if(slope_compare(interval_start,test_slope_point,interval_start,q,t,y,store,check,top)){             
            test_slope_point = q;
            /* slope condition and distance conditions help us classify the type of point we are processing - the flag_left and flag_right are passed from the process_interval and tells us the distance of the adjacent interval - if they are wide, then we need not consider any slope condition, likewise, if we are dealing with the first or last interval, we also need not consider adjancent slope conditions */
            slope_condition_left = (flag_left == 1 || interval_before == -1 || slope_compare(interval_start,q,interval_before,interval_start,t,y,store,check,top)) ? 1 : 0;
            slope_condition_right = (flag_right == 1 || interval_next == -1 || slope_compare(interval_finish,interval_next,q,interval_finish,t,y,store,check,top)) ? 1 : 0;
            distance_condition_left = interval_before == -1 || distance_check(interval_start,q,t,store,check,user_tmin);
            distance_condition_right = interval_next == -1 || distance_check(q,interval_finish,t,store,check,user_tmin);
            
            /* if all four conditions are satisfied, we have a Region I pivot, else if one of the conditions fails, then the point is placed according to which side of the horizontal half-way place it is located*/
            if(slope_condition_left + slope_condition_right + distance_condition_left + distance_condition_right == 4){
                potential_pivot[pivot_counter[0]++][0] = q;
            }
            else{
                left_right_condition = 2 - distance_compare(interval_start,q,q,interval_finish,t,store,check);
                potential_pivot[pivot_counter[left_right_condition]++][left_right_condition] = q;
            }
            
            
            /* if we are also checking through all previous end points, there is more to check */
            if(back_track == 2){  
                
                for(i = 0; i <= back_counter-1; i++){
                    /* check conditions for each of the potential interval finishes - if any of them satisfy the conditions, then add this point to the list, with the additional information stored in [*][4] of which end point was valid */
                    slope_condition_right = back_encroach[i+1] == -1 || distance_check(back_encroach[i],back_encroach[i+1],t,store,check,user_tmin) || slope_compare(back_encroach[i],back_encroach[i+1],q,back_encroach[i],t,y,store,check,top);
                    distance_condition_right = distance_check(q,back_encroach[i],t,store,check,user_tmin);
                    
                    /* if all four conditions are satisfied with the different end points */
                    if(slope_condition_left + slope_condition_right == 2){
                        if(distance_condition_left + distance_condition_right == 2){
                            potential_pivot[pivot_counter[3]][3] = q;
                            potential_pivot[pivot_counter[3]++][4] = i;
                        }
                        i = back_counter; /* stop the loop */
                    }
                }
            }
            else if(back_track == 1){
                for(i = back_counter; i >= 1; i--){
                    /* check conditions on the left for different potential interval starts */
                    slope_condition_left = back_encroach[i-1] == -1 || distance_check(back_encroach[i-1],back_encroach[i],t,store,check,user_tmin) || slope_compare(back_encroach[i],q,back_encroach[i-1],back_encroach[i],t,y,store,check,top);
                    distance_condition_left = distance_check(back_encroach[i],q,t,store,check,user_tmin);
                    /* if all four conditions are satisfied with the different end points */
                    if(slope_condition_left + slope_condition_right == 2){
                        if(distance_condition_left + distance_condition_right == 2){
                            potential_pivot[pivot_counter[3]][3] = q;
                            potential_pivot[pivot_counter[3]++][4] = i;
                        }
                        i = 0; /* stop the loop */
                    }
                }
            }
        } 
    }
    
    /* if there are no Region I points with the current interval start or finish, but there are with previous incarnations, then move these points into the good list so they can be processed and change the list of encroach pivots */
    if(pivot_counter[0] == 0 && pivot_counter[3] > 0){
        /* first, find the points with the furthest along start or end point - we can only choose pivots from the list of points with the same start/end point */
        max_back_index = 0;
        for(i = 0; i <= pivot_counter[3]-1; i++){
            if(potential_pivot[i][4] > max_back_index)   max_back_index = potential_pivot[i][4];
        }
        /* now, we put the points corresponding to that index into the potential_pivot[*][0] for use below */
        for(i = 0; i <= pivot_counter[3]-1; i++){
            if(potential_pivot[i][4] == max_back_index)   potential_pivot[pivot_counter[0]++][0] = potential_pivot[i][3];
        }
        /* also adjust the encroach pivots to remove the ones that are not needed */
        i = e;
        while(encroach_pivot[i][back_track-1] != back_encroach[max_back_index]){
            encroach_pivot[i--][back_track-1] = back_encroach[max_back_index];
        }
    }
    
    /* now that we have tested all points, we must decide which information to pass back - if there are any region I pivots (or pivots that would be region I pivots if different interval start or finishes are chosen), then choose among those points first. Otherwise, pass back the most restrictive slope points from regions II/IIIa, and II/IIIb */
    q = 0;
    if(pivot_counter[0] > 0){
        /* initialize by setting all status to 0 (we will use potential_pivot[*][1] as a tracking array - if the status is 0, the point is still eligible to be chosen, 1 the point has been chosen, and 2 the point has been rejected from further consideration */
        for(i = 0; i <= pivot_counter[0]-1; i++)   potential_pivot[i][1] = 0;
        int largest_pivot, eligible_points = pivot_counter[0];
        
        /* continue this procedure while there are still points yet to be classified */
        while(eligible_points > 0){            
            largest_pivot = -1;
            
            /* initialize by finding the highest point of all remaining eligible points */
            for(i = 0; i <= pivot_counter[0]-1; i++){
                if(potential_pivot[i][1] == 0 && (largest_pivot == -1 || y[potential_pivot[i][0]] > y[potential_pivot[largest_pivot][0]]))   largest_pivot = i;
            }
            potential_pivot[largest_pivot][1] = 1; eligible_points--;
            
            /* first, test each point to the left of the most recently added point and remove from consideration any point which fails a line of sight or distance test */
            i = largest_pivot - 1;
            if(largest_pivot >= 1)  test_slope_point = potential_pivot[largest_pivot-1][0];
            
            while(i >= 0 && potential_pivot[i][1] != 1){
                distance_condition_left = distance_check(potential_pivot[i][0],potential_pivot[largest_pivot][0],t,store,check,user_tmin);
                slope_condition_left = slope_compare(potential_pivot[i][0],potential_pivot[largest_pivot][0],test_slope_point,potential_pivot[largest_pivot][0],t,y,store,check,top);
                
                if(distance_condition_left + slope_condition_left < 2){  potential_pivot[i][1] = 2; eligible_points--; }
                if(slope_condition_left == 1)  test_slope_point = potential_pivot[i][0];
                i--;
            }
            
            /* then check all the points to the right of the most recently added point - same test as above */
            i = largest_pivot + 1;
            test_slope_point = largest_pivot + 2 <= pivot_counter[0] ? potential_pivot[largest_pivot+1][0] : 0;
            
            while(i <= pivot_counter[0]-1 && potential_pivot[i][1] != 1){
                distance_condition_right = distance_check(potential_pivot[largest_pivot][0],potential_pivot[i][0],t,store,check,user_tmin);
                slope_condition_right = slope_compare(potential_pivot[largest_pivot][0],test_slope_point,potential_pivot[largest_pivot][0],potential_pivot[i][0],t,y,store,check,top);
                
                if(distance_condition_right + slope_condition_right < 2){  potential_pivot[i][1] = 2; eligible_points--; }
                if(slope_condition_right == 1)  test_slope_point = potential_pivot[i][0];
                i++;
            }
        }
        
        q = 0;
        for(i = 0; i <= pivot_counter[0]-1; i++){
            if(potential_pivot[i][1] == 1)   pivots[p+(q++)][top] = potential_pivot[i][0]; 
        }        
    }
    else{
        if(pivot_counter[1] == 0){
            potential_pivot[0][1] = interval_start;
            (pivot_counter[1])++;
        }
        if(pivot_counter[2] == 0){
            potential_pivot[0][2] = interval_finish;
            (pivot_counter[2])++;
        }
        
        /* take preference for points that are within k*(interval width) of the respective end point */
        interval_width = 1*(t[interval_finish]-t[interval_start]);
        j = 0; k = pivot_counter[1]-1;
        if(pivot_counter[2] > 1 && distance_check(potential_pivot[0][2],interval_finish,t,store,check,interval_width)){
            for(i = 1; i <= pivot_counter[2]-1; i++){
                if(distance_check(potential_pivot[i][2],interval_finish,t,store,check,interval_width) == 0){
                    j = i; i = pivot_counter[2];
                }
            }
        }
        if(pivot_counter[1] > 1 && distance_check(interval_start,potential_pivot[pivot_counter[1]-1][1],t,store,check,interval_width)){
            for(i = pivot_counter[1]-2; i >= 0; i--){
                if(!distance_check(interval_start,potential_pivot[i][1],t,store,check,interval_width)){ k = i; i = -1; }
            }
        }
    }
    
    temp_a = potential_pivot[k][1], temp_b = potential_pivot[j][2];
    free(back_encroach);
    
    return q == 0 ? -((interval_finish-interval_start+1)*(temp_a-interval_start)+(temp_b-interval_start)) : q;
}

int phase_three(double *t, double *y, intRange *pivots, int *interval_list, int p, Range **store, intRange **check, double user_tmin, intRange *status, int iteration, int top){
    
    int q, left_encroach, right_encroach, pivot_counter = 0, p_temp, inflection_down, concave_start = 0, concave_finish, temp_pivot_count, k, alternate_counter = 0; 
    double area = 0;
    int *alternate; alternate = (int*) malloc(p*sizeof(int));
#ifdef db5
    printf("\nPHASE 3 begins\n");
#endif
    p_temp = p;
    for(q = 0; q <= p_temp-1; q++)   interval_list[q] = pivots[q][top];
#ifdef db4
    printf("PIVOTS:\n");
    for(q = 0; q <= p_temp-1; q++)  printf("%d ",pivots[q][top]);
    printf("\n");
#endif
    
    /* first, deal with the inflection intervals */
    for(q = 0; q <= p_temp-2; q++){
        
        /* if we are dealing with a 'down' inflection interval */
        if((q == p_temp-2 || slope_compare(interval_list[q],interval_list[q+1],interval_list[q+1],interval_list[q+2],t,y,store,check,top)) && (q == 0 || slope_compare(interval_list[q],interval_list[q+1],interval_list[q-1],interval_list[q],t,y,store,check,top))){
#ifdef db5
            printf("[%d %d] inflection down\t",interval_list[q],interval_list[q+1]);
#endif  
            pivots[pivot_counter++][top] = interval_list[q];  
            /* if q == p_temp-2 (last interval) set flag to 2, if q == 0 (first interval) set flag to 1 */
            find_inflection_points(t, y, user_tmin, interval_list, store, check, top, q, 1, &left_encroach, &right_encroach, q == p_temp-2 ? 2 : (q == 0 ? 1 : 0));
#ifdef db5
            printf("%d,%d\n",left_encroach,right_encroach);
#endif  
            
            /* if there is a point on the right of the interval, and the distance from the left suggested point is wide and the point fits the up slope condition */
            if(right_encroach != interval_list[q+1] && (distance_check(left_encroach,right_encroach,t,store,check,user_tmin) || (q == 0 && left_encroach == interval_list[q])) && (q == p_temp-2 || slope_compare(right_encroach,interval_list[q+1],right_encroach,interval_list[q+2],t,y,store,check,top))){
                
                inflection_down = q;
                /* if there is a point in the way that needs to be considered */
                if(interval_list[q] != left_encroach && slope_compare(interval_list[q],right_encroach,interval_list[q],left_encroach,t,y,store,check,top)){
                    pivot_counter += convex(t,y,interval_list[q],left_encroach+1,pivots,2,top,pivot_counter-1,store,check);
                }
                pivots[pivot_counter++][top] = right_encroach;
                alternate[alternate_counter++] = right_encroach;
            }
        }
        /* else if we are dealing with a 'up' inflection interval */
        else if((q == p_temp-2 || slope_compare(interval_list[q+1],interval_list[q+2],interval_list[q],interval_list[q+1],t,y,store,check,top)) && (q == 0 || slope_compare(interval_list[q-1],interval_list[q],interval_list[q],interval_list[q+1],t,y,store,check,top))){
                    
#ifdef db5        
            printf("[%d %d] inflection up\t",interval_list[q],interval_list[q+1]);
#endif
            pivots[pivot_counter++][top] = interval_list[q];
            /* find points for use by the inflection subroutine - we are looking for points that fail the slope condition on the left and are greater than user_tmin from the right */
            find_inflection_points(t, y, user_tmin, interval_list, store, check, top, q, 0, &left_encroach, &right_encroach, q == p_temp-2 ? 2 : (q == 0 ? 1 : 0));
#ifdef db5      
            printf("%d,%d\n",left_encroach,right_encroach);
#endif
                        
            /* if (0) we aren't too close to the start, and (1) the previous interval was a 'down' interval, and (2) we actually added a point in the down interval, and (3) the new value we wish to add is not valid, and (4) the horiztonal distance for the new interval addition is wider than the left-adjacent interval */
            if(q >= 2 && slope_compare(interval_list[q-1],interval_list[q],interval_list[q],interval_list[q+1],t,y,store,check,top) && slope_compare(interval_list[q-1],interval_list[q],interval_list[q-2],interval_list[q-1],t,y,store,check,top) && interval_list[q-1] != pivots[pivot_counter-2][top] && (left_encroach != interval_list[q] && !slope_compare(pivots[pivot_counter-2][top],interval_list[q],interval_list[q],left_encroach,t,y,store,check,top)) && distance_compare(pivots[pivot_counter-2][top],interval_list[q],interval_list[q],left_encroach,t,store,check)){
                /* go back and remove all the added points from the left-adjacent interval */                
                while(pivots[(pivot_counter--)-2][top] != interval_list[q-1]);
                pivots[pivot_counter++][top] = interval_list[q];
            }
            
            /* if there is a point on the left of the interval, and the distance from the right suggested point is wide and the point fits the up slope condition - ### some of these checks may not be needed */
            if(left_encroach != interval_list[q] && (distance_check(left_encroach,right_encroach,t,store,check,user_tmin) || (q == p_temp-2 && right_encroach == interval_list[q+1])) && (q == 0 || slope_compare(interval_list[q-1],interval_list[q],interval_list[q],left_encroach,t,y,store,check,top)) && (q == 0 || slope_compare(pivots[pivot_counter-2][top],interval_list[q],pivots[pivot_counter-2][top],left_encroach,t,y,store,check,top))){
                
                pivots[pivot_counter++][top] = left_encroach;
                alternate[alternate_counter++] = left_encroach;
                /* if there is a point in the way that needs to be considered */
                if(interval_list[q+1] != right_encroach && slope_compare(left_encroach,interval_list[q+1],left_encroach,right_encroach,t,y,store,check,top)){
                    pivot_counter+= convex(t,y,right_encroach,interval_list[q+1]+1,pivots,2,top,pivot_counter,store,check);
                }
            }
        }
        else{   /* else we are dealing with neither and no action needs to be taken */
            pivots[pivot_counter++][top] = interval_list[q];
        }
    }
    
    pivots[pivot_counter++][top] = interval_list[p_temp-1];
    
    p_temp = pivot_counter;
    for(q = 0; q <= p_temp-1; q++)   interval_list[q] = pivots[q][top];

#ifdef db4
    printf("PIVOTS:\n");
    for(q = 0; q <= p_temp-1; q++)  printf("%d ",pivots[q][top]);
    printf("\n");
#endif
    
    /* now deal with the concave regions */ /* ### what if we want to begin a concave region at the beginning */
    pivot_counter = 1; 
    for(q = 0; q <= p_temp-2; q++){
        /* if we hit a inflection down interval, begin the count */
        if((q >= p_temp-3 || slope_compare(interval_list[q],interval_list[q+1],interval_list[q+1],interval_list[q+2],t,y,store,check,top)) && (q == 0 || slope_compare(interval_list[q],interval_list[q+1],interval_list[q-1],interval_list[q],t,y,store,check,top))){
            pivots[pivot_counter++][top] = interval_list[q];
            concave_start = q+1;    
        }
        /* else we have hit an inflection up interval and we may begin the concave hull */
        else if(concave_start != 0 && (q >= p_temp-2 || slope_compare(interval_list[q+1],interval_list[q+2],interval_list[q],interval_list[q+1],t,y,store,check,top)) && slope_compare(interval_list[q-1],interval_list[q],interval_list[q],interval_list[q+1],t,y,store,check,top)){
            concave_finish = q;
                        
            /* first, calculate the area of the current configuration */
            area = 0;
#ifdef db6
            flop_counter += 3*(concave_finish-concave_start+1)+1;
#endif
            for(k = concave_start; k <= concave_finish-1; k++){
                area += t[interval_list[k]]*y[interval_list[k+1]] - y[interval_list[k]]*t[interval_list[k+1]];
            }
            area += t[interval_list[concave_finish]]*y[interval_list[concave_start]] - y[interval_list[concave_finish]]*t[interval_list[concave_start]];
            area = 0.5*abs_value(area);
            
#ifdef db5
            printf("Calling concave on [%d %d]\n",interval_list[concave_start],interval_list[concave_finish]);
#endif
            /* now call the concave subroutine */
            temp_pivot_count = interval_list[concave_start] != interval_list[concave_finish] ? concave(t, y, interval_list[concave_start], interval_list[concave_finish], interval_list[concave_start-1], interval_list[concave_finish+1], interval_list[concave_start+1], interval_list[concave_finish-1], pivots, store, check, alternate, alternate_counter, area, user_tmin, 2, top, pivot_counter) : -1;
            
            if(temp_pivot_count >= 0)   pivot_counter += temp_pivot_count;
            else{
                for(k = concave_start; k <= concave_finish; k++)   pivots[pivot_counter++][top] = interval_list[k];
            }
            concave_start = 0;
        }
        else if(concave_start == 0){
            pivots[pivot_counter++][top] = interval_list[q];    
        }
    }
    
    
    if(pivots[pivot_counter-1][top] != interval_list[p_temp-2]){
        pivots[pivot_counter++][top] = interval_list[p_temp-2];
    }
    pivots[pivot_counter++][top] = interval_list[p_temp-1];
    
#ifdef db4
    printf("PIVOTS:\n");
    for(q = 0; q <= pivot_counter-1; q++)  printf("%d ",pivots[q][top]);
    printf("\n");
#endif
    
    free(alternate);
    return pivot_counter;
}

int find_inflection_points(double *t, double *y, double user_tmin, int *interval_list, Range **store, intRange **check, int top, int k, int up_down, int *left_encroach, int *right_encroach, int flag){
    
    intRange* inflection_points; 
    inflection_points = (intRange*) calloc(interval_list[k+1]-interval_list[k]+1,sizeof(intRange));
    int inflection_counter = 1, q, i, j, stop_adding = 0, temp_slope_point;
    
    /* dealing with an inflection-up interval */
    if(up_down == 0){
        
        inflection_points[0][0] = interval_list[k];
        for(q = interval_list[k]+1; q <= interval_list[k+1]-1; q++){
            for(i = 0; i <= inflection_counter-1; i++){
                if((inflection_points[i][1] == 0 || !slope_compare(inflection_points[i][0],q,inflection_points[i][0],inflection_points[i][1],t,y,store,check,top)) && !slope_compare(inflection_points[i][0],q,inflection_points[i][0],interval_list[k+1],t,y,store,check,top)){
                    
                    
                    if(distance_check(inflection_points[i][0],q,t,store,check,user_tmin)){
                        inflection_points[i][1] = q; 
                    }
                    else{
                        for(j = q+1; j <= interval_list[k+1]-1; j++){
                            if(slope_compare(inflection_points[i][0],q,inflection_points[i][0],j,t,y,store,check,top)){
                                j = interval_list[k+1]+2;
                            }
                        }
                        
                        if(j != interval_list[k+1]+3){
                            for(j = i; j <= inflection_counter-1; j++){
                                inflection_points[j][0] = inflection_points[j+1][0];
                                inflection_points[j][1] = inflection_points[j+1][1];
                            }
                            i--; inflection_counter--;
                        }
                    }
                }
            }

            if(stop_adding == 0 && flag != 2 && !distance_check(q,interval_list[k+1],t,store,check,user_tmin))  stop_adding = 1;
                        
            if(stop_adding == 0 && ((inflection_counter == 1 && (flag == 1 || k == 0 || slope_compare(interval_list[k-1],interval_list[k],interval_list[k],q,t,y,store,check,top))) || (inflection_counter > 1 && slope_compare(interval_list[k],temp_slope_point,interval_list[k],q,t,y,store,check,top)))){
                inflection_points[inflection_counter++][0] = q;
                temp_slope_point = q;
            }
        }
        
        *left_encroach = inflection_points[inflection_counter-1][0];
        *right_encroach = inflection_points[inflection_counter-1][1] == 0 ? interval_list[k+1] : inflection_points[inflection_counter-1][1];
    }
    else{
        inflection_points[0][0] = interval_list[k+1];
        for(q = interval_list[k+1]-1; q >= interval_list[k]+1; q--){
            for(i = 0; i <= inflection_counter-1; i++){
                if((inflection_points[i][1] == 0 || !slope_compare(inflection_points[i][1],inflection_points[i][0],q,inflection_points[i][0],t,y,store,check,top)) && !slope_compare(interval_list[k],q,interval_list[k],inflection_points[i][0],t,y,store,check,top)){
                                        
                    if(distance_check(q,inflection_points[i][0],t,store,check,user_tmin)){
                        inflection_points[i][1] = q;
                    }
                    else{
                        for(j = q-1; j >= interval_list[k]+1; j--){
                            if(slope_compare(j,q,j,inflection_points[i][0],t,y,store,check,top)){
                                j = -2;
                            }
                        }
                        
                        if(j != -3){
                            for(j = i; j <= inflection_counter-1; j++){
                                inflection_points[j][0] = inflection_points[j+1][0];
                                inflection_points[j][1] = inflection_points[j+1][1];
                            }
                            i--; inflection_counter--;
                        }
                        
                    }
                }
            }
            
            if(stop_adding == 0 && flag != 1 && !distance_check(interval_list[k],q,t,store,check,user_tmin)) stop_adding = 1;
            if(stop_adding == 0 && ((inflection_counter == 1 && (flag == 2 || slope_compare(q,interval_list[k+1],interval_list[k+1],interval_list[k+2],t,y,store,check,top))) || (inflection_counter > 1 && slope_compare(q,interval_list[k+1],temp_slope_point,interval_list[k+1],t,y,store,check,top)))){
                inflection_points[inflection_counter++][0] = q;
                
                temp_slope_point = q;            
            }
        }
        
        *right_encroach = inflection_points[inflection_counter-1][0];
        *left_encroach = inflection_points[inflection_counter-1][1] == 0 ? interval_list[k] : inflection_points[inflection_counter-1][1];         
    }
    
    free(inflection_points);           
    return 1;
}

int concave(double *t, double *y, int start, int finish, int before, int next, int start_alt, int finish_alt, intRange *pivots, Range **store, intRange **check, int *alternate, int alternate_counter, double input_area, double user_tmin, int initial, int top, int insert){
    
    intRange3 *hull; hull = (intRange3*) malloc((finish-start+3)*sizeof(intRange3));
    int i, q, q_back;
    int k[3]; k[0] = 0; k[1] = 0; k[2] = 0; /* number of elements in the hull k[0]=forward, k[1]=backward */
    double area[3]; area[0] = 0; area[1] = 0; area[2] = 0; /* area removed by each hull area[0]=forward, area[1]=backward */
    
#ifdef db5
    printf("\tarea = %f (original)\n",input_area);
#endif
    
    /* forward direction concave banding */
    hull[0][0] = before;
    hull[1][0] = start;
    k[0] = 2;
    for(i = start+1; i <= finish; i++){ 
        /* if the slope is greater than the slope of the segment connecing the previous two points and does not cause a violation with the right-adjacent interval, then this point may be added to the list */
        if(slope_compare(hull[k[0]-2][0],hull[k[0]-1][0],hull[k[0]-2][0],i,t,y,store,check,top) && slope_compare(i,next,finish,next,t,y,store,check,top))   hull[k[0]++][0] = i; 
        
        /* if the newly added point violates the slope condition */
        while(k[0] >= 4 && !slope_compare(hull[k[0]-2][0],hull[k[0]-1][0],hull[k[0]-2][0],finish,t,y,store,check,top)){
            hull[k[0]-2][0] = hull[k[0]-1][0]; 
            k[0]--;
        }
    }  
    
    /* calculate the area of the region using the surveyor's formula */
    hull[k[0]][0] = hull[1][0]; /* for the wrap-around property of the area formula */
    for(q = 1; q <= k[0]-1; q++)   area[0] += (t[hull[q][0]]*y[hull[q+1][0]]-y[hull[q][0]]*t[hull[q+1][0]]);
    area[0] = 0.5*abs_value(area[0]);
#ifdef db5
    printf("\tarea = %f (forward)\t",area[0]);
    for(q = 1; q <= k[0]-1; q++)   printf("%d ",hull[q][0]);
    printf("\n");
#endif
    
    /* backward direction concave banding - same structure as above */
    hull[0][1] = next;
    hull[1][1] = finish;
    k[1] = 2;
    for(i = finish-1; i >= start; i--){
        if(slope_compare(i,hull[k[1]-2][1],hull[k[1]-1][1],hull[k[1]-2][1],t,y,store,check,top) && slope_compare(before,start,before,i,t,y,store,check,top))   hull[k[1]++][1] = i; 
        
        while(k[1] >= 4 && !slope_compare(start,hull[k[1]-2][1],hull[k[1]-1][1],hull[k[1]-2][1],t,y,store,check,top)){
            hull[k[1]-2][1] = hull[k[1]-1][1]; 
            k[1]--;
        }
    }  
    
    /* calculate area of backward direction */
    hull[k[1]][1] = hull[1][1]; /* for the wrap-around property of the area formula */
    for(q = 1; q <= k[1]-1; q++)   area[1] += (t[hull[q][1]]*y[hull[q+1][1]]-y[hull[q][1]]*t[hull[q+1][1]]);
    area[1] = 0.5*abs_value(area[1]);
#ifdef db5
    printf("\tarea = %f (backward)\t",area[1]);
    for(q = k[1]-1; q >= 1; q--)   printf("%d ",hull[q][1]);
    printf("\n");
#endif
    
    /* alternate concave banding if there was an inflection point added */
    for(q = 0; q <= alternate_counter-1; q++){
        if(alternate[q] == start || alternate[q] == finish){
            hull[0][2] = before;
            hull[1][2] = alternate[q] == start ? start_alt : start;
            k[2] = 2;
            for(i = hull[1][2]+1; i <= (alternate[q] == finish ? finish_alt : finish); i++){ 
                /* if the slope is greater than the slope of the segment connecing the previous two points and does not cause a violation with the right-adjacent interval, then this point may be added to the list */
                if(slope_compare(hull[k[2]-2][2],hull[k[2]-1][2],hull[k[2]-2][2],i,t,y,store,check,top) && slope_compare(i,next,finish,next,t,y,store,check,top))   hull[k[2]++][2] = i; 
                
                /* if the newly added point violates the slope condition */
                while(k[2] >= 4 && !slope_compare(hull[k[2]-2][2],hull[k[2]-1][2],hull[k[2]-2][2],finish,t,y,store,check,top)){
                    hull[k[2]-2][2] = hull[k[2]-1][2]; 
                    k[2]--;
                }
            }  
            /* calculate the area of the region using the surveyor's formula */
            hull[k[2]][2] = hull[1][2]; /* for the wrap-around property of the area formula */
            for(q = 1; q <= k[2]-1; q++)   area[2] += (t[hull[q][2]]*y[hull[q+1][2]]-y[hull[q][2]]*t[hull[q+1][2]]);
            area[2] = 0.5*abs_value(area[2]);
#ifdef db5
            printf("\tarea = %f (alternate)\t",area[2]);
            for(q = 1; q <= k[2]-1; q++)   printf("%d ",hull[q][2]);
            printf("\n");
#endif
            q = alternate_counter+1;
        }
    }
    

    
    /* depending on which configuration is the best, we place hull into pivot list */
    if(area[1] > area[0] && area[1] > input_area && area[1] > area[2]){ 
        for(q = k[1]-1; q >= 1; q--)    pivots[insert+k[1]-1-q][top] = hull[q][1];
        i = k[1]-1;
    }
    else if(area[0] >= area[1] && area[0] > input_area && area[0] > area[2]){
        for(q = 1; q <= k[0]-1; q++)    pivots[insert+q-1][top] = hull[q][0];
        i = k[0]-1;
    }
    else if(pivots[insert-2][top]+pivots[insert-1][top]!=0 && area[2]>=area[0] && area[2]>=area[1] && area[2]>input_area){
        q = -1; q_back = 0;
        while(!distance_check(pivots[insert+q-1][top],pivots[insert+q][top],t,store,check,user_tmin) && slope_compare(pivots[insert+q-1][top],pivots[insert+q][top],pivots[insert+q][top],hull[1][2],t,y,store,check,top)){
            q--; q_back++;
        }
        for(q = 1; q <= k[2]-1; q++)    pivots[insert+q-1-q_back][top] = hull[q][2];
        i = k[2]-1-q_back;
    }
    else{
        i = -1;
    }
    
    free(hull);
    return i;
}
