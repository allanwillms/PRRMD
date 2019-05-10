#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "read_data.h"
#define MAXLINE 5000

double ** read_data(char *fname,int *ncols, int *nrows) {
    /* This function reads in columns of data from fname.  It outputs the number of
     * rows and columns in fname, and returns a (double **) so that if
     * x = read_data(fname,&ncols,&nrows), then
     * x[j][i] references the i'th row and j'th column of the data file. */

    double **data;
    register int i,j,prev_was_sp;
    FILE *fptr;
    char line[MAXLINE+1];
    char c;

    if ((fptr = fopen(fname,"r")) == NULL) {
	fprintf(stderr,"read_data: error opening %s.\n",fname);
	exit(1);
    }
    /* Count the columns in the first line. */
    fgets(line,MAXLINE,fptr);
    i = 0;
    j = 0;
    prev_was_sp = 1;
    while ((c = *(line+i++)) != '\n') {
	if (isspace(c))
	    prev_was_sp = 1;
	else {
	    if (prev_was_sp) {
		prev_was_sp = 0;
		j++;
	    }
	}
    }
    *ncols = j;
    /* Count the lines. */
    i = 1;  /* Already read one line. */
    while (fgets(line,MAXLINE,fptr) != NULL) i++;
    *nrows = i;
    /* Allocate space for the data. */
    data = (double **) malloc(*ncols*sizeof(double *));
    for (i=0; i<*ncols; i++) data[i] = (double *) malloc(*nrows*sizeof(double));
    /* Rewind to start and read in the data. */
    rewind(fptr);
    for (i=0; i<*nrows; i++) {
	for (j=0; j<*ncols; j++) {
	    if (fscanf(fptr,"%lf",data[j]+i) != 1) {
		fprintf(stderr,"read_data: error reading input file %s "
			"at line %d column %d.\n",fname,i+1,j+1);
		exit(1);
	    }
	}
    }
    if (fclose(fptr) != 0) {
	fprintf(stderr,"read_data: error closing input file %s.\n",fname);
	exit(1);
    }
    return data;
}
