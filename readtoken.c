#define MAX_LINE 120
#include <stdio.h>
#include <string.h>
#include "readtoken.h"
#include "deblank.h"

/*----------------------------------------------------------------------*/
/* READTOKEN: Read a line of token [sep] value from file
 where [sep] is either = or : or ;
 INPUTS 
 fin: file pointer
 num: maximum number of values to read
 str: string token to look for
 type: data type, any valid character conversion string 
 including the leading %. eg.
 %s - string
 %5s - string with at most five characters (excluding NUL)
 %c - single character
 %5c - five characters
 %d or %i - decimal integer
 %ld - long decimal integer
 %u - unsigned decimal integer
 %o - unsigned octal integer
 %hx - unsigned short hexidecimal integer
 %f or %e or %g - real float
 %lf - double precision float
 ptr: pointer to memory to store read tokens
 OUTPUTS
 sep: separator character used
 RETURNS
 number of values read
 
 If the input num is less than or equal to zero, then no reading is done and 
 the value 0 is returned.  If num>0, the number of values read is returned 
 unless the token itself is not found in which case -1 is returned.
 
 STRINGS: (type="%s", "%9s" etc.)
 If only one string is requested (num=1) then ptr should be a pointer to
 a character.
 If you wish to request more than one string then ptr should be (a pointer 
 to) an array of pointers to characters.  Thus, for example, if the calling 
 function wishes to read in four strings, then it should have a memory 
 location defined as
 char *in_strings[4];
 for (i=0; i<4; i++) in_strings[i] = (char *) malloc(X*sizeof(char));
 and should call this function with num=4 and ptr=in_strings.  X should be 
 large enough to hold the size of each string plus the NUL character.  To 
 prevent reading in strings that are too long, type could be specified as
 type="%Ys"  where Y=X-1.
 Note that declaring in_strings as 
 char in_strings[4][X];
 is not the same thing.
 
 CHARACTERS: (type="%Xc")
 If type is specified as "%Xc" where X>1 then ptr should be a pointer to a 
 character and the read values will be put every X locations after pointer.
 
 This implementation does not support long longs.
 */
int readtoken(FILE *fin, int num, char *str, char *type,
              void *ptr, char *separator) {
    
    /* local variables */
    char line[MAX_LINE], *tok, *rs, *val, *sep;
    int i,read,ptrfac;
    
    if (num <= 0) return 0;
    
    if (strstr(type,"c") != NULL) {
        if (sscanf(type,"%%%d",&ptrfac) !=1) ptrfac = 1;
    }
    
    read = -1;
    while (read<num) {
        if (fgets(line,MAX_LINE,fin) == NULL) break;
        if ((sep = strpbrk(line,"=:;")) != NULL) {     /* has a separator */
            *separator = *sep;                        /* keep it */
            if ((tok = strtok(line,"=:;")) != NULL && /* not an empty line */
                (*deblank(tok) != '%') )  {         /* not a comment line */
                /* Read the token value if the token matches specified one */
                if (strcmp(tok,str)==0) {
                    read = 0;
                    rs = sep+1;
                    for (i=0; i<num; i++) {        /* Read up to num values */
                        /* Read next value into string */
                        if ((val = strtok(rs," ,\t")) == NULL) break;  
                        if (strstr(type,"s") != NULL) {
                            if (num==1) {
                                if (sscanf(val,type,(char *)ptr) !=1) break;
                                else read += 1;
                            }
                            else {
                                if (sscanf(val,type,*((char **)ptr+i)) !=1) break;
                                else read += 1;
                            }
                        }
                        else if (strstr(type,"c") != NULL) {
                            if (sscanf(val,type,(char *)ptr+ptrfac*i) !=1) break;
                            else read += 1;
                        }
                        else if (strstr(type,"d") != NULL || 
                                 strstr(type,"i") != NULL ||
                                 strstr(type,"u") != NULL ||
                                 strstr(type,"o") != NULL ||
                                 strstr(type,"x") != NULL) {
                            if (strstr(type,"h") != NULL) {
                                if (sscanf(val,type,(short *)ptr+i) !=1) break;
                                else read += 1;
                            }
                            else if (strstr(type,"l") != NULL) {
                                if (sscanf(val,type,(long *)ptr+i) !=1) break;
                                else read += 1;
                            }
                            else {
                                if (sscanf(val,type,(int *)ptr+i) !=1) break;
                                else read += 1;
                            }
                        }
                        else if (strstr(type,"f") != NULL || 
                                 strstr(type,"e") != NULL ||
                                 strstr(type,"g") != NULL) {
                            if (strstr(type,"l") != NULL) {
                                if (sscanf(val,type,(double *)ptr+i) !=1)
                                    break;
                                else read += 1;
                            }
                            else {
                                if (sscanf(val,type,(float *)ptr+i) !=1) break;
                                else read += 1;
                            }
                        }
                        rs = NULL;    /* So next call to strtok will work */
                    }
                    break;
                }
            }
        }
    }
    return read;
}
