#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "deblank.h"

/*----------------------------------------------------------------------*/
/* DEBLANK:  Remove white-space characters from the start and end of
             a string.
     inputs:  pointer to null-terminated string
     outputs: modifies the input string
     returns: value of input pointer

     Trailing and leading white space is removed from the string
     and the resulting string is shifted so that the first non white
     space character is now at the start of the string.
     If a null pointer is passed in, nothing happens.
*/
char *deblank(char *str) {

   /* local variables */
   int i;
   char *cptr;

   if (str == NULL) return str;
   /* Find the last non white space character and terminate the
      string there. */
   cptr = strchr(str,'\0');
   while (--cptr-str>=0 && isspace(*cptr)) ;
   *++cptr = '\0';
   /* Find the first non white space character.  The null character is
      not white space so we need not worry about going past the end. */
   cptr = str-1;
   while (isspace(*++cptr)) ;
   i = 0;
   while ((*(str+i++) = *cptr++) != '\0') ;

   return str;
}

