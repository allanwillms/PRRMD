OBJS = 	parm_reduction.o  A1OUT.o NLP.o NLP_u.o PHARM.o \
	interval_math.o readtoken.o read_data.o deblank.o tasle.o 

CFLAGS = -Wall 
LDFLAGS = -o parm_reduction -lm
CC = gcc

parm_reduction : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

parm_reduction.o : parm_reduction.h parm_red_structs.h\
                   parm_red_discr_sets.h parm_red_models.h \
		   readtoken.h read_data.h tasle.h interval_math.h

interval_math.o : interval_math.h parm_red_structs.h

A1OUT.o : parm_red_structs.h readtoken.h

PHARM.o : parm_red_structs.h readtoken.h

NLP.o : parm_red_structs.h interval_math.h readtoken.h

NLP_u.o : parm_red_structs.h interval_math.h readtoken.h

readtoken.o : readtoken.h

read_data.o : read_data.h

deblank.o : deblank.h

tasle.o : tasle.h

discr_error_test : parm_red_structs.h parm_reduction.h discr_error_test.h

.PHONY : clean
clean :
	-rm $(OBJS)
