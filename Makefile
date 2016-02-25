# simulation options

#OPT   +=  -DADAPTATIVE_TIMESTEP
#OPT   +=  -DBULIRSCH_STOER
#OPT   +=  -DEXTENDED_BODY
#OPT   +=  -DSOFTENING


CC=gcc
UBICATIONGSL=/home/lquiroga/local
#CFLAGS=-g -I.                                                                  
CFLAGS= -I. -O4 -Wall -I$(UBICATIONGSL)/include $(OPT)
LFLAGS= -lm -L$(UBICATIONGSL)/lib -lgsl -lgslcblas

%.out:%.o
	$(CC) $^ $(LFLAGS) -o $@
#       cp $@ /home/lquiroga/local/bin

clean:
	rm -rf *.o* *~ *.out

#clean_ajuste:
#        rm parametros_* puntos_*
