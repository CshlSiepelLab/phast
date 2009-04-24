include ../make-include.mk

#CFLAGS += -g -fno-inline -Wall

#EXEC = $(addprefix ${BIN}/,$(notdir ${PWD}))
PROGS = dmotif dmsample dmsimulate dmcompare dmcondition
MODULES = dmotif_phmm pssm dmotif_indel_mod
EXEC = $(addprefix ${BIN}/,${PROGS})

# assume all *.c files are source
OBJS =  $(addsuffix .o,${PROGS} ${MODULES})
OBJS_MODULES =  $(addsuffix .o,${MODULES})
HELP = $(addsuffix .help,$(basename $(wildcard *.help_src)))

%.o : %.c
# (cancels built-in rule; otherwise gets used instead if *.help missing)
.SECONDARY : ${HELP}
# (prevents *.help from being deleted as an intermediate file)

all: ${EXEC}

%.o: %.c ${HELP} ../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${EXEC}: ${BIN}/%: %.o ${OBJS_MODULES} ${PHAST}/lib/libphast.a
	${CC} ${LFLAGS} ${LIBPATH} -o $@ $< ${OBJS_MODULES} ${LIBS} ${CFLAGS}

%.help : %.help_src
	../munge-help.sh $< > $@

clean: 
	rm -f *.o ${HELP} ${EXEC} 
