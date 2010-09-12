include ../make-include.mk
PHAST := ${PHAST}/..

# assume executable name is given by directory name
PROGS = pbsDecode pbsEncode pbsScoreMatrix pbsTrain prequel
MODULES = simplex_grid pbs_code
EXEC = $(addprefix ${BIN}/,${PROGS})

# assume all *.c files are source
OBJS =  $(addsuffix .o,${PROGS} ${MODULES})
OBJS_MODULES =  $(addsuffix .o,${MODULES})
HELP = $(addsuffix .help,$(basename $(wildcard *.help_src)))

%.o : %.c
# (cancels built-in rule; otherwise gets used instead if *.help missing)
.SECONDARY : ${HELP}
# (prevents *.help from being deleted as a intermediate file)

all: ${EXEC}

%.o : %.c ${HELP} ../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${EXEC} : ${BIN}/%: %.o ${PHAST}/lib/libphast.a ${OBJS_MODULES}
	${CC} ${LFLAGS} ${LIBPATH} -o $@ $< ${OBJS_MODULES} ${LIBS} 

%.help : %.help_src
	../munge-help.sh $< > $@

clean: 
	rm -f *.o ${EXEC} ${HELP}
