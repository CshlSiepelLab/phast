include ../make-include.mk
PHAST := ${PHAST}/..

# assume every *.c file represents an executable
EXECS = $(addprefix ${BIN}/,$(basename $(wildcard *.c)))
OBJS =  $(addsuffix .o,${EXECS})
HELP = $(addsuffix .help,$(basename $(wildcard *.help_src)))

%.o : %.c
# (cancels built-in rule; otherwise gets used instead if *.help missing)
.SECONDARY : ${HELP}
# (prevents *.help from being deleted as a intermediate file)

all: ${EXECS}

%.o: %.c ${HELP} ../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${EXECS}: ${BIN}/%: %.o ${PHAST}/lib/libphast.a
	${CC} ${LIBPATH} ${LFLAGS} -o $@ $< ${LIBS} 

%.help : %.help_src
	../munge-help.sh $< > $@

clean: 
	rm -f *.o *.help ${EXECS} 
