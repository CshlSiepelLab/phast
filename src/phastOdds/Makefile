include ../make-include.mk
PHAST := ${PHAST}/..
# assume executable name is given by directory name
EXEC = $(addprefix ${BIN}/,$(notdir ${PWD}))

# assume all *.c files are source
SRCS = $(basename $(wildcard *.c))
OBJS =  $(addsuffix .o,${SRCS})
HELP = $(addsuffix .help,$(basename $(wildcard *.help_src)))

%.o : %.c
# (cancels built-in rule; otherwise gets used instead if *.help missing)
.SECONDARY : ${HELP}
# (prevents *.help from being deleted as a intermediate file)

all: ${EXEC}

%.o: %.c ${HELP} ../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${EXEC}: ${OBJS} ${PHAST}/lib/libphast.a
	${CC} ${LFLAGS} ${LIBPATH} -o $@ ${OBJS} ${LIBS} 

%.help : %.help_src
	../munge-help.sh $< > $@

clean: 
	rm -f *.o *.help ${EXEC} 
