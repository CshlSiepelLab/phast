include ../make-include.mk
PHAST := ${PHAST}/..
# assume executable name is given by directory name
EXEC = $(addprefix ${BIN}/,$(notdir ${PWD}))

# assume all *.c files are source
SRCS = $(basename $(wildcard *.c))
OBJS =  $(addsuffix .o,${SRCS})

all: ${EXEC}

%.o: %.c ../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${EXEC}: ${OBJS} ${PHAST}/lib/libphast.a
	${CC} ${LFLAGS} ${LIBPATH} -o $@ ${OBJS} ${LIBS} 

clean: 
	rm -f *.o ${EXEC} 
