include ../make-include.mk

# assume every *.c file represents an executable
EXECS = $(addprefix ${BIN}/,$(basename $(wildcard *.c)))
OBJS =  $(addsuffix .o,${EXECS})

all: ${EXECS}

%.o: %.c ../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${EXECS}: ${BIN}/%: %.o ${PHAST}/lib/libphast.a
	${CC} ${LIBPATH} ${LFLAGS} -o $@ $< ${LIBS} 

clean: 
	rm -f *.o ${EXECS} 
