include ../../make-include.mk

# assume every *.c file goes in the library
OBJS =  $(addsuffix .o,$(basename $(wildcard *.c)))

all: ${TARGETLIB} 

%.o: %.c ../../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${TARGETLIB} : ${OBJS} 
	${AR} rs ${TARGETLIB} $?

clean:
	rm -f *.o 
