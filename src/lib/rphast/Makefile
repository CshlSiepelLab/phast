include ../../make-include.mk

# assume every *.c file goes in the library
OBJS =  $(addsuffix .o,$(basename $(wildcard *.c)))

all: ${LIB}/librphast.so

%.o: %.c ../../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${LIB}/librphast.so: ${OBJS}
	R CMD SHLIB ${OBJS} ${TARGETLIB} -o $@

clean:
	rm -f *.o ${LIB}/librphast.so
