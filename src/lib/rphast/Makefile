include ../../make-include.mk

# assume every *.c file goes in the library
OBJS =  $(addsuffix .o,$(basename $(wildcard *.c)))

ALLOBJS = $(wildcard ../*/*.o)

all: ${LIB}/librphast.so

%.o: %.c ../../make-include.mk 
	$(CC) $(CFLAGS) -c $< -o $@ 

${LIB}/librphast.so: ${OBJS} ${ALLOBJS}
	R CMD SHLIB -fPIC ${LIBPATH} ${LIBS} ../*/*.o -o $@

clean:
	rm -f *.o ${LIB}/librphast.so
