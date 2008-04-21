include ../../make-include.mk

# assume every *.c file goes in the library
OBJS =  $(addsuffix .o,$(basename $(wildcard *.c)))

ALLOBJS = $(wildcard ../*/*.o)

all: ${LIB}/librphast.so

%.o: %.c ../../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

${LIB}/librphast.so: ${OBJS}
#	R CMD SHLIB ${OBJS} ${TARGETLIB} -o $@
	R CMD SHLIB ${ALLOBJS} -o $@

# FIXME: sometimes doesn't work to link to libphast.a.  Seems to be
# version dependent (either R or gcc, not sure which)

clean:
	rm -f *.o ${LIB}/librphast.so
