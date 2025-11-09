include ../../make-include.mk
PHAST := ${PHAST}/../..

# assume every *.c file goes in the library
OBJS =  $(addsuffix .o,$(basename $(wildcard *.c)))

all: ${TARGETLIB}

sharedlib: CFLAGS += -fPIC
sharedlib: all

%.o: %.c ../../make-include.mk
	$(CC) $(CFLAGS) -c $< -o $@ 

$(TARGETLIB): $(OBJS)
	$(AR) rcs $@ $^

clean:
	rm -f *.o 
