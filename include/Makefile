# Generate man pages for individual functions from .h headers.  
# Uses 'c2man'.

C2MAN = c2man			# must be in path

MAN_DEST = ../man/man3
C2MAN_FLAGS = -g -F"int f(a, b)" -n -L -o${MAN_DEST} -l n
# TODO: maybe don't want to use -g for all 

TMPFILE = /tmp/c2man.tmp
ADDENDUMFILE = c2man-addendum

#targets = $(wildcard *.h)
NAMEROOTS = lists stacks queues stringsplus gff category_map hashtable matrix
SRC = $(addsuffix .h,${NAMEROOTS})
TARGETS = $(addprefix ${MAN_DEST}/,$(addsuffix .3,${NAMEROOTS}))

all: makedocs postprocess

makedocs: 
	for file in ${SRC} ; do ${C2MAN} ${C2MAN_FLAGS} $$file; done

postprocess:
	for file in `find ${MAN_DEST} -type f -name "*.3"` ; do cat $$file ${ADDENDUMFILE} > ${TMPFILE} && mv ${TMPFILE} $$file ; done

clean:
	rm -f ${MAN_DEST}/*.3
