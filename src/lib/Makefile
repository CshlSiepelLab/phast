include ../make-include.mk
PHAST := ${PHAST}/..

CDIR = ${PWD}

PACKAGES=base msa feature hmm phylo phylo_hmm motif pcre

ifdef RPHAST
PACKAGES += rphast
endif


ifdef CLAPACKPATH		# only do below if using CLAPACK

all: phastlib

sharedlib: sharedphastlib

all: phastlib

sharedlib: sharedphastlib

endif

phastlib:
	for dir in ${PACKAGES}; do cd ${CDIR}/$$dir && ${MAKE} ; done

sharedphastlib:
	for dir in ${PACKAGES}; do cd ${CDIR}/$$dir && ${MAKE} sharedlib ; done

clean:
	for dir in ${PACKAGES}; do cd ${CDIR}/$$dir && ${MAKE} clean ; done
	rm -f ${PHAST}/lib/*
