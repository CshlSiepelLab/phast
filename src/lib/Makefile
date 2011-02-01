include ../make-include.mk
PHAST := ${PHAST}/..

CDIR = ${PWD}

PACKAGES=base msa feature hmm phylo phylo_hmm motif pcre

ifdef RPHAST
PACKAGES += rphast
endif


ifdef CLAPACKPATH		# only do below if using CLAPACK

all: liblinks phastlib

sharedlib: liblinks sharedphastlib

# by default, CLAPACK uses non-standard library names, which confuses
# the compiler.  We'll just create symbolic links in ${PHAST}/lib that
# have proper names
liblinks: ${PHAST}/lib/liblapack.a ${PHAST}/lib/libtmg.a \
	${PHAST}/lib/libblaswr.a 

${PHAST}/lib/liblapack.a: ${CLAPACKPATH}/lapack${PLAT}.a
	mkdir -p ${PHAST}/lib
	${LN} -s $< $@

${PHAST}/lib/libtmg.a: ${CLAPACKPATH}/tmglib${PLAT}.a
	mkdir -p ${PHAST}/lib
	${LN} -s $< $@

${PHAST}/lib/libblaswr.a: ${CLAPACKPATH}/blas${PLAT}.a
	mkdir -p ${PHAST}/lib
	${LN} -s $< $@

else				# no CLAPACK: links not necessary

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
