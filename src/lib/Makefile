include ../make-include.mk

CDIR = ${PWD}

PACKAGES=base msa feature hmm phylo phylo_hmm motif

ifndef SKIP_CLAPACK		# only do below if using CLAPACK

all: liblinks phastlib

# CLAPACK uses non-standard library names, which confuses the
# compiler.  We'll just create symbolic links in ${PHAST}/lib that have
# proper names
liblinks: ${PHAST}/lib/liblapack.a ${PHAST}/lib/libtmg.a \
	${PHAST}/lib/libblaswr.a 

${PHAST}/lib/liblapack.a: ${CLAPACKPATH}/lapack_${ARCH}.a
	mkdir -p ${PHAST}/lib
	${LN} -s $< $@

${PHAST}/lib/libtmg.a: ${CLAPACKPATH}/tmglib_${ARCH}.a
	mkdir -p ${PHAST}/lib
	${LN} -s $< $@

${PHAST}/lib/libblaswr.a: ${CLAPACKPATH}/blas_${ARCH}.a
	mkdir -p ${PHAST}/lib
	${LN} -s $< $@

else				# SKIP_CLAPACK case: links not necessary

all: phastlib

endif

phastlib:
	for dir in ${PACKAGES}; do cd ${CDIR}/$$dir && ${MAKE} ; done

clean:
	for dir in ${PACKAGES}; do cd ${CDIR}/$$dir && ${MAKE} clean ; done
	rm -f ${PHAST}/lib/*