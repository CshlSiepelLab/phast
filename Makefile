TMPDIR = /tmp/phast
#export CVSROOT = /projects/compbio/cvsroot
CWD = ${PWD}

all:
	@echo "Type \"make package\" to create a tarball reflecting the current state of the CVS tree."

package:
	rm -rf ${TMPDIR}
	mkdir -p ${TMPDIR}
	cd ${TMPDIR} ; cvs checkout phast 
	find ${TMPDIR}/phast -name "CVS" | xargs rm -rf
	VERSION=`cat ${TMPDIR}/phast/version | sed 's/\./_/g'` ;\
	cd ${TMPDIR} ; tar cfz ${CWD}/phast.$$VERSION.tgz phast
	rm -rf ${TMPDIR}