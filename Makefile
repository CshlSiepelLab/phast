TMPDIR = /tmp/phast
CWD = ${PWD}

all:
	@echo "Type \"make package\" to create a tarball reflecting the current state of the CVS tree."

package:
	rm -rf ${TMPDIR}
	mkdir -p ${TMPDIR}
	cd ${TMPDIR} ; cvs checkout phast 
	find ${TMPDIR}/phast -name "CVS" | xargs rm -rf
	rm -r ${TMPDIR}/phast/doc
	VERSION=`cat ${TMPDIR}/phast/version | sed 's/\./_/g'` ;\
	cd ${TMPDIR} ; tar cfz ${CWD}/phast.$$VERSION.tgz --exclude test phast
	rm -rf ${TMPDIR}
