TMPDIR = /tmp/phast
export CVSROOT = /home/cvsroot
CWD = ${PWD}

all:
	@echo "Type \"make package\" to create a tarball reflecting the current state of the CVS tree."

package:
	rm -rf ${TMPDIR}
	mkdir -p ${TMPDIR}
	cd ${TMPDIR} ; cvs checkout phast 
	rm -r `find ${TMPDIR}/phast -name "CVS"`
	VERSION = `cat ${TMPDIR}/phast/version | sed 's/\./_/g'`
	cd ${TMPDIR} ; tar cfz ${CWD}/phast.$$VERSION.tgz phast
