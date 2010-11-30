TMPDIR = /tmp/phast
CWD = ${PWD}

all:
	@echo "Type \"make package\" to create a tarball reflecting the current state of the CVS tree."
	(cd src; make DESTDIR=${DESTDIR} CLAPACKPATH=/usr/lib )

package:
	rm -rf ${TMPDIR}
	mkdir -p ${TMPDIR}
	cd ${TMPDIR} ; svn checkout http://compgen.bscb.cornell.edu/svnrepo/phast/trunk phast
	find ${TMPDIR}/phast -name ".svn" | xargs rm -rf
	rm -r ${TMPDIR}/phast/doc ${TMPDIR}/phast/src/lib/rphast ${TMPDIR}/phast/test ${TMPDIR}/phast/binary_install.sh
	VERSION=`cat ${TMPDIR}/phast/version | sed 's/\./_/g'` ;\
	cd ${TMPDIR} ; tar cfz ${CWD}/phast.$$VERSION.tgz phast
	rm -rf ${TMPDIR}

doc::
	(cat Doxyfile; echo "PROJECT_NUMBER=`cat version`") | doxygen -

install:
	(cd src; make install DESTDIR=${DESTDIR} )
