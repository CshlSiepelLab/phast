SUB = lib dless exoniphy phastCons phastOdds phastMotif phyloFit phyloBoot phyloP prequel util 

CDIR = ${PWD}
INSTALLDIR = /cluster/bin/phast/${MACHTYPE}

all: 
	mkdir -p ../bin ../lib
	for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} ; done
# for Windows builds, uncomment the following line to add .exe to all binary files generated
#	for file in `ls ../bin` ; do mv ../bin/$$file ../bin/$$file.exe ; done

install:
	@if [ "${MACHTYPE}" == "" ] ; then echo "ERROR: MACHTYPE must be defined." ; exit 1 ; fi
	mkdir -p ${INSTALLDIR}
	cp ../bin/* ${INSTALLDIR}

rphast:
	cd lib/rphast/src; ln -s `find ../.. -name "*.c" | grep -v "rphast/"` ../../../../include/*.h .

doc:
	doxygen

clean:
	@for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} clean ; done
	rm -f ../bin/* ../lib/*
