include make-include.mk

SUB = lib dless exoniphy phastCons phastOdds phastMotif phyloFit phyloBoot phyloP prequel util 

CDIR = ${PWD}

all: 
	mkdir -p ../bin ../lib
	for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} ; done	
ifeq ($(TARGETOS), Windows)
	@for file in `ls ../bin` ; do mv ../bin/$$file ../bin/$$file.exe ; done
endif

install:
ifndef DESTDIR
	DESTDIR=/
endif
	mkdir -p ${DESTDIR}/usr/bin/
	mkdir -p ${DESTDIR}/opt/phast/data/
	cp ../bin/* ${DESTDIR}/usr/bin/
	cp -R ../data/* ${DESTDIR}/opt/phast/data/

doc:
	doxygen

clean:
	@for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} clean ; done
	rm -f ../bin/* ../lib/*
