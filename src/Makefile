include make-include.mk

SUB = lib dless exoniphy phastCons phastOdds phastMotif phyloFit phyloBoot phyloP prequel util 

CDIR = ${PWD}

all: 
	mkdir -p ../bin ../lib
	for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} ; done	
ifeq ($(TARGETOS), Windows)
	@for file in `ls ../bin` ; do mv ../bin/$$file ../bin/$$file.exe ; done
else
	mkdir -p ../doc/man/
	@echo "Generating man pages..."
	@for file in `ls ../bin` ; do perl help2man.pl ../bin/$$file ../doc/man/ ; done
	@echo "Done."
endif

sharedlib:
	cd ${CDIR}/lib && ${MAKE} sharedlib 
	gcc ${CFLAGS} -shared ${LFLAGS} ${LIBPATH} -o libphast.so `find lib/ -name "*.o"` ${LIBS}
	mkdir -p ../lib/sharedlib
	mv libphast.so ../lib/sharedlib/

install:
ifndef DESTDIR
	DESTDIR=/
endif
	mkdir -p ${DESTDIR}/usr/bin/
	mkdir -p ${DESTDIR}/opt/phast/data/
	mkdir -p ${DESTDIR}/usr/share/man/man1/
	cp ../bin/* ${DESTDIR}/usr/bin/
	cp -R ../data/* ${DESTDIR}/opt/phast/data/
	cp -R ../doc/man/* ${DESTDIR}/usr/share/man/man1/

doc:
	cd ../; make doc 

clean:
	@for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} clean ; done
	rm -rf ../bin ../lib ../doc

manpages:
	@for file in `ls ../bin` ; do perl help2man.pl ../bin/$$file ; done
