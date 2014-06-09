all:	analysis_lib.a analysis_lib_test

analysis_lib_test:	analysis_lib.o analysis_lib_test.c makefile
	gcc -O3 -lm -o analysis_lib_test analysis_lib_test.c analysis_lib.a
analysis_lib.a:	analysis_lib.h analysis_lib.c makefile
	gcc -O3 -c analysis_lib.c
	ar rcs analysis_lib.a analysis_lib.o
	rm -f analysis_lib.o
