# .c->preproccesor(handle #include,#define)-compiler(c to assembly)->.s(assembly)->translate assembly to object file(particular to an OS)->.o->linker brings together object file to produce executable

all: program1

program1: matrix.o test.o
	gcc -static-libasan -Wall -pedantic     -Werror -g matrix.o  test.o  -o program1

matrix.o: matrix.c
	gcc -static-libasan  -Wall -pedantic   -Werror -g -c matrix.c -o matrix.o

test.o: test.c
	gcc -static-libasan  -Wall -pedantic  -static-libasan  -Werror -g -c test.c -o test.o


clean:
	rm -f *.o program1
