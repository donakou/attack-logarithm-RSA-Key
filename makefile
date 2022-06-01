CSLD: main.o log.o rho.o
	gcc -g main.o log.o rho.o -lgmp -o CSLD
main.o: main.c log.c rho.c
	gcc -c -Wall main.c -lgmp -o main.o
log.o: log.c
	gcc -c -Wall log.c -lgmp -o log.o
rho.o: rho.c
	gcc -c -Wall rho.c -lgmp -o rho.o
