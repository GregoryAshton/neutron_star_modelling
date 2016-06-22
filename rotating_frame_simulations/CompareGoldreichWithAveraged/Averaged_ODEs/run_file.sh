

# Compile the script and run a.out
gcc -Wall -I/usr/local/include -c Goldreich_averaged_ODE_solver.c
gcc -static Goldreich_averaged_ODE_solver.o -lgsl -lgslcblas -lm

./a.out > generic.txt 

#rm a.out *.o 

