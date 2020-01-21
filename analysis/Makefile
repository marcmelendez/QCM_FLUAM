all: strand_forces strand_structure strand_stress
strand_forces:
	gcc -Ofast -march=native -Wall strand_forces.c -o strand_forces -lm
strand_structure:
	gcc -Ofast -march=native -Wall strand_structure.c -o strand_structure -lm
strand_stress:
	gcc -Ofast -march=native -Wall strand_stress.c -o strand_stress -lm
