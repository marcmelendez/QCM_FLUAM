all: configurations lattice elastic_network
configurations:
	gcc -Wall strand_configuration.c -o strand_configuration -lm -Ofast
	gcc -Wall liposome_configuration.c -o liposome_configuration -lm -Ofast
lattice:
	gcc -Wall lattice.c -o lattice -lm -Ofast
elastic_network:
	gcc -Wall elastic_network.c -o elastic_network -lm -Ofast
	gcc -Wall elastic_network.nbody.c -o elastic_network.nbody -lm -Ofast
	gcc -Wall bond_checker.c -o bond_checker -lm -Ofast
