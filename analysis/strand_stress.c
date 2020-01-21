/* strand_stress.c:
   Calculate the stress on the QCM wall at z = 0 due to an immersed semiflexible
   polymer bead strand for a given configuration. */
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

int main(int argc, char * argv[])
{
  /* Usage message */
  if(argc < 5) {
    printf("Usage: %s <number of configurations> <number of monomers> <shear viscosity> <penetration> [<configuration file>]\n", argv[0]);
    return 0;
  }
  /* Variable declarations */
  int Nc = atoi(argv[1]); /* Number of polymer strand configurations */
  int N = atoi(argv[2]); /* Number of particles in the polymer strand */
  int eta = atof(argv[3]); /* Shear viscosity */
  int delta = atof(argv[4]); /* Penetration of oscillating flow */
  FILE * data_file = stdin; /* Name of data file (defaults to stdin) */
  char text_buffer[250];
  int nc; /* Configuration number */
  int i; /* Particle index */
  int m; /* Character index */
  int nread; /* Number of items read from data */

  double positions[3*N]; /* Monomer positions array */
  double forces[3*N]; /* Array to store the forces on each monomer */
  double stress; /* Shear stress */

  /* Open file if file name provided */
  if(argc > 5) {
    data_file = fopen(argv[5], "r");
    if(data_file == NULL) {
      fprintf(stderr, "Error: file not found: %s.\n", argv[5]);
      return -1;
    }
  }

  /* Loop over configurations */
  for(nc = 0; nc < Nc; nc++) {
    /* Read in positions and forces */
    for(i = 0; i < N; i++) {
      if(!fgets(text_buffer, 250, data_file)) {
        fprintf(stderr, "Error: end of initial conditions file reached prematurely.\n");
        exit(-1);
      }

      /* Ignore whitespace and comments */
      m = 0; while((text_buffer[m] == ' ') || (text_buffer[m] == '\t')) m++;
      if(text_buffer[m] != '#' && text_buffer[m] != '\n') {
        nread = sscanf(text_buffer, "%lf %lf %lf %lf %lf %lf",
                       &(positions[3*i]), &(positions[3*i + 1]), &(positions[3*i + 2]),
                       &(forces[3*i]), &(forces[3*i + 1]), &(forces[3*i + 2]));
      }
      else i--;
      if(nread < 6) {
        fprintf(stderr, "Error: insufficient rows for particle %d.\n", i);
        exit(-1);
      }
    }

    /* Clear stress */
    stress = 0;

    /* Calculate stress */
    for(i = 0; i < N; i++)
      stress += forces[3*i]*eta*exp(-positions[3*i + 2]/delta)*cos(positions[3*i + 2]/delta)/(4*M_PI);

    /* Output stress */
    printf("%f\n", stress);
  }

  /* Clean up */
  fclose(data_file);

  return 0;
}
