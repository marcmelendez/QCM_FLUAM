/* elastic_network.c: Output a list of elastic newtork model bonds
   This code uses a naïve N-body check. Use elastic_network.c for
   a faster method based on neighbour lists.
*/

/* Standard library */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define DIM 3 /* Dimensionality */
/* Read command (which depends on the dimensionality) */
# if DIM == 1
# define READ_DATA "%f", &pos[i][0]
# elif DIM == 2
# define READ_DATA "%f %f", &pos[i][0], &pos[i][1]
# else
# define READ_DATA "%f %f %f", &pos[i][0], &pos[i][1], &pos[i][2]
# endif

int main(int argc, char * argv[])
{
  if(argc < 8) { /* Use message */
    printf("Use: %s <N> <Rc> <K> <Lx> <Ly> <Lz> <file>\nParamters:\n", argv[0]);
    printf("\tN:\tNumber of particles to read from file,\n");
    printf("\tRc:\tCut-off radius for bonds.\n");
    printf("\tK:\tBond strength parameter.\n");
    printf("\tLx, Ly, Lz:\tBox dimensions "
           "(enter -1 for no periodic boundary conditions).\n");
    printf("\tfile:\tFilename of particle positions "
           "(Format: x y z ... by rows).\n");
    return 0;
  }

  int N = atoi(argv[1]); /* Number of particles */
  float pos[N][DIM]; /* Array of particle positions */
  float L[3]; /* Box dimensions */
  float Rc = atof(argv[2]); /* Cut-off radius */
  float K = atof(argv[3]); /* Bond strength parameter */
  FILE * data = fopen(argv[7], "r"); /* Pointer to data file */
  char line[250];
  int i, j; /* Particle indices */
  int k; /* Coordinate index */
  int nread; /* Number of fields read */
  double r; /* Distance between particles */
  double xij; /* Coordinate */

  /* Check file pointer */
  if(data == NULL) {
    printf("Error: unable to open file %s.\n", argv[7]);
    return -1;
  }

  /* Get box dimensions */
  for(k = 0; k < DIM; k++) L[k] = atof(argv[4 + k]);

  i = 0;
  while(fgets(line, 250, data)) {
    /* Read in the position of particle i */
    nread = sscanf(line, READ_DATA);
    if(nread == DIM) {
      /* Bond new particle to previous particles */
      for(j = 0; j < i; j++) {
        /* Calculate distance between particles */
        r = 0;
        for(k = 0; k < DIM; k++) {
          xij = (pos[j][k] - pos[i][k]);
          if(L[k] > 0) xij -= floorf(xij/L[k] + 0.5f)*L[k];
          r += xij*xij;
        }
        r = sqrt(r);

        /* Link particles i and j if they are close enough */
        if(r < Rc)
          printf("%d %d %f %f\n", i, j, K, r);
      }

      i++; /* Move on to next particle */
      if(i >= N) break;
    }
  }

  return 0;
}
