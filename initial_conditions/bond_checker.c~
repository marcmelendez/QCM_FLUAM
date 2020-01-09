/* bond_checker.c: Check a list of bonds */

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
  if(argc < 7) { /* Use message */
    printf("Use: %s <N> <Lx> <Ly> <Lz> <particle file> <bond file> [tolerance] [nbonds]\nParamters:\n", argv[0]);
    printf("  N:\t\tNumber of particles to read from file (set to -1 to read from the first line in the particle file),\n");
    printf("  Lx, Ly, Lz:\tBox dimensions "
           "(enter -1 for no periodic boundary conditions).\n");
    printf("  particle file:\tFilename of particle positions "
           "(Format: x y z ... by rows).\n"
           "  bond file:\tFilename containing bonds "
           "(Format: i j k r0, by rows --\n\t\t\ti, j: particle index, k: spring constant, r0: equilibrium length)\n"
           "  tolerance:\tAdmissible separation divided by r0 (default value: 2)\n"
           "  nbonds:\t\tRead in the number of bonds on the first line of the bonds file\n"
           "\t\t\t(use any value of this argument to activate)\n");
    return 0;
  }

  float L[3]; /* Box dimensions */
  FILE * positions = fopen(argv[5], "r"); /* Pointer to data file */
  FILE * bonds = fopen(argv[6], "r"); /* Pointer to data file */
  char line[250];
  int i, j; /* Particle indices */
  int l; /* Line number */
  int k; /* Coordinate index */
  float kspring, r0; /* Bond parameters (spring strength and equilibrium distance) */
  float tol; /* Tolerance (admissible separation divided by r0) */
  int nread; /* Number of fields read */
  double r, r2; /* Distance between particles and distance squared */
  double xij; /* Coordinate */

  /* Set tolerance */
  if(argc > 7) tol = atof(argv[7]);
  else tol = 2;

  /* Check file pointers */
  if(positions == NULL) {
    printf("Error: unable to open file %s.\n", argv[5]);
    return -1;
  } else if (bonds == NULL) {
    printf("Error: unable to open file %s.\n", argv[6]);
    return -1;
  }

  /* Get box dimensions */
  for(k = 0; k < DIM; k++) L[k] = atof(argv[2 + k]);

  int N = atoi(argv[1]); /* Number of particles */
  int nbonds = -1; /* Number of bonds to read (if nbonds > 0) */
  /* Read number of particles */
  if(N < 0) {
    fgets(line, 250, positions);
    nread = sscanf(line, "%d", &N);
    fprintf(stderr, "Number of particles: %d\n", N);
  }
  float pos[N][DIM]; /* Array of particle positions */

  /* Get positions data */
  i = 0;
  while(fgets(line, 250, positions)) {
    /* Read in the position of particle i */
    nread = sscanf(line, READ_DATA);
    i++; /* Move on to next particle */
    if(i >= N) break;
  }

  /* Get number of bonds from first line */
  if(argc > 8) {
    fgets(line, 250, bonds);
    nread = sscanf(line, "%d", &nbonds);
    fprintf(stderr, "Number of bonds: %d\n", nbonds);
  }

  /* Check bonds */
  l = 1;
  while(fgets(line, 250, bonds)) {
    /* Read in bond information */
    nread = sscanf(line, "%d %d %f %f", &i, &j, &kspring, &r0);

    /* Check number of elements read */
    if(nread < 4) printf("Error reading line %d (%d elements read).\n", l, nread);

    /* Calculate the distance between particles i and j */
    r2 = 0;
    for(k = 0; k < DIM; k++) {
      xij = (pos[j][k] - pos[i][k]);
      if(L[k] > 0) xij -= floorf(xij/L[k] + 0.5f)*L[k];
        r2 += xij*xij;
      }
    r = sqrt(r2);

    /* Check separation */
    if(r > tol*r0) printf("Large strain between particles %d and %d (r = %f, r0 = %f).\n", i, j, r, r0);
    if(r < r0/tol) printf("Strong compression between particles %d and %d (r = %f, r0 = %f).\n", i, j, r, r0);

# ifdef DEBUG
    /* Info read */
    if(r > tol*r0 || r < r0/tol) {
      printf("Bonded particles:\n  r_%d = (%f, %f, %f), r_%d = (%f, %f, %f)\n", i, pos[i][0], pos[i][1], pos[i][2], j, pos[j][0], pos[j][1], pos[j][2]);
      printf("Bond parameters: k = %f, r0 = %f\n", kspring, r0);
    }
# endif


    l++; /* Move on to next line */
    if(nbonds >= 0 && l > nbonds) break;
  }

  printf("Read %d lines.\n", l - 1);

  return 0;
}
