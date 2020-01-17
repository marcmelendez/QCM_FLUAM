/* strand_forces.c:
   Calculate the angular forces on every bead of a semiflexible polymer bead
   strand for a given configuration. */
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

int main(int argc, char * argv[])
{
  /* Usage message */
  if(argc < 4) {
    printf("Usage: %s <number of configurations> <number of monomers> <angular spring constant> [<configuration file>]\n", argv[0]);
    return 0;
  }
  /* Variable declarations */
  int Nc = atoi(argv[1]); /* Number of polymer strand configurations */
  int N = atoi(argv[2]); /* Number of particles in the polymer strand */
  int kang = atof(argv[3]); /* Angular bond spring constant */
  FILE * positions_file = stdin; /* Name of configuration file (defaults to stdin) */
  char text_buffer[250];
  int nc; /* Configuration number */
  int i; /* Particle index */
  int l; /* Coordinate index */
  int n; /* Bond index */
  int m; /* Character index */
  int nread; /* Number of items read from data */

  double positions[3*N]; /* Monomer positions array */
  double forces[3*N]; /* Array to store the forces on each monomer */
  double rij[3], rjk[3]; /* Separation vectors */
  double r2_1, r2_2; /* Square distances */
  double invr_1, invr_2; /* Inverse distances */
  double costheta; /* Cosine of bending angle */
  double Fmod; /* Magnitude of force vector */
  double Fi, Fk; /* Force terms */

  /* Open file if file name provided */
  if(argc > 4) {
    positions_file = fopen(argv[4], "r");
    if(positions_file == NULL) {
      fprintf(stderr, "Error: file not found: %s.\n", argv[3]);
      return -1;
    }
  }

  /* Loop over configurations */
  for(nc = 0; nc < Nc; nc++) {
    /* Read in configuration */
    for(i = 0; i < N; i++) {
      if(!fgets(text_buffer, 250, positions_file)) {
        fprintf(stderr, "Error: end of initial conditions file reached prematurely.\n");
        exit(-1);
      }

      /* Ignore whitespace and comments */
      m = 0; while((text_buffer[m] == ' ') || (text_buffer[m] == '\t')) m++;
      if(text_buffer[m] != '#' && text_buffer[m] != '\n') {
        nread = sscanf(text_buffer, "%lf %lf %lf",
                       &(positions[3*i]), &(positions[3*i + 1]), &(positions[3*i + 2]));
      }
      else i--;
      if(nread < 3) {
        fprintf(stderr, "Error: insufficient rows for particle %d.\n", i);
        exit(-1);
      }
    }

    /* Clear force vector */
    for(i = 0; i < 3*N; i++) forces[i] = 0;

    /* Calculate forces */
    for(n = 0; n < N - 2; n++) {
      /* Vectors rij and rjk */
      for(l = 0; l < 3; l++) rij[l] = positions[3*(n + 1) + l] - positions[3*n + l];
      for(l = 0; l < 3; l++) rjk[l] = positions[3*(n + 2) + l] - positions[3*(n + 1) + l];

      /* Distances */
      r2_1 = r2_2 = 0;
      for(l = 0; l < 3; l++) {
        r2_1 += rij[l]*rij[l];
        r2_2 += rjk[l]*rjk[l];
      }

      /* Inverse distances */
      invr_1 = 1/sqrt(r2_1);
      invr_2 = 1/sqrt(r2_2);

      /* Cosine of angle */
      costheta = 0; for(l = 0; l < 3; l++) costheta += rij[l]*rjk[l]; costheta *= invr_1*invr_2;

      /* Force (magnitude) */
      Fmod = 0;
      # ifdef SMALL_ANGLE_BENDING
      double theta = acos(costheta);
      if(costheta*costheta <= 1 && theta > 0)
        Fmod = kang*theta/sin(theta);
      # else
        Fmod = kang;
      # endif

      /* Normalise separation vectors */
      for(l = 0; l < 3; l++) {
        rij[l] *= invr_1;
        rjk[l] *= invr_2;
      }

      /* Force (components) */
      for(l = 0; l < 3; l++) {
        Fi = -Fmod*(rjk[l] - costheta*rij[l])*invr_1;
        Fk = -Fmod*(rij[l] - costheta*rjk[l])*invr_2;
        forces[3*n + l] += Fi;
        forces[3*(n + 1) + l] -= (Fi - Fk);
        forces[3*(n + 2) + l] -= Fk;
      }
    }

    /* Output positions and forces */
    for(i = 0; i < N; i++)
      printf("%f %f %f %f %f %f\n", positions[3*i], positions[3*i + 1], positions[3*i + 2],
                                    forces[3*i], forces[3*i + 1], forces[3*i + 2]);
    printf("\n");
  }

  return 0;
}
