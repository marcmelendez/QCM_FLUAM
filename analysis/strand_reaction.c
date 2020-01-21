/* strand_reaction.c:
   Calculate the forces on the fluid due to a semiflexible polymer bead strand
   in a Stokes boundary layer flow. */
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <complex.h>

/* Calculate the forces on each bead of the strand */
void strand_forces(int N, double kbond, double kang, double * positions, double * forces);

int main(int argc, char * argv[])
{
  /* Usage message */
  if(argc < 8) {
    printf("Usage: %s <number of configurations> <number of monomers> <harmonic bond constant> <angular spring constant> <fluid density> <fluid viscosity> <oscillation frequency> [<configuration file>]\n", argv[0]);
    return 0;
  }
  /* Variable declarations */
  int Nc = atoi(argv[1]); /* Number of polymer strand configurations */
  int N = atoi(argv[2]); /* Number of particles in the polymer strand */
  double kbond = atof(argv[3]); /* Harmonic bond spring constant */
  double kang = atof(argv[4]); /* Angular bond spring constant */
  double
  FILE * positions_file = stdin; /* Name of configuration file (defaults to stdin) */
  char text_buffer[250];
  int nc; /* Configuration number */
  int i, j; /* Particle indices */
  int m; /* Character index */
  int nread; /* Number of items read from data */

  double positions[3*N]; /* Monomer positions array */
  double forces[3*N]; /* Array to store the forces on each monomer */
  double forces_new[3*N]; /* Changed forces due to a small displacement */
  double K[N*N]; /* Connectivity matrix */
  complex x[N]; /* Horizontal amplitude phasor */
  complex F[N]; /* Force phasor */

  /* Open file if file name provided */
  if(argc > 5) {
    positions_file = fopen(argv[5], "r");
    if(positions_file == NULL) {
      fprintf(stderr, "Error: file not found: %s.\n", argv[5]);
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

    /* Calculate elastic network connectivity matrix for this configuration */
    /* Clear the connectivity matrix */
    for(i = 0; i < N*N; i++) K[i] = 0;

    /* Base strand configuration forces */
    strand_forces(N, kbond, kang, positions, forces);

    /* Loop over particles */
    for(i = 0; i < N ; i++) {
      # define DELTA_X 0.1d /* Magnitude of the displacements */

      /* Small displacement of particle i in the x direction */
      positions[3*i] += DELTA_X;

      /* Calculate the forces in the new configuration */
      strand_forces(N, kbond, kang, positions, forces_new);

      /* Calculate the elements of the connectivity matrix */
      for(j = i; j < N; j++) K[N*i + j] = (forces_new[3*j] - forces[3*j])/DELTA_X;

      /* The connectivity matrix is symmetric */
      for(j = i + 1; j < N; j++) K[N*j + i] = K[N* i + j];

      /* Undo the displacement of particle i in the x direction */
      positions[3*i] += DELTA_X;
    }
  }

  /* Calculate position phasors for all the points on the strand  */
  


  /* Clean up */
  fclose(positions_file);

  return 0;
}

/* Calculate the forces on each bead of the strand */
void strand_forces(int N, double kbond, double kang, double * positions, double * forces) {
  double rij[3], rjk[3]; /* Separation vectors */
  double r2_1, r2_2; /* Square distances */
  double invr_1, invr_2; /* Inverse distances */
  double costheta; /* Cosine of bending angle */
  double Fmod; /* Magnitude of force vector */
  double Fi, Fk; /* Force terms */

  int l; /* Coordinate index */
  int n; /* Bond index */

    /* Clear force vector */
    for(l = 0; l < 3*N; l++) forces[l] = 0;

    /* Calculate forces */

    /* Harmonic bonds */
    for(n = 0; n < N - 1; n++) {
      /* Separation vector rij */
      for(l = 0; l < 3; l++) rij[l] = positions[3*(n + 1) + l] - positions[3*n + l];

      /* Distance */
      r2_1 = 0; for(l = 0; l < 3; l++) r2_1 += rij[l]*rij[l];

      /* Inverse distance */
      invr_1 = 1/sqrt(r2_1);

      /* Forces */
      for(l = 0; l < 3; l++) {
        forces[3*n + l] += kbond*(r2_1 - 1)*rij[l]*invr_1;
        forces[3*(n + 1) + l] -= kbond*(r2_1 - 1)*rij[l]*invr_1;
      }
    }

    /* Angular springs */
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

  return;
}
