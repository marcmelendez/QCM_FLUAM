/* Generate configurations of a polymer bead strands attached to a wall */

# include <stdio.h>
# include <math.h>
# include <stdint.h>
# include <stdlib.h>

/* Uncomment to ignore excluded volume */
# define EXCLUDED_VOLUME

/*** Auxiliary functions ***/
typedef enum {false, true} bool; /* Boolean type */

/* Pseudorandom number generation */
uint64_t s[] = {12679825035178159220ull, 15438657923749336752ull}; /* PRNG seed */

# define RANDOM_MAX 18446744073709551615ull
/* 64-bit (pseudo)random integer */
uint64_t xorshift128plus(void)
{
  uint64_t x = s[0];
  uint64_t const y = s[1];
  s[0] = y;
  x ^= x << 23; // a
  x ^= x >> 17; // b
  x ^= y ^ (y >> 26); // c
  s[1] = x;
  return x + y;
}

/* Random number from a uniform distribution */
float uniform(float min, float max)
{
  return min + (xorshift128plus()/((float) RANDOM_MAX))*(max - min);
}

/* Random number from a Gaussian distribution */
float Gaussian(float mean, float stddev)
{
  double u1, u2, s0 = 2;
  static bool saved_value = false; /* Flag indicating stored number */
  static double Gnum; /* Stored Gaussian number */

  if(saved_value) {
    saved_value = false;
    return Gnum;
  }

  while(s0 >= 1)
  {
    u1 = uniform(-1, 1);
    u2 = uniform(-1, 1);
    s0 = u1*u1 + u2*u2;
  }

  Gnum = mean + stddev*u2*sqrt(-2*logf(s0)/s0);
  saved_value = true;

  return mean + stddev*u1*sqrt(-2*logf(s0)/s0);
}

/* Marsaglia's algorithm for random point on a unit sphere */
void random_axis(double * n)
{
  float x, y, s = 2;

  while(s > 1) {
    x = uniform(-1, 1);
    y = uniform(-1, 1);
    s = x*x + y*y;
  }

  n[0] = 2*x*sqrt(1 - s);
  n[1] = 2*y*sqrt(1 - s);
  n[2] = 1 - 2*s;

  return;
}

/* Random initial position of the first monomer */
void first_monomer(float linker_energy, double * q)
{
  if(linker_energy == 0.0f) {
    random_axis(q);
    return;
  }
  float stddev = 1.0/sqrt(linker_energy);
  float theta = Gaussian(0.0f, stddev);
  float phi = uniform(0.0f, 2*M_PI);
  q[0] = sin(theta)*cos(phi);
  q[1] = sin(theta)*sin(phi);
  q[2] = sqrt(1 - q[0]*q[0] - q[1]*q[1]);

  return;
}

/* Define tabula rasa rule */
# define TABULA_RASA \
    first_monomer(linker_energy, q); \
    for(k = 0; k < 3; k++) { \
      p[k] = q[k]; \
      conf[k] = q[k]; \
    } \
    i = 1;


/* Integration step function declaration */
void integration_step(double * q, double * p, double gamma, double dt);

int main(int argc, char * argv[])
{
  /* Usage message */
  if(argc < 4) {
    printf("Usage: %s <number of configurations> <number of monomers> <persistence length>\n\n"
           "Distances are all measured in monomer diameters.\n", argv[0]);
    return 0;
  }

  /* Variables */
  int n; /* Configuration number */
  int i, j; /* Monomer indices */
  int k; /* Coordinate index */
  double q[3], p[3]; /* Tracking particle position and momentum */
  int nconfigurations = atoi(argv[1]); /* Number of configurations to generate */
  int nmonomers = atoi(argv[2]); /* Number of monomers in the DNA strand */
  double conf[3*nmonomers]; /* Configuration of strand + liposome */
  double lp = atof(argv[3]); /* Persistence length */
  double gamma = 2.0/lp; /* Friction coefficient */
  /* Assume kT = M = 1 */
  float linker_energy = 0.0; /* Energy of angular spring linking DNA strand to wall */
  double dt = 0.01; /* Time step */
  double rij[3]; /* Displacement vector */
  double r2; /* Distance squared */

  /* Loop over configurations */
  for(n = 0; n < nconfigurations; n++) {
    /* Initial conditions */
    TABULA_RASA

    /* Generate the positions of the monomers */
    while(i < nmonomers) {
      /* Move tracking particle position */
      integration_step(q, p, gamma, dt);

      /* Check distance to previous monomer */
      for(k = 0; k < 3; k++) rij[k] = q[k] - conf[3*(i-1) + k];
      r2 = 0; for(k = 0; k < 3; k++) r2 += rij[k]*rij[k];
      if(r2 >= 1) { /* New monomer position */
        for(k = 0; k < 3; k++) conf[3*i + k] = q[k];

        /* Check for overlap with wall */
        if(conf[3*i + 2] < .5) {
          /* Start over */
          TABULA_RASA
          continue;
        }

        # ifdef EXCLUDED_VOLUME
        /* Check for overlaps with other monomers */
        for(j = 0; j < i - 1; j++) {
          for(k = 0; k < 3; k++) rij[k] = conf[3*j + k] - conf[3*i + k];
          r2 = 0; for(k = 0; k < 3; k++) r2 += rij[k]*rij[k];
          if(r2 < 1) {
            /* Start over */
            j = -1;
            TABULA_RASA
            break;
          }
        }
        # endif
        /* If we broke out of the previous loop, we need to return to the beginning of the while loop */
        if(j < 0) continue;

        i++; /* Move on to next monomer */
      }
    }

    /* Output result (only DNA strands) */
    for(i = 0; i < nmonomers; i++)
      printf("%f\t%f\t%f\n", conf[3*i], conf[3*i + 1], conf[3*i + 2]);
    printf("\n");
  }

  return 0;
}

/*** Integration scheme for the Langevin stochastic differential equation ***

  The function below implements the algorithm in [1] for the special case in
  which the force vanishes.

    q(t + dt) = q(t) + b v(t) dt + b F(q(t)) dt^2/(2m)
                  + b zeta(t + dt) dt/(2m),
    v(t + dt) = v(t) + (F(q(t)) + F(q(t + dt))) dt/(2m)
                  - gamma/m (q(t + dt) - q(t)) + zeta(t + dt)/m.

  The step dt represents the difference between two consecutive times, q and v
  stand for the coordinates and velocities, respectively, F(q(t)) is the force
  at time t, m the particle mass, gamma the friction coefficient and zeta[i]
  a random number that satisfies

    <zeta[i]> = 0,
    <zeta[i] zeta[j]> = 2 gamma kB T dt delta[i, j],

  where kB refers to Boltzmann's constant and delta[i, j] to the Kronecker delta
  function. Furthermore, the constant b equals

    b = 1/(1 + gamma dt/(2m)).

  WARNING: Please remember that in this code, F = 0.

  -----
  References:

  [1] N. Gronbech-Jensen, and O. Farago: "A simple and effective Verlet-type
  algorithm for simulating Langevin dynamics", Molecular Physics (2013).
  http://dx.doi.org/10.1080/00268976.2012.760055 */

void integration_step(double * q, double * p, double gamma, double dt)
{
  int i; /* Index */
  double stddev = sqrt(2.0*gamma*dt); /* Standard deviation of distribution */
  double zeta[3]; /* Gaussian random numbers */

  /* Generate Gaussian random numbers */
  for(i = 0; i < 3; i++) zeta[i] = Gaussian(0, stddev);

  for(i = 0; i < 3; i++) {
    p[i] += gamma*q[i] + zeta[i];
    q[i] += (p[i] - gamma*q[i])*dt/(1 + 0.5*gamma*dt);
    p[i] -= gamma*q[i];
  }

//  printf("%f %f %f\n", q[0], q[1], q[2]);

  return;
}


