/* elastic_network.c: Output a list of elastic newtork model bonds */

/* Standard library */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

int main(int argc, char * argv[])
{
  if(argc < 8) { /* Use message */
    printf("Usage: %s <N> <Rc> <K> <Lx> <Ly> <Lz> <file> [DIM] \nParameters:\n", argv[0]);
    printf("\tN:\tNumber of particles to read from file (set to -1 to read from the top of the input file).\n"
           "\tRc:\tCut-off radius for bonds.\n"
           "\tK:\tBond strength parameter.\n"
           "\tLx, Ly, Lz:\tBox dimensions "
           "(enter -1 for no periodic boundary conditions).\n"
           "\tfile:\tFilename of particle positions "
           "(Format: x y z ... by rows).\n"
           "\tDIM:\tDimensionality of space (1, 2 or 3).\n");
    return 0;
  }

  int DIM; /* Dimensionality */
  if(argc == 9) DIM = atoi(argv[8]); else DIM = 3;
  int N = atoi(argv[1]); /* Number of particles */
  FILE * data = fopen(argv[7], "r"); /* Pointer to data file */
  char line[250]; /* Line buffer */
  int nread; /* Number of fields read */
  if(N < 0) {
    fgets(line, 250, data);
    nread = sscanf(line, "%d", &N);
    fprintf(stderr, "Number of particles: %d\n", N);
  }
  float pos[N][DIM]; /* Array of particle positions */
  float L[3]; /* Box dimensions */
  int periodicL[3] = {0, 0, 0}; /* Flag indicating the periodic directions */
  float Rc = atof(argv[2]); /* Cut-off radius */
  float K = atof(argv[3]); /* Bond strength parameter */
  int i, j; /* Particle indices */
  int k, l; /* Coordinate indices */
  int m, nm[3]; /* Cell indices */
  double r; /* Distance between particles */
  double r2; /* Distance squared */
  double rij[DIM]; /* Displacement */
  int ncells = 1; /* Number of neighbour cells */
  int llist[N]; /* Linked list for particles */
  int n[DIM]; /* Number of cells in each direction */
  float cellsize[DIM]; /* Size of the cell in each direction */
  int cellcoord[DIM], ncellcoord[DIM]; /* Cell (integer) coordinates */
  int cellidx, tmpidx; /* Cell indices */

  /* Set dimensions of box */
  for(k = 0; k < DIM; k++) {
    if(atof(argv[4 + k]) > 0) {
      L[k] = atof(argv[4 + k]);
      periodicL[k] = 1;
    }
  }

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
    switch(DIM) {
      case 1:
        nread = sscanf(line, "%f", &pos[i][0]);
        break;
      case 2:
        nread = sscanf(line, "%f %f", &pos[i][0], &pos[i][1]);
        break;
      case 3:
        nread = sscanf(line, "%f %f %f", &pos[i][0], &pos[i][1], &pos[i][2]);
        break;
      default:
        printf("Error: invalid dimensionality.\n");
        return -1;
    }
    if(nread == DIM) {
      /* Set box size if not specified */
      for(k = 0; k < DIM; k++) {
        if(periodicL[k] == 0) {
          if(fabs(pos[i][k]) > L[k]) L[k] = 2*fabs(pos[i][k]);
        }
      }
      i++; /* Move on to next particle */
      if(i >= N) break;
    }
  }

  for(k = 0; k < DIM; k++) {
    if(L[k] < 3*Rc) L[k] = 3*Rc;
  }

  /* Have we read the positions of all the particles? */
  if(i < N) {
    printf("Error: end of file reached prematurely.\n");
    return -1;
  }

  /* Calculate the number of cells and cell size */
  for(k = 0; k < DIM; k++) {
    n[k] = (int) floorf(L[k]/Rc);
    cellsize[k] = L[k]/((float) n[k]);
    ncells *= n[k];
  }
  int hlist[ncells]; /* Head list for cells */

  /* Initialise head list */
  for(m = 0; m < ncells; m++) hlist[m] = -1;

  /* Build the linked list */
  for(i = 0; i < N; i++) {
    /* Cell coordinates of position */
    for(k = 0; k < DIM; k++) {
      cellcoord[k] = (int) floorf((pos[i][k] + .5*L[k])/cellsize[k]);
      if(cellcoord[k] == n[k]) cellcoord[k] = 0;
    }

    /* Get cell index */
    cellidx = 0; /* Clear cell index */
    for(k = DIM; k > 0; k--) {
      tmpidx = cellcoord[k - 1];
      for(l = k - 1; l > 0; l--)
        tmpidx *= n[l - 1];
      cellidx += tmpidx;
    }

    /* Add particle to linked list */
    llist[i] = hlist[cellidx];
    hlist[cellidx] = i;
  }

  /* Loop over all particle indices working out which ones are in bond range Rc */
  for(i = 0; i < N; i++) {
    /* Get position cell number */
    for(k = 0; k < DIM; k++) {
      cellcoord[k] = (int) floorf((pos[i][k] + .5*L[k])/cellsize[k]);
    }
    /* Loop over neighbour cells */
    nm[0] = nm[1] = nm[2] = -1;
    for(m = 0; m < pow(3,DIM); m++) {
      /* Neighbour cell coordinates */
      for(k = 0; k < DIM; k++) ncellcoord[k] = cellcoord[k] + nm[k];

      /* Periodic boundary conditions */
      for(k = 0; k < DIM; k++) {
        if(ncellcoord[k] < 0) ncellcoord[k] += n[k];
        if(ncellcoord[k] >= n[k]) ncellcoord[k] -= n[k];
      }

      /* Determine neighbour cell index */
      cellidx = 0; /* Clear cell index */
      for(k = DIM; k > 0; k--) {
        tmpidx = ncellcoord[k - 1];
        for(l = k - 1; l > 0; l--)
          tmpidx *= n[l - 1];
        cellidx += tmpidx;
      }

      /* Go through linked lists checking distances */
      j = hlist[cellidx];
      while(j >= 0) {
        if(j >= i) {j = llist[j]; continue;} /* No bonds from a particle to itself */

        /* Displacement */
        for(k = 0; k < 3; k++) rij[k] = pos[j][k] - pos[i][k];

        /* Periodic boundary conditions */
        for(k = 0; k < DIM; k++) {
          if(periodicL[k] == 1)
            rij[k] -= floorf(rij[k]/L[k] + .5)*L[k];
        }

        /* Distance */
        r2 = 0; for(k = 0; k < DIM; k++) r2 += rij[k]*rij[k];
        r = sqrt(r2);

        /* Output bond */
        if(r <= Rc) {
          printf("%d\t%d\t%f\t%f\n", i, j, K, r);
        }

        j = llist[j];
      }

      /* Next neighbour cell */
      nm[0]++;
      if(nm[0] > 1) {nm[0] = -1; nm[1]++;}
      if(nm[1] > 1) {nm[1] = -1; nm[2]++;}
    }
  }

  return 0;
}
