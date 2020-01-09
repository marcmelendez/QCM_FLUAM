/***************************  lattice.c  ****************************

  Generation of Bravais lattices. The program calculates positions

         R = i a1 + j a2 + k a3

  within a simulation box of side length L and inserts the basis at
  each point of the lattice. */

/***** Standard libraries *****/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***** Help message *****/
void help_message(char name[])
{
  printf("-- Bravais lattices --\n\nUsage: %s [options]\n\n", name);
  printf("Options:\n"
         "--type -t type\t\t\tLattice type (sc, bcc, fcc, dia, hcp, sq or tri).\n"
         "--number -N integer\t\tNumber of lattice nodes to generate.\n"
         "--length -L number\t\tBox side length.\n"
         "-Lx -Ly -Lz number\t\tBox length, width and height.\n"
         "--radius -r number\t\tParticle radius.\n"
         "--colour --color -c integer\tParticle colour (rgb integer).\n"
         "--keep-aspect -k\t\tKeep aspect ratio.\n"
         "--fill -f number\t\tFill the box with particles separated by a distance = <number>.\n"
         "-b file\t\t\t\tUser defined basis file.\n"
         "-v file\t\t\t\tUser defined lattice vectors file.\n"
  );
  return;
}


/***** Bravais lattice types *****/
typedef enum {sc,bcc,fcc,dia,hcp,sq,tri} lattice;
typedef enum {false=0, true} bool;
#define MAXBASIS 1000 /* Maximum number of user basis elements */

/***** Auxiliary functions *****/

/* Hash for the last four letters in a word */
unsigned int hash(char word[])
{
  int i = 0; /* Index */
  unsigned long hash = 0; /* Hash value */

  /* Calculate hash */
  while(word[i] != '\0') {
    hash = (hash << 8) + word[i];
    i++;
  }

  return hash;
}


/***** Main function *****/

int main(int argc, char * argv[])
{
  int N = 1000, nx, ny, nz; /* Number of lattice nodes */
  int node = 0; /* Lattice node number */
  float V; /* Box volume */
  int ncells; /* Number of unit cells */
  int i, j, k, l, d, n; /* Indices */
  lattice type = sc; /* Lattice type */
  bool keepaspect = false; /* Do not keep aspect ratio (if unnecessary) */
  bool fillbox = false; /* No need to fill the box */
  float separation = -1; /* Minimum separation between two atoms */
  float cellsize = -1; /* Side length of unit cell */
  float stretchfactor[3], minstretch; /* Stretch factors */
  float L[3] = {1, 1, 1}; /* Box dimensions */
  float radius = -1; /* Particle radius */
  int colour = -1; /* Particle colour */
  float e[3][3]; /* Lattice vectors */
  float basis[8][3]; /* Basis positions */
  int nbasis = 1; /* Number of elements in the basis */
  float pbasis[MAXBASIS][5]; /* User basis elements */
  int npbasis = -1; /* Number of user basis elements */
  int nvectors = -1; /* Number of user lattice vectors */
  float r[3]; /* Position */
  FILE * basisfile, * vectorsfile; /* User basis and vectors files */
  char buffer[250]; /* Buffer to store file data */
  int count; /* Number of elements read */

  /*** Use message ***/
  if(argc < 2) {
    help_message(argv[0]);
    return 0;
  }

  /* Clear all the lattice vector components */
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      e[i][j] = 0;

  /* Clear personalised basis elements */
  for(j = 0; j < MAXBASIS; j++)
    for(k = 0; k < 5; k++)
      pbasis[j][k] = 0;

  /*** Read command-line options ***/
  for(i = 1; i < argc; i++) {
    switch(hash(argv[i])) {
      case 0x74797065: /* --type, -t Lattice type */
      case 0x2D74:
        i++;
        if(i < argc)
        switch(hash(argv[i])) {
          case 0x7363: /* Simple cubic (sc) */
            type = sc;
            break;
          case 0x626363: /* Body-centred cubic (bcc) */
            type = bcc;
            break;
          case 0x666363: /* Face-centred cubic (fcc) */
            type = fcc;
            break;
          case 0x646961: /* Diamond structure (dia)*/
            type = dia;
            break;
          case 0x686370: /* Hexagonal close-packed (hcp) */
            type = hcp;
            break;
          case 0x7371: /* Square (sq)*/
            type = sq;
            break;
          case 0x747269: /* Triangular (tri) */
            type = tri;
            break;
          default: /* Unrecognised lattice type */
          fprintf(stderr, "Error: Unrecognised lattice type: %s. "
                          "Using simple cubic.\n", argv[i]);
        }
        break;
      case 0x6D626572: /* --number, -N Number of lattice nodes */
      case 0x2D4E:
        i++;
        if(i < argc) N = atoi(argv[i]);
        break;
      case 0x6E677468: /* --length, -L Box side length */
      case 0x2D4C:
        i++;
        if(i < argc) L[0] = L[1] = L[2] = atof(argv[i]);
        break;
      case 0x2D4C78: /* -Lx Box length */
        i++;
        if(i < argc) L[0] = atof(argv[i]);
        break;
      case 0x2D4C79: /* -Ly Box width */
        i++;
        if(i < argc) L[1] = atof(argv[i]);
        break;
      case 0x2D4C7A: /* -Lz Box height */
        i++;
        if(i < argc) L[2] = atof(argv[i]);
        break;
      case 0x64697573: /* --radius, -r Particle radius */
      case 0x2D72:
        i++;
        if(i < argc) radius = atof(argv[i]);
        break;
      case 0x6C6F7572: /* --colour, --color, -c Particle colour */
      case 0x6F6C6F72:
      case 0x2D63:
        i++;
        if(i < argc) colour = atoi(argv[i]);
        break;
      case 0x70656374: /* --keep-aspect, -k Keep aspect ratio */
      case 0x2D6B:
        keepaspect = true;
        break;
      case 0x66696C6C: /* --fill, -f fill the box */
      case 0x2D66:
        i++;
        if(i < argc) separation = atof(argv[i]);
        fillbox = true;
        break;
      case 0x61736973: /* -b User basis. */
      case 0x2D62:
        i++;
        if(i < argc) {
          basisfile = fopen(argv[i],"r");
          if(basisfile == NULL) {
            printf("Error: file not found (%s).\n", argv[i]);
            return -1;
          }

          /* Read data from file */
          npbasis = 0;
          while(fgets(buffer, 250, basisfile) && npbasis < MAXBASIS) {
            count = sscanf(buffer, "%f %f %f %f %f", &pbasis[npbasis][0],
                                   &pbasis[npbasis][1], &pbasis[npbasis][2],
                                   &pbasis[npbasis][3], &pbasis[npbasis][4]);
            if(count > 0) npbasis++;
          }
          fclose(basisfile);
        }
        break;
      case 0x746F7273: /* --vectors -v User lattice vectors. */
      case 0x2D76:
        i++;
        if(i < argc) {
          vectorsfile = fopen(argv[i],"r");
          if(vectorsfile == NULL) {
            printf("Error: file not found (%s).\n", argv[i]);
            return -1;
          }

          /* Read data from file */
          nvectors = 0;
          while(fgets(buffer, 250, vectorsfile) && nvectors < 3) {
            count = sscanf(buffer, "%f %f %f", &e[nvectors][0], &e[nvectors][1],
                                               &e[nvectors][2]);
            if(count > 0) nvectors++;
          }
          fclose(vectorsfile);
        }
        break;
      default: /* Ignore unrecognised options */
        fprintf(stderr, "Unrecognised option: %s\n", argv[i]);
    }
  }

  /*** Definition of the lattice vectors (if undefined by user) ***/
  if(nvectors < 0) {
    /* Write non-zero components of lattice vectors */
    switch(type) {
      case sc:
      case bcc:
      case fcc:
      case dia:
        e[0][0] = e[1][1] = e[2][2] = 1;
        break;
      case hcp:
        e[0][0] = 1;
        e[1][0] = 0.5;  e[1][1] = sqrt(3)/2;
        e[2][2] = 2*sqrt(6)/3;
        break;
      case sq:
        e[0][0] = e[1][1] = 1;
        L[2] = 1; /* Lz does not contribute to the volume */
        break;
      case tri:
        e[0][0] = 1;
        e[1][0] = 0.5;  e[1][1] = sqrt(3)/2;
        L[2] = 1; /* Lz does not contribute to the volume */
        break;
    }
  }

  /*** Definition of the basis positions ***/
  for(i = 0; i < 8; i++)
    for(j = 0; j < 3; j++)
      basis[i][j] = 0; /* Clear all the basis positions */

  /* Non-zero basis positions */
  switch(type)
  {
    case bcc:
      cellsize = 2*separation/sqrt(3);
      nbasis = 2;
      basis[1][0] = basis[1][1] = basis[1][2] = 0.5;
      break;
    case fcc:
      cellsize = 2*separation/sqrt(2);
      nbasis = 4;
      basis[1][0] = basis[1][1] = 0.5;
      basis[2][0] = basis[2][2] = 0.5;
      basis[3][1] = basis[3][2] = 0.5;
      break;
    case dia:
      cellsize = 4*separation/sqrt(2);
      nbasis = 8;
      basis[1][0] = basis[1][1] = 0.5;
      basis[2][0] = basis[2][2] = 0.5;
      basis[3][1] = basis[3][2] = 0.5;
      basis[4][0] = basis[4][1] = basis[4][2] = 0.25;
      basis[5][0] = basis[5][1] = 0.75; basis[5][2] = 0.25;
      basis[6][0] = basis[6][2] = 0.75; basis[6][1] = 0.25;
      basis[7][1] = basis[7][2] = 0.75; basis[6][0] = 0.25;
      break;
    case hcp:
      cellsize = separation;
      nbasis = 2;
      basis[1][0] = .5; basis[1][1] = 0.25; basis[1][2] = sqrt(6.)/3;
      break;
    default:
      cellsize = separation;
      break;
  }

  /* Total number of unit cells */
  ncells = ceil(N/(1.f*nbasis));

  /* Box volume */
  V = L[0]*L[1]*L[2];

  /* Unit cells on the side of the box */
  if(fillbox) {
    nx = (int) ceil(L[0]/(cellsize*e[0][0]));
    ny = (int) ceil(L[1]/(cellsize*e[1][1]));
    ncells = nx*ny;
    N = ncells*nbasis;
  }
  else {
    nx = ceil(sqrt(ncells/V)*L[0]);
    ny = ceil(sqrt(ncells/V)*L[1]);
  }

  if(type == sq || type == tri) { /* 2D lattices */
    nz = 1;
  }
  else { /* 3D lattices */
    if(fillbox) {
      nz = (int) ceil(L[2]/(cellsize*e[2][2]));
      ncells *= nz;
      N = ncells*nbasis;
    }
    else {
      nz = ceil(pow(ncells/V,1/3.)*L[2]);
    }
  }

  /* Stretch factors */
  stretchfactor[0] = L[0]/(nx*e[0][0]);
  stretchfactor[1] = L[1]/(ny*e[1][1]);
  stretchfactor[2] = L[2]/(nz*e[2][2]);

  if(keepaspect) {
    minstretch = (stretchfactor[0] < stretchfactor[1])?
                     stretchfactor[0]:
                     stretchfactor[1];
    if(type != sq || type != tri)
      minstretch = (minstretch < stretchfactor[2])?
                     minstretch:
                     stretchfactor[2];
    for(i = 0; i < 3; i++) stretchfactor[i] = minstretch;
  }

  /* Output header */
  printf("#  Bravais lattice\n#\n"
         "# x \t\t y \t\t z \t\t r \t c #\n"
         "#---\t\t---\t\t---\t\t---\t---#\n");

  /* Positions in the Bravais lattice */
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++) {
      for(k = 0; k < nz; k++) {
        for(l = 0; l < nbasis; l++) { /* Loop over the elements in the basis */
          if(node < N) {
            for(d = 0; d < 3; d++) { /* Loop over dimensions */
              r[d] = -L[d]/2. + stretchfactor[d]*(i*e[0][d] + j*e[1][d]
                                                  + k*e[2][d] + basis[l][d]);
            }

            /* Wrap-around (triangular and hexagonal closed-packed lattices) */
            if(type == tri || type == hcp)
              if(r[0] > L[0]/2) r[0] -= L[0];

            /* No third dimension for 2D lattices */
            if(type == sq || type == tri)
              r[2] = 0;

            /* Output position */
            if(npbasis > 0) {
              for(n = 0; n < npbasis; n++) { /* Output user basis */
                printf("%f\t%f\t%f", r[0] + stretchfactor[0]*pbasis[n][0],
                                     r[1] + stretchfactor[1]*pbasis[n][1],
                                     r[2] + stretchfactor[2]*pbasis[n][2]);
                if(pbasis[n][3] > 0) printf("\t%f",
                                            stretchfactor[0]*pbasis[n][3]);
                if(pbasis[n][4] > 0) printf("\t%d", (int) pbasis[n][4]);
                printf("\n");
              }
            }
            else { /* Output lattice node */
              printf("%f\t%f\t%f", r[0], r[1], r[2]);
              if(radius >= 0) printf("\t%f", radius);
              if(colour >= 0) printf("\t%d", colour);
              printf("\n");
            }
            node++;
          }
        }
      }
    }
  }

  printf("\n#");
  return 0;
}
