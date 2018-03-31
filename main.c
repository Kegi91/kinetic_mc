#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "rng.h"

#define pi 3.14159265358979323846
#define k_b 8.6173303e-5 // In electron volts

typedef struct Defect {
  double x, y, z;
  int t; // Defect type 1 == interstitial, -1 == vacancy, 0 == deleted
} dfct;

dfct *gen_randpos(float sigma1, float sigma2, int size);
double dist(dfct v1, dfct v2);
int check_defects(dfct *defects, double r_rec, int size);
int check_single_defect(dfct *x, double r_rec, int idx, int size);
void update_cum_func(dfct *defects, double *R, double rate_i, double rate_v, int size);
int find_dfct(double *R, double value, int size);
int find_dfct_binary(double *R, double value, int size);
void single_simulation(
  double T, double r_rec, int size,
  double t_lim, int write_freq, int rank
);

int main(int argc, char **argv) {
  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  seed_mt(time(NULL)+rank);

  if (argc < 5) {
    fprintf(stderr, "\n%s\n\n", "Wrong cmd line parameters");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  double T = atof(argv[1]);
  double r_rec = atof(argv[2]);
  int size = atoi(argv[3]);
  double t_lim = atof(argv[4]);

  single_simulation(T, r_rec, size, t_lim, 1e4, rank);

  MPI_Finalize();

  return 0;
}

void single_simulation(
  double T, double r_rec, int size, double t_lim, int write_freq, int rank
) {

  double t = 0, rate_i, rate_v, *R, u;
  double theta, phi;
  int N = size, i, rm = -1;
  long iteration = 0, i_jump = 0, v_jump = 0;

  char fname[20];
  sprintf(fname, "output/%03d.dat", rank);

  FILE *f_out;
  f_out = fopen(fname, "w");

  R = malloc(N*sizeof(double));
  if (R == NULL) {
    fprintf(stderr, "\n%s\n\n", "Error allocating memory");
    exit(1);
  }

  dfct *x = gen_randpos(80, 40, N);

  rate_i = 1.717 * exp(-1.37/(k_b*T));
  rate_v = 0.001282 * exp(-0.1/(k_b*T));

  N -= check_defects(x, r_rec, size);

  while (t < t_lim) {
    // Profiling showed update_cum_func to use
    // the most cpu time so let's call it only
    // when necessary by checking if last iteration
    // removed any defects
    if (rm) update_cum_func(x, R, rate_i, rate_v, size);

    // Select random defect randomly
    // weighted by cum_func R
    u = rand_mt();
    i = find_dfct_binary(R, u*R[size-1], size);

    // if (x[i].t == 0) {
    //   fprintf(stderr, "ERROR\n");
    // }

    if (x[i].t == 1) {
      i_jump++;
    } else {
      v_jump++;
    }

    // // Checking that find_dfct_binary works
    // if (i != find_dfct(R, u*R[size-1], size)) {
    //   fprintf(stderr, "\n%s\n\n", "ERROR");
    //   exit(1);
    // }

    // Generate random spherical
    // coordinate angles phi and theta
    theta = acos(1-2*rand_mt());
    phi = 2*pi*rand_mt();

    // Move defect i by 2.35å in the
    // randomly generated direction
    x[i].x += 2.35*sin(theta)*cos(phi);
    x[i].y += 2.35*sin(theta)*sin(phi);
    x[i].z += 2.35*cos(theta);

    // Check if there's overlap with moved defect
    // and any of the other type defects and fix it
    rm = check_single_defect(x, r_rec, i, size);
    N -= rm;

    t += (-log(rand_mt())/R[size-1]);
    iteration++;

    if (iteration % write_freq == 0) {
      fprintf(f_out, "%12.1lf%6d%11ld%11ld\n", t, N, i_jump, v_jump);
    }
  }

  free(x);
  free(R);

  fclose(f_out);
}

// Generating normal distributed
// random positions using Box-Muller generator
dfct *gen_randpos(float sigma1, float sigma2, int size) {
  double u1, u2, phi, r, x, y, z, sigma;
  dfct *array = malloc(size*sizeof(dfct));

  if (array == NULL) {
    fprintf(stderr, "\n%s\n\n", "Error allocating memory");
    exit(1);
  }

  for (int i = 0; i<size; i++) {
    u1 = rand_mt();
    u2 = rand_mt();

    phi = 2*pi*u1;
    r = sqrt(-2*log(u2));

    if (i < size/2) {
      sigma = sigma1;
    } else {
      sigma = sigma2;
    }

    x = r*cos(phi) * sigma;
    y = r*sin(phi) * sigma;

    // z \in [0,70] å
    z = rand_mt()*70;

    array[i].x = x;
    array[i].y = y;
    array[i].z = z;

    if (i < size/2) {
      array[i].t = 1;
    } else {
      array[i].t = -1;
    }
  }

  return array;
}

double dist(dfct v1, dfct v2) {
  double x = v1.x-v2.x;
  double y = v1.y-v2.y;
  double z = v1.z-v2.z;

  return sqrt(x*x + y*y + z*z);
}

int check_defects(dfct *defects, double r_rec, int size) {
  int removed = 0, t1, t2;
  double distance;

  // All interstitials
  for (int i = 0; i<size/2; i++) {
    // All vacancies
    for (int v = size/2; v<size; v++) {
      t1 = defects[i].t;
      t2 = defects[v].t;
      distance = dist(defects[i], defects[v]);

      if (t1 != 0 && t2 != 0 && distance < r_rec) {
        defects[i].t = 0;
        defects[v].t = 0;
        removed += 2;
      }
    }
  }

  return removed;
}

int check_single_defect(dfct *x, double r_rec, int idx, int size) {
  int lower_lim, upper_lim;

  if (idx < size/2) {
    lower_lim = size/2;
    upper_lim = size;
  } else {
    lower_lim = 0;
    upper_lim = size/2;
  }

  for (int i = lower_lim; i<upper_lim; i++) {
    if (x[i].t != 0 && dist(x[i], x[idx]) < r_rec) {
      x[i].t = 0;
      x[idx].t = 0;

      return 2;
    }
  }

  return 0;
}

void update_cum_func(
  dfct *defects, double *R, double rate_i, double rate_v, int size
) {
  if (defects[0].t == 1) {
    R[0] = rate_i;
  } else if (defects[0].t == -1) {
    R[0] = rate_v;
  } else {
    R[0] = 0;
  }

  for (int i = 1; i<size; i++) {
    if (defects[i].t == 1) {
      R[i] = R[i-1] + rate_i;
    } else if (defects[i].t == -1) {
      R[i] = R[i-1] + rate_v;
    } else {
      R[i] = R[i-1];
    }
  }
}

int find_dfct(double *R, double value, int size) {

  // Simply iterate over all elements of
  // R to find_dfct to find the match
  for (int i = 0; i<size; i++) {
    if (value < R[i]) return i;
  }

  return size-1;
}

int find_dfct_binary(double *R, double value, int size) {
  if (R[0] > value) {
    return 0;
  } else if (R[size-2]<value) {
    return size-1;
  }

  int l = 0, r = size-1, c;

  // Binary search over R to find i
  // such that R[i-1] < value <= R[i+1]
  while (r-l > 1) {
    c = (l+r)/2;

    if (R[c] > value) {
      r = c;
    } else {
      l = c;
    }
  }

  return r;
}
