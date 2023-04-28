#include "Error.h"
#include "GaussianElimination.h"
#include "KNearestNeighbors.h"
#include "Point.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern inline int ExplainError(int error, const char *hint);

static int Usage(const char *argv0, int error) {
  fprintf(
      stdout,
      "Usage: %s Lx Ly [OPTIONS]\n"
      "DESCRIPTION:\n"
      "\tUsing input data calculate values on a rectangle [0, Lx] x [0, Ly] with\n"
      "\tthe grid which consists Nx vertical and Ny horizontal lines\n"
      "\n"
      "OPTIONS:\n"
      "\tLx -- width of desired rectangle\n"
      "\tLy -- length of desired rectangle\n"
      "\t-h, --help\n"
      "\t\tProduce this help message\n"
      "\t--Nx <unsigned> ( = 100 )\n"
      "\t\tAmount of splitting of a desired rectangle along X axis\n"
      "\t--Ny <unsigned> ( = Nx )\n"
      "\t\tAmount of splitting of a desired rectangle along Y axis\n"
      "\t-k <unsigned>, --k_neighbors <unsigned> ( = 0 )\n"
      "\t\tIf specified, the value of the interpolant at each evaluation\n"
      "\t\tpoint will be computed using only this many nearest data points.\n"
      "\t\tAll the data points are used by default.\n"
      "\n"
      "INPUT FORMAT:\n"
      "\t# N\n"
      "\tx1 y1 f(x1, y1)\n"
      "\tx2 y2 f(x2, y2)\n"
      "\t . . .\n"
      "\txN yN f(xN, yN)\n"
      "\n"
      "OUTPUT FORMAT:\n"
      "\tConsists (Nx + 1) x (Ny + 1) rows of three columns: x, y, approximate value in (x, y)\n",
      argv0);
  return error;
}

static inline double RadialBasisFunction(double distance) {
  // `thin_plate_spline` is used as it does not require shape parameter
  return distance * log(pow(distance, distance));
}

double Distance(Point l, Point r) { return sqrt(pow(r.x - l.x, 2.) + pow(r.y - l.y, 2.)); }

static inline void PrintPoints(FILE *to, const Point *points, int N) {
  for (int i = 0; i < N; i++) {
    fprintf(to, "%20.15lf %20.15lf\n", points[i].x, points[i].y);
  }
}

static inline void PrintDoubleArray(FILE *to, const double *values, int N) {
  for (int i = 0; i < N; i++) {
    fprintf(to, "%20.15lf\n", values[i]);
  }
}

static inline void PrintGrid(FILE *to, int Nx, int Ny, const Point *grid) {
  for (int i = 0; i <= Nx; i++) {
    for (int j = 0; j <= Ny; j++) {
      fprintf(to, "( %5.5lf, %5.5lf )%c", grid[i * (Nx + 1) + j].x, grid[i * (Nx + 1) + j].y,
              " \n"[j == Ny]);
    }
  }
  fprintf(to, "\n");
}

static inline void PrintResult(FILE *to, int N, const Point *points, const double *values) {
  for (int i = 0; i < N; i++) {
    fprintf(to, "%20.15lf %20.15lf %20.15lf\n", points[i].x, points[i].y, values[i]);
  }
}

static Point *ConstructGrid(int Nx, int Ny, double Lx, double Ly) {
  Point *grid;
  const double x_length = Lx / Nx;
  const double y_length = Ly / Ny;
  grid = malloc((Nx + 1) * (Ny + 1) * sizeof(*grid));
  if (!grid)
    return grid;

  for (int i = 0; i <= Nx; ++i) {
    for (int j = 0; j <= Ny; ++j) {
      grid[i * (Nx + 1) + j].x = i * x_length;
      grid[i * (Nx + 1) + j].y = j * y_length;
    }
  }
  return grid;
}

static int ReadInputData(FILE *from, int N, Point *points, double *values) {
  for (int i = 0u; i < N; i++) {
    if (fscanf(from, "%lf%lf%lf", &points[i].x, &points[i].y, values + i) != 3) {
      return InputError;
    }
  }
  return Success;
}

static int ComputeWeights(const Point *points, const double *values, int N, double *matrix,
                          double *weights) {
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      matrix[i * N + j] = RadialBasisFunction(Distance(points[i], points[j]));
    }
  }

  return GaussMaxCol(matrix, N, values, weights);
}

static int SpatialInterpolation(int input_N, const Point *input_points, const double *input_values,
                                int desired_N, const Point *desired_points, double *result_values) {
  int error = Success;
  double value;
  double *weights = malloc(input_N * sizeof(*weights));
  double *matrix = malloc(input_N * input_N * sizeof(*matrix));
  if (!matrix || !weights) {
    error = ExplainError(NotEnoughMemory, "spatial interpolation");
    goto cleanup_spatial;
  }

  if (ComputeWeights(input_points, input_values, input_N, matrix, weights) != 0) {
    error = ExplainError(LogicError, "matrix is singular");
    goto cleanup_spatial;
  }

  for (int i = 0; i < desired_N; ++i) {
    value = 0;
    for (int k = 0; k < input_N; ++k) {
      value += weights[k] * RadialBasisFunction(Distance(desired_points[i], input_points[k]));
    }
    result_values[i] = value;
  }

cleanup_spatial:
  free(weights);
  free(matrix);
  return error;
}

static int SpatialInterpolationUsingNearestNeighbors(int input_N, const Point *input_points,
                                                     const double *input_values, int k_neighbors,
                                                     int desired_N, const Point *desired_points,
                                                     double *result_values) {
  Point current;
  double value;
  int error = Success;
  double *weights = malloc(k_neighbors * sizeof(*weights));
  double *matrix = malloc(k_neighbors * k_neighbors * sizeof(*matrix));
  Point *neighbors = malloc(k_neighbors * sizeof(*neighbors));
  double *neighbors_values = malloc(k_neighbors * sizeof(*neighbors_values));
  DistanceAndIndex *distances = malloc(input_N * sizeof(*distances));
  if (!matrix || !weights || !neighbors || !distances || !neighbors_values) {
    error = ExplainError(NotEnoughMemory, "spatial interpolation with neighbors");
    goto cleanup_neighbors;
  }

  for (int i = 0; i < desired_N; ++i) {
    current = desired_points[i];
    KNearestNeighbors(input_points, input_N, current, k_neighbors, distances);
    for (int j = 0; j < k_neighbors; ++j) {
      neighbors[j] = input_points[distances[j].index];
      neighbors_values[j] = input_values[distances[j].index];
    }
    if (ComputeWeights(neighbors, neighbors_values, k_neighbors, matrix, weights) != 0) {
      error = ExplainError(LogicError, "RBF matrix is singular");
      goto cleanup_neighbors;
    }
    value = 0;
    for (int k = 0; k < k_neighbors; ++k) {
      value += weights[k] * RadialBasisFunction(Distance(current, neighbors[k]));
    }
    result_values[i] = value;
  }

cleanup_neighbors:
  free(weights);
  free(matrix);
  free(neighbors);
  free(distances);
  free(neighbors_values);
  return error;
}

int main(int argc, char **argv) {
  double Lx = 0, Ly = 0;
  int Nx = 100, Ny = 0;
  int N = 0;
  int error = Success;
  Point *input_points = NULL;
  double *input_values = NULL;
  Point *grid = NULL;
  double *grid_values = NULL;
  int k_neighbors = 0;
  struct option long_options[] = {{"help", no_argument, 0, 'h'},
                                  {"Nx", required_argument, 0, 0},
                                  {"Ny", required_argument, 0, 0},
                                  {"k_neighbors", required_argument, 0, 'k'},
                                  {0, 0, 0, 0}};
  int option_index = 0, c = 0;
  int is_Ny_specified = 0;

  while (1) {
    c = getopt_long(argc, argv, "hk:", long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      if (strcmp("Nx", long_options[option_index].name) == 0) {
        if (!(sscanf(optarg, "%d", &Nx) == 1 && Nx >= 0))
          return ExplainError(IncorrectUsage, "Nx expected to be positive int");
      } else if (strcmp("Ny", long_options[option_index].name) == 0) {
        if (!(sscanf(optarg, "%d", &Ny) == 1 && Ny >= 0))
          return ExplainError(IncorrectUsage, "Ny expected to be positive int");
        is_Ny_specified = 1;
      }
      break;
    case 'k':
      if (!(sscanf(optarg, "%d", &k_neighbors) == 1 && k_neighbors >= 0))
        return ExplainError(IncorrectUsage, "neighbours expected to be positve int");
      break;
    case 'h':
      return Usage(argv[0], Success);
    default:
      return ExplainError(IncorrectUsage, "unexpected paramer");
    }
  }

  if (optind + 2 > argc) {
    return ExplainError(IncorrectUsage, "expected Lx and Ly");
  }

  if (!(sscanf(argv[optind++], "%lf", &Lx) == 1 && sscanf(argv[optind++], "%lf", &Ly) == 1)) {
    return ExplainError(InputError, "length or width");
  }

  if (!is_Ny_specified)
    Ny = Nx;

  if (fscanf(stdin, "# %d", &N) != 1)
    return ExplainError(InputError, "N");

  if (N < 2)
    return ExplainError(LogicError, "N should be greater than 1");

  input_points = malloc(N * sizeof(*input_points));
  input_values = malloc(N * sizeof(*input_values));
  if (!input_points || !input_values) {
    error = ExplainError(NotEnoughMemory, "input");
    goto cleanup_main;
  }

  if (ReadInputData(stdin, N, input_points, input_values) != Success) {
    error = ExplainError(InputError, "data");
    goto cleanup_main;
  }

  grid = ConstructGrid(Nx, Ny, Lx, Ly);
  grid_values = malloc((Nx + 1) * (Ny + 1) * sizeof(*grid_values));
  if (!grid || !grid_values) {
    error = ExplainError(NotEnoughMemory, "grid");
    goto cleanup_main;
  }

  if (k_neighbors != 0) {
    k_neighbors = k_neighbors > N ? N : k_neighbors;
    SpatialInterpolationUsingNearestNeighbors(N, input_points, input_values, k_neighbors,
                                              (Nx + 1) * (Ny + 1), grid, grid_values);
  } else
    SpatialInterpolation(N, input_points, input_values, (Nx + 1) * (Ny + 1), grid, grid_values);

  PrintResult(stdout, (Nx + 1) * (Ny + 1), grid, grid_values);

cleanup_main:
  free(input_points);
  free(input_values);
  free(grid);
  free(grid_values);
  return error;
}
