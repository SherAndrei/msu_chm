#include "Error.h"
#include "GaussianElimination.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern inline int ExplainError(int error, const char *hint);

static int Usage(const char *argv0, int error) {
  fprintf(
      stdout,
      "Usage: %s Lx Ly [Nx] [Ny]\n"
      "DESCRIPTION:\n"
      "\tUsing input data calculate values on a rectangle [0, Lx] x [0, Ly] with\n"
      "\tthe grid which consists Nx vertical and Ny horizontal lines\n"
      "\n"
      "OPTIONS:\n"
      "\tdouble Lx -- width of desired rectangle\n"
      "\tdouble Ly -- length of desired rectangle\n"
      "\tunsigned Nx (default = 100) -- amount of splitting of a desired rectangle along X axis\n"
      "\tunsigned Ny (default = Nx) -- amount of splitting of a desired rectangle along Y axis\n"
      "\n"
      "INPUT FORMAT:\n"
      "\tN\n"
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

typedef struct PointType {
  double x;
  double y;
} Point;

static inline double RadialBasisFunction(double distance) {
  // `thin_plate_spline` is used as it does not require shape parameter
  return distance * log(pow(distance, distance));
}

static inline double Distance(Point l, Point r) {
  return sqrt(pow(r.x - l.x, 2.) + pow(r.y - l.y, 2.));
}

static inline void PrintPoints(FILE *to, const Point *points, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    fprintf(to, "%20.15lf %20.15lf\n", points[i].x, points[i].y);
  }
}

static inline void PrintDoubleArray(FILE *to, const double *values, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    fprintf(to, "%20.15lf\n", values[i]);
  }
}

static inline void PrintGrid(FILE *to, unsigned Nx, unsigned Ny, const Point *grid) {
  for (unsigned i = 0; i <= Nx; i++) {
    for (unsigned j = 0; j <= Ny; j++) {
      fprintf(to, "( %5.5lf, %5.5lf )%c", grid[i * (Nx + 1) + j].x, grid[i * (Nx + 1) + j].y,
              " \n"[j == Ny]);
    }
  }
  fprintf(to, "\n");
}

static inline void PrintResult(FILE *to, unsigned N, const Point *points, const double *values) {
  for (unsigned i = 0; i < N; i++) {
    fprintf(to, "%20.15lf %20.15lf %20.15lf\n", points[i].x, points[i].y, values[i]);
  }
}

static Point *ConstructGrid(unsigned Nx, unsigned Ny, double Lx, double Ly) {
  Point *grid;
  const double x_length = Lx / Nx;
  const double y_length = Ly / Ny;
  grid = malloc((Nx + 1) * (Ny + 1) * sizeof(*grid));
  if (!grid)
    return grid;

  for (unsigned i = 0; i <= Nx; ++i) {
    for (unsigned j = 0; j <= Ny; ++j) {
      grid[i * (Nx + 1) + j].x = i * x_length;
      grid[i * (Nx + 1) + j].y = j * y_length;
    }
  }
  return grid;
}

static int ReadInputData(FILE *from, unsigned N, Point *points, double *values) {
  for (unsigned i = 0u; i < N; i++) {
    if (fscanf(from, "%lf%lf%lf", &points[i].x, &points[i].y, values + i) != 3) {
      return InputError;
    }
  }
  return Success;
}

static int SpatialInterpolation(unsigned input_N, const Point *input_points,
                                const double *input_values, unsigned desired_N,
                                const Point *desired_points, double *result_values) {
  double value;
  double *weights = malloc(input_N * sizeof(*weights));
  double *matrix = malloc(input_N * input_N * sizeof(*matrix));
  if (!matrix || !weights) {
    free(weights);
    free(matrix);
    return ExplainError(NotEnoughMemory, "spatial interpolation");
  }

  for (unsigned i = 0; i < input_N; ++i) {
    for (unsigned j = 0; j < input_N; ++j) {
      matrix[i * input_N + j] = RadialBasisFunction(Distance(input_points[i], input_points[j]));
    }
  }

  if (GaussMaxCol(matrix, input_N, input_values, weights) != 0) {
    free(weights);
    free(matrix);
    return ExplainError(LogicError, "matrix is singular");
  }

  for (unsigned i = 0; i < desired_N; ++i) {
    value = 0;
    for (unsigned k = 0; k < input_N; ++k) {
      value += weights[k] * RadialBasisFunction(Distance(desired_points[i], input_points[k]));
    }
    result_values[i] = value;
  }

  free(weights);
  free(matrix);
  return Success;
}

// TODO: add neighbors int, optional
//       If specified, the value of the interpolant at each
//       evaluation point will be computed using only this
//       many nearest data points. All the data points are used by default.
// TODO: add smoothing float or (P,) array_like, optional
//       Smoothing parameter. The interpolant perfectly fits the data
//       when this is set to 0. For large values, the interpolant approaches
//       a least squares fit of a polynomial with the specified degree. Default is 0.
// TODO: add getopt to distinct parameters

int main(int argc, char **argv) {
  double Lx, Ly;
  unsigned Nx, Ny;
  unsigned N;
  Point *input_points;
  double *input_values;
  Point *grid;
  double *grid_values;

  if (!(argc >= 3 && argc <= 5))
    return Usage(argv[0], IncorrectUsage);

  if (!(sscanf(argv[1], "%lf", &Lx) == 1 && sscanf(argv[2], "%lf", &Ly) == 1)) {
    return ExplainError(InputError, "length or width");
  }

  if (argc >= 4) {
    if (sscanf(argv[3], "%u", &Nx) != 1) {
      return ExplainError(InputError, "Nx");
    }
  } else {
    Nx = 100;
  }

  if (argc == 5) {
    if (sscanf(argv[4], "%u", &Ny) != 1) {
      return ExplainError(InputError, "Ny");
    }
  } else {
    Ny = Nx;
  }

  if (fscanf(stdin, "%u", &N) != 1)
    return ExplainError(InputError, "N");

  if (N < 2u)
    return ExplainError(LogicError, "N should be greater than 1");

  input_points = malloc(N * sizeof(*input_points));
  input_values = malloc(N * sizeof(*input_values));
  if (!input_points || !input_values) {
    free(input_points);
    free(input_values);
    return ExplainError(NotEnoughMemory, "input");
  }

  if (ReadInputData(stdin, N, input_points, input_values) != Success) {
    free(input_points);
    free(input_values);
    return ExplainError(InputError, "data");
  }

  grid = ConstructGrid(Nx, Ny, Lx, Ly);
  grid_values = malloc((Nx + 1) * (Ny + 1) * sizeof(*grid_values));
  if (!grid || !grid_values) {
    free(input_points);
    free(input_values);
    free(grid);
    free(grid_values);
    return ExplainError(NotEnoughMemory, "grid");
  }

  SpatialInterpolation(N, input_points, input_values, (Nx + 1) * (Ny + 1), grid, grid_values);

  PrintResult(stdout, (Nx + 1) * (Ny + 1), grid, grid_values);

  free(input_points);
  free(input_values);
  free(grid);
  free(grid_values);
}
