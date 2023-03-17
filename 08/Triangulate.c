#include "Error.h"

#include <stdio.h>

static int Usage(const char *argv0, int error) {
  fprintf(
      stdout,
      "Usage: %s Lx Ly [Nx] [Ny]\n"
      "DESCRIPTION:\n"
      "\tFor a given rectangle with sides Lx, Ly construct its triangulation\n"
      "\n"
      "OPTIONS:\n"
      "\tdouble Lx -- width of desired rectangle\n"
      "\tdouble Ly -- length of desired rectangle\n"
      "\tunsigned Nx (default = 100) -- amount of splitting of a desired rectangle along X axis\n"
      "\tunsigned Ny (default = Nx) -- amount of splitting of a desired rectangle along Y axis\n",
      argv0);
  return error;
}

static inline void PrintVertexNumbersAndTheirCoordinates(double Lx, double Ly, unsigned Nx, unsigned Ny, FILE *out) {
  const double x_length = Lx / Nx;
  const double y_length = Ly / Ny;
  unsigned vertex_number = 0;
  for (unsigned i = 0; i <= Nx; ++i)
    for (unsigned j = 0; j <= Ny; ++j)
      fprintf(out, "%u %lf %lf\n", ++vertex_number, i * x_length, j * y_length);
}

static inline void PrintTriangleNumberToVertexNumbers(unsigned Nx, unsigned Ny, FILE *out) {
  unsigned triangle_number = 0;
  // for each square print two triangles
  for (unsigned i = 1; i <= Nx; ++i) {
    for (unsigned j = 1; j <= Ny; ++j) {
      fprintf(out, "%u %u %u %u\n", ++triangle_number,
        (Ny + 1) * (i - 1) + j,
        (Ny + 1) * (i - 1) + j + 1,
        (Ny + 1) * i + j);
      fprintf(out, "%u %u %u %u\n", ++triangle_number,
        (Ny + 1) * (i - 1) + j + 1,
        (Ny + 1) * i + j,
        (Ny + 1) * i + j + 1);
    }
  }
}

static inline void PrintInnerEdgeToVertexNumbers(unsigned Nx, unsigned Ny, FILE *out) {
  unsigned inner_edge_number = 0;
  // for each square print its diagonal, bottom and right edge
  // except those on the outer edge of a rectangle
  for (unsigned i = 1; i <= Nx; ++i) {
    for (unsigned j = 1; j <= Ny; ++j) {
      fprintf(out, "%u %u %u\n", ++inner_edge_number,
        (Ny + 1) * (i - 1) + j + 1,
        (Ny + 1) * i + j);
      if (i != Nx) {
        fprintf(out, "%u %u %u\n", ++inner_edge_number,
          Ny * i + j + 1,
          Ny * i + j);
      }
      if (j != Ny) {
        fprintf(out, "%u %u %u\n", ++inner_edge_number,
          Ny * (i - 1) + j,
          Ny * i + j);
      }
    }
  }
}

static inline void PrintOuterEdgeToVertexNumbers(unsigned Nx, unsigned Ny, FILE *out) {
  unsigned outer_edge_number = 0;
  for (unsigned i = 0; i < Nx; ++i) {
    fprintf(out, "%u %u %u\n", ++outer_edge_number,
      (Ny + 1) * i + 1,
      (Ny + 1) * (i + 1) + 1);
    fprintf(out, "%u %u %u\n", ++outer_edge_number,
      (Ny + 1) * i + Ny + 1,
      (Ny + 1) * (i + 1) + Ny + 1);
  }
  for (unsigned j = 0; j < Ny; ++j) {
    fprintf(out, "%u %u %u\n", ++outer_edge_number,
      j + 1,
      j + 2);
    fprintf(out, "%u %u %u\n", ++outer_edge_number,
      (Ny + 1) * Nx + j + 1,
      (Ny + 1) * Nx + j + 2);
  }
}

static void PrintTriangulation(double Lx, double Ly, unsigned Nx, unsigned Ny, FILE *out) {
  // amount of vertices
  fprintf(out, "%u\n", (Nx + 1) * (Ny + 1));
  // amount of triangles
  fprintf(out, "%u\n", 2 * Nx * Ny);
  // amount of inner edges
  fprintf(out, "%u\n", 3 * Nx * Ny - Ny - Nx);
  // amount of outer edges
  fprintf(out, "%u\n", 2 * Nx + 2 * Ny);

  PrintVertexNumbersAndTheirCoordinates(Lx, Ly, Nx, Ny, out);
  PrintTriangleNumberToVertexNumbers(Nx, Ny, out);
  PrintInnerEdgeToVertexNumbers(Nx, Ny, out);
  PrintOuterEdgeToVertexNumbers(Nx, Ny, out);
}

int main(int argc, const char *argv[]) {
  double Lx, Ly;
  unsigned Nx, Ny;

  if (!(argc >= 3 && argc <= 5))
    return Usage(argv[0], IncorrectUsage);

  if (!(sscanf(argv[1], "%lf", &Lx) == 1 && sscanf(argv[2], "%lf", &Ly) == 1)) {
    fprintf(stderr, "error: parsing input parameters failed\n");
    return InputError;
  }

  if (argc >= 4) {
    if (sscanf(argv[3], "%u", &Nx) != 1) {
      fprintf(stderr, "error: parsing Nx failed\n");
      return InputError;
    }
  } else {
    Nx = 100;
  }

  if (argc == 5) {
    if (sscanf(argv[4], "%u", &Ny) != 1) {
      fprintf(stderr, "error: parsing Ny failed\n");
      return InputError;
    }
  } else {
    Ny = Nx;
  }

  PrintTriangulation(Lx, Ly, Nx, Ny, stdout);
}
