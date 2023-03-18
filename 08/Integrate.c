#include "Error.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int Usage(const char *argv0, int error, const char *hint) {
  switch (error) {
  case IncorrectUsage:
    fprintf(stdout,
            "Usage: %s\n"
            "DESCRIPTION:\n"
            "\tCalculate the approximate value of the 2-dimensional definite integral\n"
            "\n"
            "EXPECTED INPUT FORMAT:\n"
            "\t<amount of vertices>\n"
            "\t<amount of triangles>\n"
            "\t<amount of inner edges>\n"
            "\t<amount of outer edges>\n"
            "(for each vertex)\n"
            "\t<vertex number>:<x y> (vertex coordinates)\n"
            "\t. . .\n"
            "(for each triangle)\n"
            "\t<triangle number>:<i j k> (vertex numbers)\n"
            "\t. . .\n"
            "(for each inner edge)\n"
            "\t<inner edge number>:<m n> (vertex numbers)\n"
            "\t. . .\n"
            "(for each outer edge)\n"
            "\t<outer edge number>:<m n> (vertex numbers)\n"
            "\t. . .\n",
            argv0);
    break;
  case InputError:
    fprintf(stderr, "error: incorrect input, cannot parse %s\n", hint);
    break;
  case NotEnoughMemory:
    fprintf(stderr, "error: cannot allocate memory\n");
    break;
  case Success:
    break;
  }
  return error;
}

typedef struct PointType {
  double x1;
  double x2;
} Point;

static inline Point Middle(Point l, Point r) {
  Point res;
  res.x1 = (l.x1 + r.x1) / 2.;
  res.x2 = (l.x2 + r.x2) / 2.;
  return res;
}

static inline double TriangleArea(Point A, Point B, Point C) {
  // using determinant method
  return (1. / 2.) * fabs((A.x1 * (B.x2 - C.x2) + B.x1 * (C.x2 - A.x2) + C.x1 * (A.x2 - B.x2)));
}

static inline double TriangleQuadrature(Point A, Point B, Point C, double (*f)(Point)) {
  return (TriangleArea(A, B, C) / 3.) * (f(Middle(A, B)) + f(Middle(B, C)) + f(Middle(C, A)));
}

typedef struct TriangleVertexNumbersType {
  unsigned A_pos;
  unsigned B_pos;
  unsigned C_pos;
} TriangleVertexNumbers;

static double Integrate(const Point *vertices, const TriangleVertexNumbers *vertex_numbers,
                        unsigned n_triangles, double (*f)(Point)) {
  double result = 0;
  Point A, B, C;
  for (unsigned i = 0; i < n_triangles; ++i) {
    A = vertices[vertex_numbers[i].A_pos];
    B = vertices[vertex_numbers[i].B_pos];
    C = vertices[vertex_numbers[i].C_pos];
    result += TriangleQuadrature(A, B, C, f);
  }
  return result;
}

static double Polynom(Point p) { return pow(p.x1, 4) + pow(p.x1 * p.x2, 2) + pow(p.x2, 4); }

static int FscanfTriangleVertexNumbers(FILE *in, unsigned n_vertices, unsigned n_triangles,
                                       TriangleVertexNumbers *vertex_numbers) {
  TriangleVertexNumbers *current;
  unsigned dummy;
  (void)n_vertices;
  for (unsigned i = 0; i < n_triangles; ++i) {
    current = vertex_numbers + i;
    if (fscanf(in, "%u%u%u%u", &dummy, &current->A_pos, &current->B_pos, &current->C_pos) != 4)
      return 1;
    if (--current->A_pos >= n_vertices || --current->B_pos >= n_vertices ||
        --current->C_pos >= n_vertices)
      return 1;
  }
  return 0;
}

static int FscanfVertices(FILE *in, unsigned n_vertices, Point *vertices) {
  Point *current;
  unsigned dummy;
  for (unsigned i = 0; i < n_vertices; ++i) {
    current = vertices + i;
    if (fscanf(in, "%u%lf%lf", &dummy, &current->x1, &current->x2) != 3)
      return 1;
  }
  return 0;
}

int main(int argc, const char *argv[]) {
  double result = 0;
  unsigned dummy;
  unsigned n_vertices;
  unsigned n_triangles;
  TriangleVertexNumbers *vertex_numbers = NULL;
  Point *vertices = NULL;
  FILE *in = stdin;

  if (argc != 1)
    return Usage(argv[0], IncorrectUsage, "");

  if (fscanf(in, "%u", &n_vertices) != 1)
    return Usage(argv[0], InputError, "amount of vertices");

  if (fscanf(in, "%u", &n_triangles) != 1)
    return Usage(argv[0], InputError, "amount of triangles");

  if (fscanf(in, "%u%u", &dummy, &dummy) != 2)
    return Usage(argv[0], InputError, "amount of inner and outer edges");

  vertices = malloc(n_vertices * sizeof(*vertices));
  vertex_numbers = malloc(n_triangles * sizeof(*vertex_numbers));
  if (!vertex_numbers || !vertices) {
    free(vertex_numbers);
    free(vertices);
    fclose(in);
    return Usage(argv[0], NotEnoughMemory, "");
  }

  if (FscanfVertices(in, n_vertices, vertices) != 0) {
    free(vertex_numbers);
    free(vertices);
    fclose(in);
    return Usage(argv[0], InputError, "vertices");
  }

  if (FscanfTriangleVertexNumbers(in, n_vertices, n_triangles, vertex_numbers) != 0) {
    free(vertex_numbers);
    free(vertices);
    fclose(in);
    return Usage(argv[0], InputError, "vertex numbers");
  }

  result = Integrate(vertices, vertex_numbers, n_triangles, Polynom);
  printf("%.14lf\n", result);

  free(vertex_numbers);
  free(vertices);
  fclose(in);
}
