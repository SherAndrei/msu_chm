#include "Error.h"

#include <math.h>
#include <stdio.h>

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

static double IntegrateTriangle(Point A, Point B, Point C, double (*f)(double, double)) {
  return (1. / 3.) * f(A.x1, A.x2) + f(B.x1, B.x2) + f(C.x1, C.x2);
}

typedef struct TriangleType {
  Point A;
  Point B;
  Point C;
} Triangle;

static double Integrate(const Triangle *restrict triangles, unsigned amount_of_trinagles) {
  double result = 0;
  for (unsigned i = 0; i < amount_of_trinagles; ++i) {
  }
  return result;
}

static double Polynom(double x1, double x2) { return pow(x1, 4) + pow(x1 * x2, 2) + pow(x2, 4); }

typedef struct TriangleVertexNumbersType {
  unsigned A_pos;
  unsigned B_pos;
  unsigned C_pos;
} TriangleVertexNumbers;

static int FscanfTriangleVertexNumbers(FILE* in, unsigned amount_of_triangles, TriangleVertexNumbers* vertex_numbers)
{
  TriangleVertexNumbers* current;
  for (unsigned i = 0; i < amount_of_triangles; ++i) {
    current = vertex_numbers + i;
    if (fscanf(in, "%*u%u%u%u", &current->A_pos, &current->B_pos, &current->C_pos) != 4)
      return 1;
  }
  return 0;
}

static int FscanfVertvies(FILE* in, unsigned amount_of_vertices, const TriangleVertexNumbers* vertex_numbers, Point* vertices)
{
  Point* current;
  for (unsigned i = 0; i < amount_of_vertices; ++i) {
    current = vertices + i;
    if (fscanf(in, "%*u%lf%lf", &current->x1, &current->x2) != 3)
      return 1;
  }
  return 0;
}

int main(int argc, const char *argv[]) {
  double result = 0;
  unsigned amount_of_vertices;
  unsigned amount_of_triangles;
  TriangleVertexNumbers *vertex_numbers = NULL;
  Point* vertices = NULL;
  FILE *in = stdin;

  if (argc != 1)
    return Usage(argv[0], IncorrectUsage);

  if (fscanf(in, "%u", &amount_of_vertices) != 1)
    return Usage(argv[0], InputError, "amount of vertices");

  if (fscanf(in, "%u", &amount_of_triangles) != 1)
    return Usage(argv[0], InputError, "amount of triangles");

  if (fscanf(in, "%*u%*u") != 2)
    return Usage(argv[0], InputError, "amount of inner and outer edges");

  vertices = (Point*) malloc(amount_of_vertices * sizeof(Point));
  vertex_numbers =
      (TriangleVertexNumbers *)malloc(amount_of_triangles * sizeof(TriangleVertexNumbers));
  if (!vertex_numbers || !vertices) {
    free(vertex_numbers);
    free(vertices);
    fclose(in);
    return Usage(argv[0], NotEnoughMemory, "");
  }

  if (FscanfVertices(in, amount_of_vertices, vertices) != 0)
    return Usage(argv[0], InputError, "vertices");

  // TODO: add check for correct vertex number
  if (FscanfTriangleVertexNumbers(in, amount_of_triangles, vertex_numbers) != 0)
    return Usage(argv[0], InputError, "vertex numbers");

  result = Integrate(vertices, vertex_numbers, amount_of_triangles, Polynom);
  printf("%.14lf\n", result);

  free(vertex_numbers);
  free(vertices);
  fclose(in);
}
