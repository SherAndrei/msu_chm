#include "Error.h"

#include <stdio.h>

static int Usage(const char *argv0, int error) {
  fprintf(
    stdout,
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
  return error;
}

int main(int argc, const char *argv[])
{
  if (argc != 1)
    return Usage(argv[0], IncorrectUsage);
}
