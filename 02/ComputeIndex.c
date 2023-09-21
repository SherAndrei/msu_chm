
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// (i_1, i_2, ..., i_{2^{p-1}})
// (2^p+1-i_1, i_1, 2^p+1-i_2, i_2, .., 2^p+1-i_{2^{p-1}}, i_{2^{p-1}})
// 1
// 2 1
// 3 2 4 1
// 6 3 7 2 5 4 8 1

int main(int argc, char** argv) {
  const double eps = 1e-15;
  unsigned N;
  unsigned* indexes;
  unsigned* tmp;
  if (!(argc == 2 && sscanf(argv[1], "%u", &N) == 1 && fabs(log2(N) - floor(log2(N))) < eps)) {
    printf("Usage: %s N\n"
           "\tunsigned N should be a power of two\n", argv[0]);
    return 0;
  }

  tmp = malloc(N * sizeof(*indexes));
  indexes = malloc(N * sizeof(*indexes));

  if (!tmp || !indexes) {
    fprintf(stderr, "not enough memory\n");
    free(tmp);
    free(indexes);
    return 1;
  }

  *tmp = 1u;
  for (unsigned tmp_size = 1; tmp_size < N; tmp_size *= 2) {
    for (unsigned i = 0, j = 0; i < 2 * tmp_size; i += 2, j++) {
      indexes[i + 0] = 2 * tmp_size + 1 - tmp[j];
      indexes[i + 1] = tmp[j];
    }
    memcpy(tmp, indexes, sizeof(*tmp) * 2 * tmp_size);
  }

  for (unsigned k = 1; k <= N; k++) {
    printf("%u%c", indexes[k-1], " \n"[k == N]);
  }

  free(tmp);
  free(indexes);
}
