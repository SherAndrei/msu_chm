#include "KNearestNeighbors.h"

#include <stdlib.h>

static int ByDistance(const void *p_lhs, const void *p_rhs) {
  const double lhs = ((const DistanceAndIndex *)p_lhs)->distance;
  const double rhs = ((const DistanceAndIndex *)p_rhs)->distance;
  return (lhs > rhs) - (lhs < rhs);
}

void KNearestNeighbors(const Point *all_points, int N, Point target, int k,
                       DistanceAndIndex *distances, Point *neighbors) {
  for (int i = 0; i < N; ++i) {
    distances[i].distance = Distance(target, all_points[i]);
    distances[i].index = i;
  }

  // FIXME: partial sort should be here but im lazy to implement it
  qsort(distances, N, sizeof(*distances), ByDistance);

  for (int i = 0; i < k; ++i) {
    neighbors[i] = all_points[distances[i].index];
  }
}
