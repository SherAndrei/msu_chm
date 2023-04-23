#pragma once

#include "Point.h"

typedef struct DistanceAndIndexType {
  double distance;
  int index;
} DistanceAndIndex;

// Functions expects:
// 1. all_points is array of size N
// 2. distances is array of size N
// Function guarantees first k distantases is the nearest neighbors
void KNearestNeighbors(const Point *all_points, int N, Point target, int k,
                       DistanceAndIndex *distances);
