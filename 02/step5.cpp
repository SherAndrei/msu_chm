#include "step.h"

unsigned m() { return 2; }

void Step(int k, double h, double A, double& el, double& next)
{
    if (k <= 1) {
        el = 1.0;
        next = 1.0 - A * h;
        return;
	}
    double prom = 2 * next - 0.5 * el/(1.5 + h * A);
    el = next;
    next = prom;
}
