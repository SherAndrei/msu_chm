#include "step.h"

void Step(int k, double h, double A, double& el, double& next)
{
	if (k <= 1) {
		el = 1.;
		next = 1. - A * h;
		return;
	}
	double prom = 4. * next + (2. * h * A - 3.) * el;
	el = next;
	next = prom;
}
