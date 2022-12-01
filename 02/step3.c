#include "step.h"

void Step(int k, double h, double A, double* el, double* next)
{
	if (k == 0) {
		*el = 0.;
		*next = 1.;
		return;
	}
	*el = *next;
	*next *= ((1. - h * A / 2.)/(1.0 + h * A / 2.));
}
