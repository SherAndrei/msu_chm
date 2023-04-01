### Task

Find a numerical vector solution of system of `m` equations $F(x)=0,\ x=(x^1,...,x^m)^T$ with an accuracy of `eps` using Newton's method.

#### Build
```bash
$ make
```

#### Solve system of equations (Newton.pdf 2)
See help
```bash
$ ./Newton.out
```
To solve system of 3 equations
$$\begin{cases}
  x^2+y^2+z^2=10 \\
  -x^2+(y-3)^2=z \\
  x+2-\frac{y}{3}+z=0
\end{cases}$$
witn precision of $1e-10$ and with initital approximation $(0, 0.5, 1)$ use:
```bash
$ ./Newton 1e-10 3
   -3.01007258012015
   -0.16475754761303
    0.95515339758247
```
To solve next system of equations with initial approximation $(0.5,0.5)$ instead
$$\begin{cases}
  x^{2}-y^{2}=1 \\
  x^{2}+y^{2}=15
\end{cases}$$
modify `SystemOfEquations.c` in the following way:
```c
void InitialApproximation(double *x, unsigned m) {
  (void)m;
  x[0] = 0.5;
  x[1] = 0.5;
}

void F(double *f, const double *x, unsigned m) {
  (void)m;
  f[0] = pow(x[0], 2) - pow(x[1], 2) - 1;
  f[1] = pow(x[0], 2) + pow(x[1], 2) - 15;
}

void dF(double *jacobian, const double *x, unsigned m) {
  (void)m;
  J(0, 0) = 2. * x[0];
  J(0, 1) = -2. * x[1];
  J(1, 0) = 2. * x[0];
  J(1, 1) = 2. * x[1];
}
```
Then rebuild and run:
```bash
$ make
$ ./Newton.out 1e-10 2
    2.82842712474619
    2.64575131106459
```
