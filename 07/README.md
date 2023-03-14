### Task

Calculate the approximate value of the 1-dimensional definite integral using the method of quadrature formulas of Gauss and Simspon

#### Build
```bash
$ make
```

#### Integral
See help
```bash
$ ./Gauss.out
```
```bash
$ ./Simpson.out
```
To find out result of the $\int_0^\pi\cos 100dx$ dividing the segment in 1000 parts:
1. In `main.c` file change main function in the following way
  ```diff
  - result = Integral(a, b, sin, N);
  + result = Integral(a, b, Cos100, N);
  ```
2. Run program
  ```bash
  $ ./Simpson.out 0 3.14159265359 1000
          2.399769e-13
  ```
  ```bash
  $ ./Gauss.out 0 3.14159265359 1000
          2.400901e-13
  ```
