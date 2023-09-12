## Solve numerically $y'=-Ay$, $y(0)=1$

### Build
```
make
```

### Run

Seek help using program named `./a.out`
```bash
$ ./a.out
Usage: ./a.out N X A
        unsigned N: amount of segment divisions
        double X: length of the segment
        double A: task parameter
```

$$\frac{y_{k+1} - y_{k}}{h} + Ay_{k} = 0, y_0 = 1$$
```bash
$ ./ExplicitEuler.out 10 1 1 > out.txt
```

$$\frac{y_{k+1} - y_{k}}{h} + Ay_{k+1} = 0, y_0 = 1$$
```bash
$ ./ImplicitEuler.out 10 1 1 > out.txt
```

$$\frac{y_{k+1} - y_{k}}{h} + A\frac{y_{k+1} + y_{k}}{2} = 0, y_0 = 1$$
```bash
$ ./ImplicitAdams.out 10 1 1 > out.txt
```

$$\frac{y_{k+1} - y_{k-1}}{2h} + Ay_k = 0, y_0 = 1, y_1 = 1 - Ah$$
```bash
$ ./Alternating.out 10 1 1 > out.txt
```

$$\frac{1.5y_{k} - 2y_{k-1} + 0.5y_{k-2}}{h}+Ay_k= 0, y_0 = 1, y_1 = 1 - Ah.$$
```bash
$ ./AlphaStable.out 10 1 1 > out.txt
```

$$\frac{-0.5y_{k+2}+2y_{k+1}-1.5y_k}{h} + Ay_k = 0, y_0 = 1, y_1 = 1 - Ah.$$
```bash
$ ./NotAlphaStable.out 10 1 1 > out.txt
```

### Plot
```bash
$ gnuplot -c compare.gnuplot out.txt out.png
```
