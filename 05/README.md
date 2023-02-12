#### Build
```bash
$ make
```

#### Generate input file
See help
```bash
$ ./Generate.out -h
```
To generate points with 100 points between 2. and 5. using equally distributed nodes and save to file "input.txt" use
```bash
$ ./Generate.out 100 2 5 -d > input.txt
```
To generate points with 20 points between -1 and 1 using Chebyshev nodes and save to file "data.txt" use
```bash
$ ./Generate.out 20 -1 1 -c > data.txt
```

#### Interpolation polynom
See help
```bash
$ ./InterpolationPolynom.out -h
```
By default `InterpolationPolynom.out` reads input from `stdin`. To use generated into file "intput.txt" input use
```bash
$ ./InterpolationPolynom.out < input.txt
```
To save interpolation result to a file "output.txt" use
```bash
$ ./InterpolationPolynom.out < input.txt > output.txt
```
To plot result from file "output.txt" with `gnuplot` into "result.png" use
```bash
$ gnuplot -c cmp.gnuplot output.txt result.png
```
