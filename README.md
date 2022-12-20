# amsc_project
## Project for AMSC 22/23: Implementation of FFT algorithm

### After cloning the directory you can use the 'make' command to build the executable, it will be found in the build directory.
## To run the file ./fft use the flags:

```
-r Recursive FFT
```

```
-t Iterative FFT (Cooley-Tukey Algorithm)
```

```
-p Parallel FFT with MPI
```

```
-i Inverse FFT
```

```
-s <file.txt> will read the input from a txt file you created. 
   It should be located in the build the directory and the format should be:

(1.,0.)
(2.1,1.5)
(6.4,10.2)

Otherwise the input is fixed.
```
