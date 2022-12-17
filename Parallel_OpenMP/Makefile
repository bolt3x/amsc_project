#Usage: 
#make:       #compile and link binary
#make clean  #clean the object and executable files

.PHONY: all objects fft_transform execute clean

CC = g++ #compiler for c++ projects

STANDARD = -std=c++20     #specify the compiling standard
OPEN = -fopenmp

all: execute

execute: fft_transform
	@echo "Executing program with #8 processes..."
	./fft_transform 8 4

fft_transform: main_OpenMP.o objects
	@echo "Linking FFT object files..."
	$(CC) $(STANDARD) $(OPEN) -o fft_transform main_OpenMP.o FFT_OpenMP.o utilities_OpenMP.o

objects:  FFT_OpenMP.cpp utilities_OpenMP.cpp fft_OpenMP.hpp utilities_OpenMP.hpp
	@echo "Compiling source (FFT_OpenMP.cpp and utilities_OpenMP.cpp) files"
	$(CC) $(STANDARD) $(OPEN) -c FFT_OpenMP.cpp utilities_OpenMP.cpp -Wall

main.o: main_OpenMP.cpp fft_OpenMP.hpp utilities_OpenMP.hpp
	@echo "Compiling main_OpenMP.cpp files..."
	$(CC) $(STANDARD) $(OPEN) -c main_OpenMP.cpp  -Wall


clean :
	@echo "Cleaning up..."
	rm main_OpenMP.o FFT_OpenMP.o utilities_OpenMP.o fft_transform
