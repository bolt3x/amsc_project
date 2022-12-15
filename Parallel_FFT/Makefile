#Usage: 
#make:       #compile and link binary
#make clean  #clean the object and executable files

.PHONY: all objects fft_transform execute clean

CC = mpic++ #compiler for c++ projects
EX = mpiexec

STANDARD = -std=c++20    #specify the compiling standard


all: execute

execute: fft_transform
	@echo "Executing program with #8 processes..."
	$(EX) -n 2 fft_transform 8

fft_transform: main.o objects
	@echo "Linking FFT object files..."
	$(CC) $(STANDARD) -o fft_transform main.o FFT.o utilities.o

objects:  FFT.cpp utilities.cpp fft.hpp utilities.hpp
	@echo "Compiling source (FFT.cpp and utilities.cpp) files"
	$(CC) $(STANDARD) -c FFT.cpp utilities.cpp -Wall

main.o: main.cpp fft.hpp utilities.hpp
	@echo "Compiling main.cpp files..."
	$(CC) $(STANDARD) -c main.cpp  -Wall


clean :
	@echo "Cleaning up..."
	rm main.o FFT.o utilities.o fft_transform
