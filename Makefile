#Usage: 
#make:       #compile and link binary
#make clean  #clean the object and executable files

.PHONY: all objects main.o clean

CC = g++ #compiler for c++ projects

STANDARD = -std=c++20    #specify the compiling standard

all: fft_transform

fft_transform: main.o objects
	@echo "Linking FFT object files..."
	$(CC) $(STANDARD)  -o fft_transform main.o FFT.o utilities.o

objects:  FFT.cpp utilities.cpp fft.hpp utilities.hpp
	@echo "Compiling source (FFT.cpp and utilities.cpp) files"
	$(CC) $(STANDARD) -c FFT.cpp utilities.cpp -Wall

main.o: main.cpp fft.hpp utilities.hpp
	@echo "Compiling main.cpp files..."
	$(CC) $(STANDARD) -c main.cpp  -Wall


clean :
	@echo "Cleaning up..."
	rm main.o FFT.o utilities.o fft_transform
