#Usage: 
#make:       #compile and link binary
#make clean  #clean the object and executable files

.PHONY: all FFT.o utilities.o main.o clean

CC = g++ #compiler for c++ projects

STANDARD = -std=c++20    #specify the compiling standard

all: fft_transform

fft_transform: main.o utilities.o FFT.o
	@echo "Linking FFT object files..."
	$(CC) $(STANDARD)  -o fft_transform main.o FFT.o utilities.o

FFT.o:  FFT.cpp fft.hpp utilities.hpp
	@echo "Compiling source_1 (FFT.cpp) files"
	$(CC) $(STANDARD) -c FFT.cpp -Wall

utilities.o: utilities.cpp utilities.hpp
	@echo "Compiling source_2 (utilities.cpp) files"
	$(CC) $(STANDARD) -c utilities.cpp -Wall

main.o: main.cpp
	@echo "Compiling main.cpp files..."
	$(CC) $(STANDARD) -c main.cpp  -Wall


clean :
	@echo "Cleaning up..."
	rm main.o FFT.o utilities.o fft_transform
