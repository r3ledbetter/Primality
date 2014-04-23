GeneratePrime: GeneratePrime.cpp
	g++ -fopenmp GeneratePrime.cpp -o GeneratePrime -lmpfr -lgmp

clean:
	rm -f GeneratePrime
