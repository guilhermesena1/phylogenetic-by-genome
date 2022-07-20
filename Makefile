all : phylo

phylo: src/phylo.cpp
	g++ -O3 -Wall -std=c++11 -o phylo src/phylo.cpp -fopenmp

clean:
	rm phylo
