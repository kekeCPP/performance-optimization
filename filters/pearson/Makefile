# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -Wunused -Wall -Wunused

all: pearson pearson_par

pearson: vector dataset analysis pearson.cpp
	$(CXX) $(CXXFLAGS) pearson.cpp vector.o dataset.o analysis.o -o pearson

pearson_par: vector dataset analysis pearson_par.cpp
	$(CXX) $(CXXFLAGS) pearson_par.cpp vector.o dataset.o analysis.o -o pearson_par

analysis: vector analysis.hpp analysis.cpp
	$(CXX) $(CXXFLAGS) -c analysis.cpp -o analysis.o

dataset: vector dataset.hpp dataset.cpp
	$(CXX) $(CXXFLAGS) -c dataset.cpp -o dataset.o

vector: vector.hpp vector.cpp
	$(CXX) $(CXXFLAGS) -c vector.cpp -o vector.o

clean:
	rm -rf pearson *.o *.dSYM 2> /dev/null
