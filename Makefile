all: compile $(wildcard *.py)
	./integrator
	gnuplot plot.plt

compile: $(wildcard *.cpp)
	g++ -O2 *.cpp -o integrator

clean:
	rm -rf result/*
