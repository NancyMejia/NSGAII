compile:
	g++ -std=c++11 -o output Individual.cpp WFG2.cpp WFG8.cpp main.cpp NSGA2.cpp

clean:
	rm output

execute:
	./output WFG2
	./output WFG8

plot:
	python plot.py