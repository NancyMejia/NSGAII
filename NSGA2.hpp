#ifndef __NSGA2_H__
#define __NSGA2_H__

#include "WFG2.hpp"
#include "WFG8.hpp"
#include <vector>
#include <string>
using namespace std;

class NSGA2
{
public:
	// Constructor de la clase
	NSGA2(string function, int obj_num, int pop_size, int mutation_num, int crossover_num);
	// Funcion principal
	vector<Individual*> solve(int generations);
	// parametros de la clase
	string function;
	int obj_num;
	int pop_size;
	int mutation_num;
	int crossover_num;
	WFG2 *wfg2;
	WFG8 *wfg8;
private:
	// Funciones que se usan dentro de NSGA2
	vector<Individual*> generatePopulation(string function, int pop_size);
	void evaluatePopulation(vector<Individual*> *population);
	vector<vector<Individual*> > fastNonDominatedSort(vector<Individual*> *population);
	void crowdingDistance(vector<Individual*> *front);
	vector<Individual*> selection(vector<Individual*> population);
	void mutation(vector<Individual*> *population);
	
	// Funciones auxiliares NSGA2
	Individual* tournament(Individual*ind1, Individual*ind2);
	void crossoverInd(Individual *parent1, Individual *parent2, Individual *ind1, Individual*ind2);
	void mutateInd(Individual* ind);

	// Funciones auxiliares
	void printPopulation(vector<Individual*> population);
	bool dominated(Individual* ind1, Individual* ind2);
	double next_double( const double bound );
	void appendPopulations(vector<Individual*> *parents,vector<Individual*> front);
	vector<Individual*> mergePopulations(vector<Individual*> population, vector<Individual*> pop_crossover);

};

#endif