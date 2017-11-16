#include "Individual.hpp"
#include "WFG2.hpp"
#include "WFG8.hpp"
#include "NSGA2.hpp"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
using namespace std;

/**
 * Mejia Juarez Nancy Arlette
 * Noviembre 2017
 */


/**
 * FUncion principal con la cual se probara el algoritmo NSGA2
 * Recibe como pametro el nombre de la funcion que se desea resolver
 * 
 */
int main(int argc, char  **argv)
{
	srand(time(NULL));
	if(argc < 2)
	{
		cout << "ERROR: El algoritmo debe recibir como parametro el nombre de la funcion a resolver WFG2/WFG8" << endl;
		return -1;
	}
	string function = argv[1];
	cout << "funcion: "<< function << endl;
	int pop_size = 100;
	int obj_num = 2;
	int mutation_num = 50;
	int crossover_num = 10;
	int generations = 350;
	int trials = 5;
	NSGA2 *nsga2 =  new NSGA2(function, obj_num, pop_size, mutation_num, crossover_num);
	for(int trial = 0; trial < trials; trial++)
	{
		char filename_solution[180];
		vector<Individual*> solution = nsga2->solve(generations);
		//Archivo en el que se escribiran las soluciones
		ofstream outputfile;
		if (function == "WFG2")
		{
			strcpy(filename_solution, "results/WFG2_solution");
		}
		else if(function == "WFG8")
		{
			strcpy(filename_solution, "results/WFG8_solution");
		}

		strcat(filename_solution, to_string(trial).c_str());
		strcat(filename_solution, ".txt");
		cout << filename_solution << endl;
		
		outputfile.open(filename_solution);
		
		for(int i = 0; i < solution.size();i++)
		{
			for(int j = 0; j < obj_num; j++)
			{
				outputfile << solution.at(i)->objectives.at(j)  << " ";
			}
			outputfile << endl;
		}		
		// Cerramos el archivo
		outputfile.close();
		
	}	
	delete nsga2;
	return 0;
}