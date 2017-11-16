#include "Individual.hpp"
#include "NSGA2.hpp"
#include "WFG2.hpp"
#include "WFG8.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <climits>
#include <string>
#include <vector>
using namespace std;


/**
 * Tarea 02 Larga - Mejia Juarez Nancy Arlette
 * Implementacion del algoritmo NSGA-II 
 *
 *
 * Archivo con todas las funciones necesarias para resolver las funciones WFG2 Y WFG8 con 
 * el algoritmo NSGA-II.
 */

/**
 * Constructor de la clase NSGA, recibe como parametros:
 * El nombre de la funcion que se quiere resolver
 * El numero de objetivo de  la funcion
 * El tamanio de la poblacion
 * El tamanio de mutaciones
 * El tamanio de cruza
 */
NSGA2::NSGA2(string function, int obj_num, int pop_size, int mutation_num, int crossover_num)
{
	this->function = function;
	this->obj_num = obj_num;
	this->pop_size = pop_size;
	this->mutation_num = mutation_num;
	this->crossover_num = crossover_num;
	if(function == "WFG2")
	{
		this->wfg2 = new WFG2();
		this->wfg2->init(obj_num);
		this->wfg8 = NULL;
	}
	else if(function == "WFG8")
	{
		this->wfg8 = new WFG8();
		this->wfg8->init(obj_num);
		this->wfg2 = NULL;
	}
}

/**
 * Funcion auxiliar que genera la poblacion inicial. 
 * Esta se genera de acuerdo a la funcion que se quiere resolver.
 * Recibe como parametros
 * Nombre de la funcion a resolver 
 * Tamanio de la poblacion que se quiere obtener
 */
vector<Individual*> NSGA2::generatePopulation(string function, int pop_size)
{
	vector<Individual*> population;
	for (int i = 0; i < pop_size; i++)
	{
		Individual *new_individual = new Individual();
		vector<double> variables;
		if(function == "WFG2")
		{
			variables = this->wfg2->WFG_2_thru_7_random_soln(this->wfg2->k, this->wfg2->l);
		}
		else if(function == "WFG8")
		{
			variables = this->wfg8->WFG_8_random_soln(this->wfg8->k, this->wfg8->l);
		}
		new_individual->init(variables);
		population.push_back(new_individual);
	}
	return population;
}

/**
 * Funcion auxiliar que permite evaluar cada uno de los individuos en 
 * la poblacion, se oibtiene como resultado el vector de objetivos
 * @param population -- vector con la poblacion que se desea evaluar
 */
void NSGA2::evaluatePopulation(vector<Individual*> *population)
{
	for (int i = 0; i < this->pop_size; i++)
	{
		vector<double> new_objectives;
		for(int j = 0; j < this->obj_num; j++)
		{
			new_objectives.push_back(0.0);
		}
		if(this->function == "WFG2")
		{
			this->wfg2->evaluate(population->at(i)->variables, new_objectives);
		}
		else if(this->function == "WFG8")
		{
			this->wfg8->evaluate(population->at(i)->variables, new_objectives);
		}
		population->at(i)->objectives = new_objectives;
	}	
}

/**
 * FUncion auxiliar para debuggeo
 * @param population -- poblacion que se desea imprimir
 */
void NSGA2::printPopulation(vector<Individual*> population)
{
	for(int i = 0; i < population.size(); i++)
	{
		cout << "distance " << population.at(i)->distance << endl;
		cout << "rank " << population.at(i)->rank << endl;
		cout << "objectives " ;
		for(int j = 0; j < population.at(i)->objectives.size(); j++)
		{
			cout << population.at(i)->objectives.at(j) << " ";
		}
		cout << endl;
		
		cout << "variables ";
		for(int j = 0; j < population.at(i)->variables.size(); j++)
		{
			cout << population.at(i)->variables.at(j) << " ";
		}
		cout << endl;
	}
}

/**
 * Funcion auxiliar que recibe como parametro dos individyos y se verifica la dominancia
 * @return bool si es que lo domina o no   
 */
bool NSGA2::dominated(Individual* ind1, Individual* ind2)
{
	bool dominated = false;
	if((ind1->objectives.at(0) <= ind2->objectives.at(0) && ind1->objectives.at(1) < ind2->objectives.at(1))  \
		|| (ind1->objectives.at(0) < ind2->objectives.at(0) && ind1->objectives.at(1) <= ind2->objectives.at(1)))
	{
		dominated = true;
	}
	return dominated;
}

/**
 * Funcion auxuliar que obtiene el conjunto de frentes usando el 
 * algoritmo fast non dominated sort.
 * En el primer paso se obtiene el primer frente y despues el resto de los frentes.
 */
vector<vector<Individual*> > NSGA2::fastNonDominatedSort(vector<Individual*> *population)
{

	vector<vector<Individual*> > fronts;
	vector<Individual*> front;
	// first front
	for(int i = 0; i < population->size(); i++)
	{
		population->at(i)->dom_count = 0;
		population->at(i)->dom_set.empty();

		for(int j =0; j< population->size(); j++)
		{
			if(dominated(population->at(i), population->at(j)))
			{
				population->at(i)->dom_set.push_back(population->at(j));
			}
			else if(dominated(population->at(j), population->at(i)))
			{
				population->at(i)->dom_count += 1;
			}
		}
		if(population->at(i)->dom_count == 0)
		{
			population->at(i)->rank = 0;
			front.push_back(population->at(i));
		}
	}
	fronts.push_back(front);
	int curr = 0;
	while(front.size() !=0 )
	{
		vector<Individual*> Q;
		for(int i = 0; i < front.size(); i++)
		{
			for(int j = 0; j < front.at(i)->dom_set.size(); j++)
			{
				front.at(i)->dom_set.at(j)->dom_count -= 1;
				if(front.at(i)->dom_set.at(j)->dom_count == 0)
				{
					front.at(i)->dom_set.at(j)->rank = curr+1;
					Q.push_back(front.at(i)->dom_set.at(j));
				}
			}
		}
		if(Q.size() != 0)
		{
			fronts.push_back(Q);
		}
		curr += 1;
		front = Q;
	}


	// writing fronts on file
	/*
	for(int i = 0; i < fronts.size(); i++)
	{
		char filename_front[180];
		if(this->function == "WFG2")
		{
			strcpy(filename_front,"results/WFG2_fronts_");
		}
		else
		{
			strcpy(filename_front,"results/WFG8_fronts_");
		}
		ofstream fronts_res;
		strcat(filename_front, to_string(i).c_str());
		strcat(filename_front, ".txt");
		cout << filename_front << endl;
		fronts_res.open(filename_front);
		vector<Individual*> front = fronts.at(i);
		for(int j = 0; j < front.size(); j++)
		{
			fronts_res << front.at(j)->objectives.at(0) << " " << front.at(j)->objectives.at(1) << endl;
		}
		fronts_res.close();
	}
	cout << fronts.size() << endl;
	*/
	return fronts;
}

/**
 * Funcion auxiliar que servira para ordenar en orden ascendente
 * el vector de objetivos de acuerdo al valor del primer objetivo
 * @param  ind1 - primer individuo a comparar 
 * @param  ind2 - segundo individuo a comparar
 * @return      - orden  
 */
bool cmp0(const Individual* ind1, const Individual*ind2)
{
	return ind1->objectives.at(0) < ind2->objectives.at(0);
}

/**
 * Funcion auxiliar que servira para ordenar en orden ascendente
 * el vector de objetivos de acuerdo al valor del segundo objetivo
 * @param  ind1 - primer individuo a comparar 
 * @param  ind2 - segundo individuo a comparar
 * @return      - orden  
 */
bool cmp1(const Individual* ind1, const Individual*ind2)
{
	return ind1->objectives.at(1) < ind2->objectives.at(1);
}

/**
 * Funcion auxiliar que servira para ordenar en orden ascendente 
 * La poblacion de acuerdo a su valor de distancai
 * @param  ind1 - primer individuo a comparar
 * @param  ind2 - segundo individuo a comparar
 * @return      - orden
 */
bool cmp3(const Individual* ind1, const Individual*ind2)
{
	return ind1->distance > ind2->distance;
}


//** Using a uniform random distribution, generate a number in [0,bound]. ***
double NSGA2::next_double( const double bound = 1.0 )
{
  assert( bound > 0.0 );

  return bound * rand() / static_cast< double >( RAND_MAX );
}

/**
 * Funcion que se encarga de realiza el calculo de la distancia
 * esto se hace con el fin de onservar la diversidad en el frente
 * @param front - elementos del frente a los cuales se les quiere calcular la
 * distancia
 */
void NSGA2::crowdingDistance(vector<Individual*> *front)
{
	const double MAX_DIST = 1.0e14;
	int N = front->size();
	for(int i = 0; i < N; i++)
	{
		front->at(i)->distance = 0.0;
	}
	for(int m = 0; m < this->obj_num; m++)
	{
		// Se ordena el frente de acuerdo al objetivo
		if(m == 0)
		{
			sort(front->begin(), front->end(), cmp0);
		}
		else if(m == 1)
		{
			sort(front->begin(), front->end(), cmp1);
		}
		// Se establecen las distancias de los extremos a un valor maximo
		front->at(0)->distance = MAX_DIST;
		front->at(N-1)->distance = MAX_DIST;
		double max = front->at(N-1)->objectives.at(m);
		double min = front->at(0)->objectives.at(m);
		double div = max - min;
		for(int i = 1; i < N-1; i++)
		{
			double objective_minus1 = front->at(i-1)->objectives.at(m); 
			double objective_plus1 = front->at(i+1)->objectives.at(m);			
			front->at(i)->distance += (objective_plus1 - objective_minus1)/div;
		}
	}
}

/**
 * Funcion que implementa un torneo binario entre dos individuos de la poblacion
 * @param  ind1 - Primer individuo que se quiere comparar
 * @param  ind2 - Segundo individuo que se quiere comparar 
 * @return      - individuo ganador 
 */
Individual* NSGA2::tournament(Individual*ind1, Individual*ind2)
{
	bool flag = dominated(ind1, ind2);
	if (dominated(ind1, ind2))
	{
		return ind1;
	}
	else if(dominated(ind2, ind1))
	{
		return ind2;
	}
	if(ind1->distance > ind2->distance)
	{
		return ind1;
	}
	if(ind2->distance > ind1->distance)
	{
		return ind2;
	}
	if(next_double() <= 0.5)
	{
		return ind1;
	}
	else
	{
		return ind2;
	}
}

/**
 * Funcion que realiza la cruza entre dos individuos y los almacena en las variables
 * child1, child2
 * @param parent1 - Primer padre de los individuos
 * @param parent2 - Segundo padre de los individuos
 * @param ind1    - Primer hijo generado
 * @param ind2    - Segundo hijo generado
 */
void NSGA2::crossoverInd(Individual *parent1, Individual *parent2, Individual *ind1, Individual*ind2)
{
	double pcross = 0.5;
	double EPS = 1.0e-14;
	double y1, y2, ylow, yupper, rnd, alpha, beta, betaq;
	double c1, c2;
	double eta_c = 10;
	int nvar;
	if(this->function == "WFG2")
	{
		nvar = this->wfg2->k + this->wfg2->l;
	}
	else if(this->function == "WFG8")
	{
		nvar = this->wfg8->k + this->wfg8->l;
	}
	if(next_double() <= pcross)
	{
		for(int i = 0; i < nvar; i++)
		{
			if(next_double() <= 0.5)
			{
				if(fabs(parent1->variables.at(i) - parent2->variables.at(i)) > EPS)
				{	
					if(parent1->variables.at(i) < parent2->variables.at(i))
					{
						y1 = parent1->variables.at(i);
						y2 = parent2->variables.at(i);
					}
					else
					{
						y1 = parent2->variables.at(i);
						y2 = parent1->variables.at(i);
					}

					ylow = 0.0;
					yupper = 2.0 * (i+1);
					rnd = next_double();
					beta = 1.0 + (2.0*(y1- ylow)/(y2-y1));
					alpha = 2.0 - pow(beta, -(eta_c+1.0));
					if(rnd <= (1.0/alpha))
					{
						betaq = pow((rnd*alpha), (1.0/(eta_c+1.0)));
					}
					else
					{
						betaq = pow((1.0/(2.0 - rnd*alpha)), (1.0/(eta_c+1.0))); 
					}
					c1 = 0.5 * ((y1+y2) - betaq*(y2-y1));
					beta = 1.0 + (2.0*(yupper -y2)/(y2-y1));
					alpha = 2.0 - pow(beta, -(eta_c+1.0));
					if(rnd <= (1.0/alpha))
					{
						betaq = pow((rnd*alpha), (1.0/(eta_c+1.0)));
					}
					else
					{
						betaq = pow((1.0/(2.0 - rnd*alpha)), (1.0/(eta_c+1.0))); 
					}
					c2 = 0.5 * ((y1+y2)+betaq*(y2-y1));
					if(c1 < ylow)
						c1 = ylow;
					if (c2 < ylow)
						c2 = ylow;
					if (c1 > yupper)
						c1 = yupper;
					if (c2 > yupper)
						c2 = yupper;
					if(next_double() <= 0.5)
					{
						ind1->variables.at(i) = c2;
						ind2->variables.at(i) = c1;
					}
					else
					{
						ind1->variables.at(i) = c1;
						ind2->variables.at(i) = c2;
					}
				}
				else
				{
					ind1->variables.at(i) = parent1->variables.at(i);
					ind2->variables.at(i) = parent2->variables.at(i);
				}
			}
			else
			{
				ind1->variables.at(i) = parent1->variables.at(i);
				ind2->variables.at(i) = parent2->variables.at(i);
			}
		}
	}
	else
	{
		for(int i = 0; i < nvar; i++)
		{
			ind1->variables.at(i) = parent1->variables.at(i);
			ind2->variables.at(i) = parent2->variables.at(i);
		}
	}
}


/**
 * Funcion que se encarga de hacer la seleccion de los individuos 
 * dada la poblacion de los padres. 
 * Para ello utiliza un torneo y ademas aplica la cruza
 */
vector<Individual*> NSGA2::selection(vector<Individual*> population)
{

	vector<Individual*> selected = generatePopulation(this->function, this->pop_size);
	int a1[this->pop_size];
	int a2[this->pop_size];
	int rnd, rnd2, temp;
	for(int i = 0; i < this->pop_size; i++)
	{
		a1[i] = i;
		a2[i] = i;
	}

	for(int i = 0; i < this->pop_size; i++)
	{
		rnd = rand() % (this->pop_size - i) + i;
		temp = a1[rnd];
		a1[rnd] = a1[i];
		a1[i] = temp;
		rnd2 = rand() % (this->pop_size - i) + i;
		temp = a2[rnd2];
		a2[rnd2] = a2[i];
		a2[i] = temp;
	}
	for(int i = 0; i < this->pop_size; i+=4)
	{
		Individual *parent1 = tournament(population.at(a1[i]), population.at(a1[i+1]));
		Individual *parent2 = tournament(population.at(a1[i+2]), population.at(a1[i+3]));
		crossoverInd(parent1, parent2, selected.at(i), selected.at(i+1));
		Individual *parent3 = tournament(population.at(a2[i]), population.at(a2[i+1]));
		Individual *parent4 = tournament(population.at(a2[i+2]), population.at(a2[i+3]));
		crossoverInd(parent3, parent4, selected.at(i+2), selected.at(i+3));
	}
	return selected;
}

/**
 * Funcion auxiliar que permite concatenar dos poblaciones, a la primera poblacion se
 * le concatenaran los elementos de la segunda poblacion
 * @param parents - poblacion de padres
 * @param front   - frente que se desea concatenar
 */
void NSGA2::appendPopulations(vector<Individual*> *parents,vector<Individual*> front)
{
	for(int i = 0;i < front.size(); i++)
	{
		parents->push_back(front.at(i));
	}
}

/**
 * Funcion auxiliar que permite crear una poblacion nueva que consiste en la union
 * de la poblacion de los padres con la de los hijos.
 * @param population - poblacion de padres 
 * @param crossover - poblacion de hijos
 */
vector<Individual*> NSGA2::mergePopulations(vector<Individual*> population, vector<Individual*> pop_crossover)
{
	for(int i= 0;  i < pop_crossover.size(); i++)
	{
		population.push_back(pop_crossover.at(i));
	}
	return population;
}

/**
 * Funcion que se encarga de mutar un individuo de acuerdo a la probabilidad
 * de mutacion que se establecio.
 * @param ind - individuo que se desea mutar.
 */
void NSGA2::mutateInd(Individual* ind)
{
	double pmut_r = 0.10;
	double eta_m = 50;
	double y, ylow, yupper, delta1, delta2, deltaq, val, xy, rnd, mut_pow;
	int nvar;
	if(this->function == "WFG2")
	{
		nvar = this->wfg2->k + this->wfg2->l;
	}
	else if(this->function == "WFG8")
	{
		nvar = this->wfg8->k + this->wfg8->l;
	}
	for(int j = 0; j < nvar; j++)
	{
		if(next_double() <= pmut_r)
		{
			y = ind->variables.at(j);
			ylow = 0.0;
			yupper = 2.0 * (j+1);
			delta1 = (y-ylow)/(yupper-ylow);
			delta2 = (yupper-y)/(yupper-ylow);
			rnd = next_double();
			mut_pow = 1.0/(eta_m+1.0);
			if(rnd <= 0.5)
			{
				xy = 1.0 - delta1;
				val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
				deltaq = pow(val, mut_pow)- 1.0;
			}
			else
			{
				xy = 1.0 - delta2;
				val = 2.0 *(1.0-rnd)+ 2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
				deltaq = 1.0 - pow(val, mut_pow);
			}
			y = y + deltaq*(yupper-ylow);
			if (y < ylow)
				y = ylow;
			if (y > yupper)
				y = yupper;
			ind->variables.at(j) = y;
		}
	}
}

/**
 * Funcion que se encarga de aplicarle la mutacion a cada uno de los individuos de la poblacion
 * @param population - poblacion que se desea mutar
 */
void NSGA2::mutation(vector<Individual*> *population)
{
	for(int i = 0; i < population->size(); i++)
	{
		mutateInd(population->at(i));
	}
}


/**
 * Funcion principal que se encarga de obtener la solucion al problema multiobjetivo
 * Recibe como parametro el numero de generaciones que se usara como criterio de paro.
 * Devuelve el vector con el frente encontrado
 */
vector<Individual*> NSGA2::solve(int generations)
{
	vector<Individual*> solution;
	vector<Individual*> population = generatePopulation(this->function, this->pop_size);
	evaluatePopulation(&population);
	vector<vector<Individual*> > fronts = fastNonDominatedSort(&population);
	for(int i = 0;i < fronts.size(); i++)
	{
		crowdingDistance(&fronts.at(i));		
	}
	vector<Individual*> selected = selection(population);
	mutation(&selected);
	for(int generation = 0; generation < generations; generation++)
	{
		evaluatePopulation(&selected);
		vector<Individual*> union_individuals = mergePopulations(population, selected);
		fronts = fastNonDominatedSort(&union_individuals);

		vector<Individual*>  parents;
		int index_front = 0;
		vector<Individual*> frontL = fronts.at(index_front);
		int counter = frontL.size();
		while (counter < this->pop_size)
		{
			crowdingDistance(&frontL);	
			appendPopulations(&parents, frontL);	
			index_front+=1;
			frontL = fronts.at(index_front);
			counter +=  frontL.size();
		}	
		int k = this->pop_size - parents.size();
		if(k > 0)
		{
			crowdingDistance(&frontL);
			sort(frontL.begin(), frontL.end(), cmp3);
			for(int i = 0; i < k ; i++)
			{
				parents.push_back(frontL.at(i));
			}
		}
		population = selected;
		selected = selection(parents);
		mutation(&selected);
	}
	evaluatePopulation(&selected);
	for(int i = 0;i < selected.size(); i++)
	{
		solution.push_back(selected.at(i));
	}
	return solution;
}
