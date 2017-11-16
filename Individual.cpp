#include "Individual.hpp"
#include <iostream>
#include <vector>
using namespace std;

/**
 * Implementacion de las funciones de la clase Individuo
 *
 * Mejia Juarez Nancy Arlette
 * Noviembre 2017
 */

/**
 * FUncion que inicializa el vector de variables de un individuo
 * @param variables vector de variables
 */
void Individual::init(vector<double>variables)
{
	this->variables = variables;
	this->dom_count = 0;
	this->rank = -1;
	this->distance = 0.0;
}

/**
 * Se sobre escribio el constructor de la clase
 */
Individual::Individual()
{
	this->dom_count = 0;
	this->rank = -1;
	this->distance = 0.0;
}
