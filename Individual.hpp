#ifndef __INDIVIDUAL__
#define __INDIVIDUAL__

#include <vector>
#include <string>
using namespace std;

/**
 * Se creo una clase para representar cada uno de los elementos
 * en la poblacion
 *
 * Mejia Juarez Nancy Arlette
 * Noviembre 2017
 */
class Individual
{
public:
	Individual();
	// Metodo que permite inicializar el inviduo dado el vector de variables
	void init(vector<double>variables);

	// Parametros de la clase
	vector<double> variables;
	vector<double> objectives;
	vector<Individual*> dom_set;
	int dom_count;
	int rank;
	double distance;
};

#endif