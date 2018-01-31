#include "Indiv.h"

//Constructor for first generation of individuals
//Initializes number of mutations as expected under mutation selection balance
Indiv::Indiv(const vector<double> &Ud, const vector<double> &sd, const vector<double> &Ub, const vector<double> &sb, int numDelTypes, int numBenTypes){

	// Initialize random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed += time(NULL));
	
	//Create number of mutations at mutation-selection balance (Poisson distributed with mean U/s)
	vector<double> expMeanDelMut(numDelTypes);
	vector<double> expMeanBenMut(numBenTypes);
	
	numAllDelMutations = 0;
	numAllBenMutations = 0;
	
	
	//loop through all mutation types and add expected number of each type to numDelMutations/numBenMutations vector
	//also calc total number of del and ben mutations
	for(int i = 0; i < numDelTypes; i ++){
	
		if(sd[i] == 0){
			expMeanDelMut[i] = 0;
		}
		else{
			expMeanDelMut[i] = Ud[i]/sd[i];
		}
		
		//Set Indiv num of del mut
		int nTemp = gsl_ran_poisson(r, expMeanDelMut[i]);
		
		//numDelMutations.push_back(nTemp);
		
		//For MA paper (with Pepijn) I start simulation with 0 deleterious mutations
		numDelMutations.push_back(0);
		numAllDelMutations += nTemp;
	}
	for(int i = 0; i < numBenTypes; i ++){
		if(sb[i] == 0){
			expMeanDelMut[i] = 0;
		}
		else{
			expMeanBenMut[i] = Ub[i]/sb[i];
		}
		
		//Set Indiv num of ben mut
		int nTemp = gsl_ran_poisson(r, expMeanBenMut[i]);
		numBenMutations.push_back(nTemp);
		numAllBenMutations += nTemp;
	}
	
		
	//adjust seed by adding time
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//Set Indiv age
	age = 1;

	//calculate fitness
	calcFitness(sd, sb, numDelTypes, numBenTypes);
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);

}



//Constructor for individuals after first generation
//Individual is given number of mutations, age and fitness
Indiv::Indiv(const vector<int> &DelMutations, const vector<int> &BenMutations, int Age, const vector<double> &sd, const vector<double> &sb, int numDelTypes, int numBenTypes){
	
	//initialize total mutations as 0
	numAllDelMutations = 0;
	numAllBenMutations = 0;
	
	//loop through all mutation types and add to vector numDelMutations
	//also calc total number of del and ben mut
	for(int i = 0; i < numDelTypes; i++){
		numDelMutations.push_back(DelMutations[i]);
		numAllDelMutations += numDelMutations[i];
	}
	for(int i = 0; i < numBenTypes; i++){
		numBenMutations.push_back(BenMutations[i]);
		numAllBenMutations += numBenMutations[i];
	}
	
	//Set Indiv age
	age = Age;
	
	//calculate fitness
	calcFitness(sd, sb, numDelTypes, numBenTypes);
}



//Calculate fitness of individual given vector of selection coefficients for different mutation types
void Indiv::calcFitness(const vector<double> &sd, const vector<double> &sb, int numDelTypes, int numBenTypes){
	
	//Initialize Indiv fitness as 1
	fitness = 1;
	
	//given selection coefficients, calc multiplicative fitness
	for(int i = 0; i < numDelTypes; i++){
		fitness *= pow( 1 - sd[i], numDelMutations[i]);
	}
	//loop through all ben mut types
	for(int i = 0; i < numBenTypes; i++){
		fitness *= pow( 1 + sb[i], numBenMutations[i]);
	}
}




