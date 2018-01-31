#include "Popul.h"

//Default constructor puts all private variables as 1
Popul::Popul(){

	numDelMutTypes = 1;
	numBenMutTypes = 1;
	sizePop = 3;
	
	vector<double> defaultUd(numDelMutTypes,1);
	vector<double> defaultsd(numDelMutTypes,1);
	vector<double> defaultUb(numDelMutTypes,1);
	vector<double> defaultsb(numDelMutTypes,1);
	
	//create sizePop=1 individual by default
	for(int i = 0; i < sizePop; i++){
		Indiv newIndividual(defaultUd, defaultsd, defaultUb, defaultsb, numDelMutTypes, numBenMutTypes);
		listIndividuals.push_back(newIndividual);
	}
	
}

//Constructor for first generation
//Creates list of new individuals given mutation and selection values
Popul::Popul(const vector<double> &Ud, const vector<double> &sd, const vector<double> &Ub, const vector<double> &sb, int numDelTypes, int numBenTypes, int N){
	
	//intialize class private variables
	numDelMutTypes = numDelTypes;
	numBenMutTypes = numBenTypes;
	sizePop = N;
	
	//create sizePop number of new individuals and add to class private variable listIndividuals
	for(int i = 0; i < sizePop; i++){
		Indiv newIndividual(Ud, sd, Ub, sb, numDelTypes, numBenTypes);
		listIndividuals.push_back(newIndividual);
	}
	
	//calculate data for population
	findDelLLC();
	calcMeanDel();
	calcMeanBen();
	findHighestFitness();
	findMeanFitness();
}


//Constructor that creates list of individuals given vector of Indiv
//May not be used?
Popul::Popul(vector<Indiv> &newIndividuals){
	
	listIndividuals = newIndividuals;
	
	findDelLLC();
	calcMeanDel();
	calcMeanBen();
	findHighestFitness();
	findMeanFitness();

}


//Add new mutations to existing individuals in population
//doMutation require Ud and Ub to add mutations
//doMutation also requires sd and sb because it changes number of mutations in individuals, so need to recalculate fitness as well, which requires selection coefficients
//currently doMutation makes new individuals after adding mutations
//	it requires sd and sb to recalcualte individual fitness and highest fitness in population (selection uses individual fitness / highest fitness)
//	sd and sb can change over time and this is implemented HERE; fitness MUST be recalculates based on updated sd and sb EVEN if no mutations are addded
void Popul::doMutation(const vector<double> &Ud, const vector<double> &sd, const vector<double> &Ub, const vector<double> &sb){

	// Initialize random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed += time(NULL));
	
	//vector to hold list of Indiv after mutation
	vector<Indiv> listIndividualsAfterMutation;
	
	//loop through all sizePop individuals
	for(int i = 0; i < sizePop; i ++){
		//vectors to hold list of new mutations for each type
		//reset for each individual
		vector<int> tempNewDelMut;
		vector<int> tempNewBenMut;
		
		//loop through all mut types and add mutation based on their mut rate
		for(int j = 0; j < numDelMutTypes; j ++){
			//generate new mutation		
			int newMut = gsl_ran_poisson(r, Ud[j]);
			//Add new mutation to existing number of ind i
			tempNewDelMut.push_back( listIndividuals[i].numDelMutations[j] + newMut );
		}
		//do same for del mut as ben mut
		for(int j = 0; j < numBenMutTypes; j ++){
			int newMut = gsl_ran_poisson(r, Ub[j]);
			tempNewBenMut.push_back( listIndividuals[i].numBenMutations[j] + newMut );
		}
		
		//create newIndividual, which is the same as listIndividuals[i] except number of mutations has changed based on tempNewDelMut, tempNewBenMut OR selection has changed
		//utilizes the second constructor in class Indiv
		Indiv newIndividual(tempNewDelMut, tempNewBenMut, listIndividuals[i].age, sd, sb, numDelMutTypes, numBenMutTypes);
	
		//adjust seed by adding time
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
		
		//Add newIndividual to vector of Indiv after mutation
		listIndividualsAfterMutation.push_back(newIndividual);
	}
	
	//Update listIndividuals to new list AFTER mutations added
	listIndividuals = listIndividualsAfterMutation;
	
	//Need to recalculate highest fitness based on sd and sb since new individuals may have changed in fitness due to new mutations or altered sd, sb
	findHighestFitness();
	
	/*
	//Update data after mutation (commented out since not needed)
	//We do NOT need to re-calc data just for mutation, we should only do so after selection
	findDelLLC();
	calcMeanDel();
	calcMeanBen();
	findMeanFitness();
	*/
	
	//adjust seed by adding time
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);

}

//select offspring for next generation base on their fitness relative to the highest fitness in population
//input parameter N allows population size of next gen to change; will only select N offspring for next gen
void Popul::doSelection(int N){
	

	// Initialize random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed += time(NULL));
	
	//vector to hold list of Indiv after selection
	vector<Indiv> listIndividualsAfterSelection;
	//count of offspring that pass selection
	int numOffspring = 0;
	
	//while loop continues until numOffspring < N, stops when numOffspring == N
	//We will choose N individuals regardless of what sizePop is
	while(numOffspring < N){
		//choose ind k (ind 0 to sizePop-1)
		int k = gsl_rng_uniform_int(r, sizePop);
		
		//calc fitness of individual k relative to highest fitness in population
		double kRelativeFitness = listIndividuals[k].fitness / highestFitness;
		
		//random number from 0 to 1 as threshold to pass selection
		double requiredFitness = gsl_ran_flat(r, 0, 1);
		
		//Only increase numOffspring if pass selection
		if( kRelativeFitness > requiredFitness ){
			
			//Add newIndividual to vector of Indiv after selection
			listIndividualsAfterSelection.push_back(listIndividuals[k]);
			//increase offspring count by 1
			numOffspring++;
		}
	}
	
	//Create next generation of population
	//Need to update listIndividuals and N since we only selected N individuals, class private popSize would need to change
	makeNextGen(listIndividualsAfterSelection, N);
	
	//adjust seed by adding time
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
	
}

//Change listIndividuals to given listIndividualsNextGen and re-calc data
//Usually done AFTER selection
void Popul::makeNextGen(const vector<Indiv> &listIndividualsNextGen, int N){
	
	//clear existing listIndividuals
	listIndividuals.clear();
	
	//update listIndividuals with listIndividualsNextGen
	listIndividuals = listIndividualsNextGen;
	
	//Update sizePop as given N since we only chose N individuals for next gen during selection
	sizePop = N;
		
	//Update data after mutation
	findDelLLC();
	calcMeanDel();
	calcMeanBen();
	findHighestFitness();
	findMeanFitness();
	

}



//find highest fitness in population
void Popul::findHighestFitness(){

	highestFitness = listIndividuals[0].fitness;

	for(int i = 1; i < sizePop; i++){
		double tempFitness = listIndividuals[i].fitness;
		
		if( tempFitness > highestFitness){
			highestFitness = tempFitness;
		}
	}
}

//find mean fitness of population
void Popul::findMeanFitness(){

	double totalFitness = listIndividuals[0].fitness;

	for(int i = 1; i < sizePop; i++){
		totalFitness += listIndividuals[i].fitness;
	}
	
	meanFitness = totalFitness/sizePop;
}



//Finds LLC for deleterious mutations
void Popul::findDelLLC(){
	
	DelLLC_Types.clear();
	numDelLLC_Types.clear();
	
	vector<int> tempDelLLC_Types;
	
	//Initialize LLC as first individual
	DelLLC_All = listIndividuals[0].numAllDelMutations;
	numDelLLC_All = 1;
	
	//initialize DelLLC_Types, numDelLLC_Types and tempDelLLC_Types using ind 0
	for(int i = 0; i < numDelMutTypes; i++){			
		DelLLC_Types.push_back( listIndividuals[0].numDelMutations[i] );
		numDelLLC_Types.push_back( 1 );
		
		tempDelLLC_Types.push_back( listIndividuals[0].numDelMutations[i] );
	}
		
	//i starts as individual 1, since using ind 0 as initial LLC
	for(int i = 1; i < sizePop; i++){

		//Find DelLLC_All
		int tempDelLLC_All;
		tempDelLLC_All = listIndividuals[i].numAllDelMutations;
		
		//Define new LLC if ind i has less mut than current LLC; reset count of LLC
		if(tempDelLLC_All < DelLLC_All){
			DelLLC_All = tempDelLLC_All;
			numDelLLC_All = 1;
		}
		//increase count if another ind has same # of mut as LLC
		else if (tempDelLLC_All == DelLLC_All){
			numDelLLC_All += 1;
		}
		else{
			//Nothing happens when tempDelLLC_All > DelLLC_All
		}
		
		
		//Find DelLLC for each mutation type
		//Need to go through each mutation type (j) per individual (i)
		for(int j = 0; j < numDelMutTypes; j++){
		
			tempDelLLC_Types[j] = listIndividuals[i].numDelMutations[j];
			
			//Define new LLC if ind i has less mut than current LLC; reset count of LLC
			if( tempDelLLC_Types[j] < DelLLC_Types[j]){
				DelLLC_Types[j] = tempDelLLC_Types[j];
				numDelLLC_Types[j] = 1;
			}
			//increase count if another ind has same # of mut as LLC
			else if (tempDelLLC_Types[j] == DelLLC_Types[j]){
				numDelLLC_Types[j] += 1;
			}
			else{
				//Nothing happens when tempDelLLC_Types[j] > DelLLC_Types[j]
			}
		}
	}	
}





//Calculate mean number of total del mutations in population
void Popul::calcMeanDel(){

	meanDel_Types.clear();
	
	//Initialize total del using first individual (element 0)
	double totalMut_All = listIndividuals[0].numAllDelMutations;
	vector<double> totalMut_Types( listIndividuals[0].numDelMutations.begin(), listIndividuals[0].numDelMutations.end() );
	
	//loop through ind all individuals to get sum of mutations (start at ind 1 since initialized with ind 0)
	for(int i = 1; i < sizePop; i++){
		totalMut_All += listIndividuals[i].numAllDelMutations;
		
		for(int j = 0; j < numDelMutTypes; j++){
			totalMut_Types[j] += listIndividuals[i].numDelMutations[j];
		} 
	}
		
	//calc mean for total mut
	meanDel_All = totalMut_All/sizePop;
	//calc mean for each type of mut
	for(int i = 0; i < numDelMutTypes; i++){
		meanDel_Types.push_back( totalMut_Types[i]/sizePop );
	}
	
}


//Calculate mean number of total del mutations in population
void Popul::calcMeanBen(){

	meanBen_Types.clear();
	
	//Initialize total ben using first individual (element 0)
	double totalMut_All = listIndividuals[0].numAllBenMutations;
	vector<double> totalMut_Types( listIndividuals[0].numBenMutations.begin(), listIndividuals[0].numBenMutations.end() );

	//loop through ind all individuals to get sum of mutations (start at ind 1 since initialized with ind 0)
	for(int i = 0; i < sizePop; i++){
		totalMut_All += listIndividuals[i].numAllBenMutations;
		
		for(int j = 0; j < numBenMutTypes; j++){
			totalMut_Types[j] += listIndividuals[i].numBenMutations[j];
		} 
	}
	
	//calc mean for total mut	
	meanBen_All = totalMut_All/sizePop;
	//calc mean for each type of mut	
	for(int i = 0; i < numBenMutTypes; i++){
		meanBen_Types.push_back( totalMut_Types[i]/sizePop );
	}
	

}







