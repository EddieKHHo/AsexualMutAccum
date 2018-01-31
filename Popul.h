#include "Indiv.h"


class Popul {

	private:
		
		//variables that are initialized as parameters, must be updated after every event (mutation, selection)
		
		//list of individuals in population
		vector<Indiv> listIndividuals;
		//population size
		int sizePop;
		//number of types of mutations for del and ben mut
		int numDelMutTypes;
		int numBenMutTypes;
		
		
		//variables that are calculated based on current population, does NOT need to be updated after every event (mutation, selection)
		//If not updated by calling their respective functions, then data will be based on the previous update (i.e. be careful about when data functions are called)		
		
		//highest and mean fitness of population
		double highestFitness;
		double meanFitness;
		
		//mean number of mutations either as total or each type
		double meanDel_All;
		double meanBen_All;
		vector<double> meanDel_Types;
		vector<double> meanBen_Types;
		
		//LLC based on total mutations of by each type and the number of individuals within each LLC
		int DelLLC_All;
		int numDelLLC_All;
		vector<int> DelLLC_Types;
		vector<int> numDelLLC_Types;

	public:
		
		//Default constructor
		Popul();
		
		//Constructor for first generation of population
		Popul(const vector<double> &Ud, const vector<double> &sd, const vector<double> &Ub, const vector<double> &sb, int numDelTypes, int numBenTypes, int N);
		
		//second constructor, not used
		Popul(vector<Indiv> &newIndividuals);
		
		//add mutation to population; requires U for generating mutations; requires s since new Indiv are created and need to calculate fitness
		void doMutation(const vector<double> &Ud, const vector<double> &sd, const vector<double> &Ub, const vector<double> &sb);
		
		//select N number of offspring for next generation; ASSUMES fitness is updated after mutation
		void doSelection(int N);
		
		//create next generation of population by updating listIndividuals (class private variable)
		void makeNextGen(const vector<Indiv> &listIndividualsNextGen, int N);
		
		//find fitness of individual with highest fitness
		void findHighestFitness();
		//find mean fitness of population
		void findMeanFitness();
		//calc mean # del mutation of each type and for whole individual
		void calcMeanDel();
		//calc mean # ben mutation of each type and for whole individual		
		void calcMeanBen();
		//find # of mutations in LLC for each mut type and for total sum of mutations
		//also find number of individuals that fall into those classes
		void findDelLLC();
	
	
	//Allow World to access Popul private variables
	friend class World;
};