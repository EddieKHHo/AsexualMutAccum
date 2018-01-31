#include "IncludeFiles.h"


// Declaration of the class "Individual"
class Indiv {

	private:
		
		//Number of del and ben mutations		
		vector<int> numDelMutations;
		vector<int> numBenMutations;
		
		//Total number of mutations in individual
		int numAllDelMutations;
		int numAllBenMutations;
		
		//Age of individual indicates length of time it hasn't gained new mutation (i.e. number of gen within a certain class)
		int age;
		
		//fitness of individual
		double fitness;

	public:
		
		//Constructor to initialize individual with expected num of mutations at mutation-selection balance
		Indiv(const vector<double> &Ud, const vector<double> &sd, const vector<double> &Ub, const vector<double> &sb, int numDelTypes, int numBenTypes);
		//Constructor to initialize individual given num of mutations, age and fitness
		Indiv(const vector<int> &DelMutations, const vector<int> &BenMutations, int Age, const vector<double> &sd, const vector<double> &sb, int numDelTypes, int numBenTypes);
		
		//Calculate fitness of individual based on all mutation types
		void calcFitness(const vector<double> &sd, const vector<double> &sb, int numDelTypes, int numBenTypes);
		
		
	
	
	//Allow World and Popul to access Indiv private variables
	friend class Popul;
	friend class World;

};