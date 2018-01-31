#include "World.h"

//Default Constructor sets default parameters
World::World(){
	
	NUM_GEN = 1;
	N = 2;
	UD.push_back(1);
	SD.push_back(1);
	UB.push_back(1);
	SB.push_back(1);
	NUM_DEL_TYPES = 1;
	NUM_BEN_TYPES = 1;
	
	//default no fluctuation
	N_F = 1;		
	N_DELTA_A = 1; 
	N_DELTA_B = 1;
	N_PROB_A = 1;
	
	//default no fluctuation
	SD_F.push_back(1);		
	SD_DELTA_A.push_back(1);
	SD_DELTA_B.push_back(1);
	SD_PROB_A.push_back(1);
	
	D_N_SDEL0 = 0;
	N_SDEL0_F = 1;
}

//set basic parameters for simulation (minimum required to run it)
void World::setParBasic(int ng, int n, const vector<double> &ud, const vector<double> &sd, const vector<double> &ub, const vector<double> &sb, int ndt, int nbt){
	
	//get basic parameters
	NUM_GEN = ng;
	N = n;
	UD = ud;
	SD = sd;
	UB = ub;
	SB = sb;
	NUM_DEL_TYPES = ndt;
	NUM_BEN_TYPES = nbt;
	
	//reset some parameters that change when num mut types gets altered (it is default with one mut type)
	vector<double> tempVecDel(NUM_DEL_TYPES,1);
	vector<double> tempVecBen(NUM_BEN_TYPES,1);
	SD_F = tempVecDel;
	SD_DELTA_A = tempVecDel;
	SD_DELTA_B = tempVecDel;
	SD_PROB_A = tempVecDel;
}

//set parameters for fluctuating N
void World::setParFlucN(double nf, double nda, double ndb, double npa){
	N_F = nf;		
	N_DELTA_A = nda; 
	N_DELTA_B = ndb;
	N_PROB_A = npa;
}

//set parameters for fluctuating SDel
void World::setParFlucSDel(vector<double> sdf, vector<double> sdda, vector<double> sddb, vector<double> sdpa){
	SD_F = sdf;		
	SD_DELTA_A = sdda;
	SD_DELTA_B = sddb;
	SD_PROB_A = sdpa;
}

void World::setParMisc(double nsdf0, double dnsd0){
	N_SDEL0_F = nsdf0;
	D_N_SDEL0 = dnsd0;
}

void World::setParSave(string fdir, string fn){
	FILE_NAME = fn;
	FILE_DIR = fdir;
}

//create population and sequence of environments
void World::createWorld_V0(){
	
	//save newPopulation in class private variable Population
	Popul newPopulation(UD, SD, UB, SB, NUM_DEL_TYPES, NUM_BEN_TYPES, N);
	Population = newPopulation;

	//create sequence of env for N
	createStochFlucN(N_F, N_DELTA_A, N_DELTA_B, N_PROB_A);
	
	//create sequence of env for SDel
	createStochFlucSDel(SD_F, SD_DELTA_A, SD_DELTA_B, SD_PROB_A);
}

//create population and sequence of environments
//World_V1 allows for altering of covariance between fluctuations in DelMut type 0 and population size
//All other mutation types will fluctuate independent of DelMut0 and population size
void World::createWorld_V1(){
	
	//save newPopulation in class private variable Population
	Popul newPopulation(UD, SD, UB, SB, NUM_DEL_TYPES, NUM_BEN_TYPES, N);
	Population = newPopulation;
	
	//create sequence of env for N and SDel
	createStochFlucNSDel_V1(N_SDEL0_F, N_DELTA_A, N_DELTA_B, N_PROB_A, SD_F, SD_DELTA_A, SD_DELTA_B, SD_PROB_A, D_N_SDEL0);
}


//run the simulation using parameters given in Constructor
void World::runWorld(int replicate){

	REPLICATE = replicate;
	
	//start at generation 0
	Generation = 0;
	
	//run from generation 1 to NUM_Gen
	while(Generation < NUM_GEN){		
		
		//Add mutation to Population and recalculates fitness
		Population.doMutation(UD, sequenceSDel[Generation], UB, SB);
						
		//Selection
		Population.doSelection( sequenceN[Generation] );
		
		//display data in shell	
		if ( (Generation+1) % 100 == 0){
			displayData();
		}
		//increase gen counter
		Generation++;
	
		
		saveData();
		
	}
	
	
}

//display data in terminal
void World::displayData(){
	
	cout << "Gen = " << Generation+1 << "\n";
	cout << "N = " << Population.sizePop << "\n";
	cout << "sequenceN[Generation] = " << sequenceN[Generation] << "\n";
	
	cout << "sDel = ";
	for(int i = 0; i < NUM_DEL_TYPES; i++){
	cout << sequenceSDel[Generation][i] << ", ";
	}
	cout << "\n";
	
	cout << "UDel = ";
	for(int i = 0; i < NUM_DEL_TYPES; i++){
	cout << UD[i] << ", ";
	}
	cout << "\n";
	
	cout << "meanFitness: " << Population.meanFitness << "\n";
	cout << "meanDel_All: " << Population.meanDel_All << "\n";
	cout << "meanDel_Types: " << Population.meanDel_Types[0] << ", " << Population.meanDel_Types[1] << "\n";
	//cout << "meanBen_All: " << Population.meanBen_All << "\n";
	cout << "DelLLC_All: " << Population.DelLLC_All << " (" << Population.numDelLLC_All << ")" << "\n";
	cout << "DelLLC_Types: " << Population.DelLLC_Types[0] << ", " << Population.DelLLC_Types[1] << "\n";
	cout << "numDelLLC_Types: " << Population.numDelLLC_Types[0] << ", " << Population.numDelLLC_Types[1] << "\n";

	
	/*
	for(int i = 0; i < Population.numDelMutTypes; i++){	
		cout << "meanDel_Types[" << i << "]: "<< Population.meanDel_Types[i] << "\n";		
		cout << "DelLLC_Types[" << i << "]: "<< Population.DelLLC_Types[i] << " (" << Population.numDelLLC_Types[i] << ")" << "\n";
	}
	*/
	cout << "\n\n";	
}

void World::saveData(){
	
	int rep = REPLICATE;
	
	stringstream sstr;
	sstr << rep;
	//file name starts are "model results" and sstr (=whichRounds)
	string fileName(FILE_DIR);
	fileName.append(FILE_NAME);
	fileName.append(sstr.str());
	fileName.append(".txt");
	
	if (Generation == 1){
	
		ofstream modelResults;
		modelResults.open(fileName.c_str(), ios::out);
	
		modelResults 	<< "FILE_NAME,";
		modelResults    << "NUM_GEN,N,NUM_DEL_TYPES,NUM_BEN_TYPES,";
		
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< "UD[" << i << "],";
		}
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< "SD[" << i << "],";
		}
		for(int i = 0; i < NUM_BEN_TYPES; i++){
			modelResults	<< "UB[" << i << "],";
		}
		for(int i = 0; i < NUM_BEN_TYPES; i++){
			modelResults	<< "SB[" << i << "],";
		}
		
		modelResults	<< "N_SDEL0_F,D_N_SDEL0,";
		modelResults	<< "N_F,N_DELTA_A,N_DELTA_B,N_PROB_A,";
		
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< "SD_F[" << i << "],";
		}
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< "SD_DELTA_A[" << i << "],";
		}
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< "SD_DELTA_B[" << i << "],";
		}
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< "SD_PROB_A[" << i << "],";
		}
	
		modelResults 	<< "Gen,";
		modelResults 	<< "meanFitness,meanDel_All,DelLLC_All,numDelLLC_All,";
		modelResults	<< "meanDel_0,meanDel_1,";
		modelResults	<< "DelLLC_0,DelLLC_1,numDelLLC_0,numDelLLC_1" << "\n";
		
	
		modelResults.close();
	
	}
	

	// open the text file in which the data will be saved
	ofstream modelResults;
	modelResults.open(fileName.c_str(), ios::out | ios::app);
	
	if(modelResults.is_open()){
		// Write info modelResults
		modelResults	<< FILE_NAME << ",";
		modelResults    << NUM_GEN << "," << N << "," << NUM_DEL_TYPES << "," << NUM_BEN_TYPES << ",";
	
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< UD[i] << ",";
		}
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< SD[i] << ",";
		}
		for(int i = 0; i < NUM_BEN_TYPES; i++){
			modelResults	<< UB[i] << ",";
		}
		for(int i = 0; i < NUM_BEN_TYPES; i++){
			modelResults	<< SB[i] << ",";
		}
	
		modelResults	<< N_SDEL0_F << "," << D_N_SDEL0 << ",";
		modelResults	<< N_F << "," << N_DELTA_A << "," << N_DELTA_B << "," << N_PROB_A  << ",";
	
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< SD_F[i] << ",";
		}
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< SD_DELTA_A[i] << ",";
		}
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< SD_DELTA_B[i] << ",";
		}
		for(int i = 0; i < NUM_DEL_TYPES; i++){
			modelResults	<< SD_PROB_A[i] << ",";
		}
		
		modelResults 	<< Generation << ",";
		modelResults 	<< Population.meanFitness << "," << Population.meanDel_All << "," << Population.DelLLC_All << "," << Population.numDelLLC_All << ",";
		modelResults	<< Population.meanDel_Types[0] << "," << Population.meanDel_Types[1] << ",";
		modelResults	<< Population.DelLLC_Types[0] << "," << Population.DelLLC_Types[1] << "," << Population.numDelLLC_Types[0] << "," << Population.numDelLLC_Types[1] << "\n";
	}
	else{
		cout << "CANNOT OPEN modelResults FILE!!!!!" << "\n";	
	}
	modelResults.close();


}

void World::createStochFlucN(double nf, double nda, double ndb, double npa){
	
	//clear sequence created before (Ex. default settings)
	sequenceN.clear();
	
	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//initialize first env as env A
	double currentEnv = nda * N;
	
	sequenceN.push_back(currentEnv);

	for(int i = 0; i < NUM_GEN; i++){
		
		//Binomial variable: Fluctuate = {0,1} means {remain the same, fluctuate}
		//Using 1-fn because fn denotes amount of autocorrelation/no change; we want prob of success for changing env = 1-f
		double Fluctuate = gsl_ran_binomial (r, 1-nf, 1);
		
		//Env fluctuate when Fluctuate == 1
		//Update currentEnv to result of fluctuation (possible to fluctuate to same env)
		if( Fluctuate == 1){
		
			double UnifRV = gsl_ran_flat(r, 0.0, 1.0);	
			
			//if UnifRV <= npa (prob of env A), sequenceN in this generation is env A (modifier dan)
			//otherwise it is modifier dbn
			if( UnifRV <= npa ){
				currentEnv = nda * N;
			}
			else{
				currentEnv = ndb * N;
			}
		}
		//if Fluctuate = 0, currentEnv not change
		else{
			//do nothing if currentEnv not change
		}
		
		sequenceN.push_back( currentEnv );
		
		//adjust seed
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
	}
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
}



void World::createStochFlucSDel(vector<double> sdf, vector<double> sdda, vector<double> sddb, vector<double> sdpa){
	
	//clear sequence created before (Ex. default settings)
	sequenceSDel.clear();

	for(int i = 0; i < NUM_GEN; i++){
		vector<double> tempVec(NUM_DEL_TYPES, 0);
		sequenceSDel.push_back(tempVec);
	}

	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	//initialize first env as env A
	vector<double> currentEnv;
	multiplyVector_Double(currentEnv, SD, sdda);
	
	sequenceSDel[0] = currentEnv;

	for(int i = 0; i < NUM_GEN; i++){
		
		//allow all j mut types to fluctuate in selection
		for(int j = 0; j < NUM_DEL_TYPES; j++){
		
			//Binomial variable: Fluctuate = {0,1} means {remain the same, fluctuate}
			//Using 1-sdf because sdf denotes amount of autocorrelation/no change; we want prob of success for changing env = 1-sdf
			double Fluctuate = gsl_ran_binomial (r, 1-sdf[j], 1);
		
			//Env fluctuate when Fluctuate == 1
			//Update currentEnv to result of fluctuation (possible to fluctuate to same env)
			if( Fluctuate == 1){
		
				double UnifRV = gsl_ran_flat(r, 0.0, 1.0);	
			
				//if UnifRV <= sdpa (prob of env A), sequenceSDel in this generation is env A (modifier dan)
				//otherwise it is modifier dbn
				if( UnifRV <= sdpa[j] ){
					currentEnv[j] = sdda[j] * SD[j];
				}
				else{
					currentEnv[j] = sddb[j] * SD[j];
				}
			}
			//if Fluctuate = 0, currentEnv[j] not change
			else{
				//do nothing currentEnv[j] remains the same
			}
		}
		
		//add currentEnv to sequenceSDel[i]
		sequenceSDel[i] = currentEnv;
		
		//adjust seed
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
	}
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);
}







//create stochastic fluctuation for population size and sDel mutation 
//can adjust correlation in fluctuation between sDel0 and N 
//all other sDel types fluctuate independent of sDel0 and N
void World::createStochFlucNSDel_V1 ( double nsd0f, double nda, double ndb, double npa, vector<double> sdf, vector<double> sdda, vector<double> sddb, vector<double> sdpa , double dnsd0){
	
	//clear sequences created before
	sequenceN.clear();
	sequenceSDel.clear();
	
	//envRecord will keep record of env at each gen in regards to N and sDel0 (env defined below at envProb)
	//it is currently NEVER used, if I need to print it out or something, it can be used though
	vector<int> envRecord_V1;
	
	//let tempVec be selection for all mut types in env A
	vector<double> tempVec;
	multiplyVector_Double(tempVec, SD, sdda);
	
	//Initialize sequenceN and sequenceSDel as env A and envRecord as env 1 for all generations
	for(int i = 0; i < NUM_GEN; i ++){
		sequenceN.push_back(nda * N);
		sequenceSDel.push_back(tempVec);
		envRecord_V1.push_back(1);
	}
	
	
	//numEnv needs to be a const, so I can initialize arrays with size [NumEnv];
	const int			numEnv = 4;
	vector <double>		envProb(numEnv,0);
	vector <double> 	bin(numEnv,0);
	double 				tempNum = 0;
	double 				currentEnv;
	
	//Assign probabilities of each environment occurring
	envProb[0] = (npa)*(sdpa[0]) + dnsd0;
	envProb[1] = (npa)*(1-sdpa[0]) - dnsd0;
	envProb[2] = (1-npa)*(sdpa[0]) - dnsd0;
	envProb[3] = (1-npa)*(1-sdpa[0]) + dnsd0;
	
	//Calculate the boundary for the bins that will categorize environments along interval 0 to 1
	for(int i = 0; i < numEnv; i++){
		tempNum	+= envProb[i];
		bin[i] = tempNum;
	}
	
	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed);

	//Generate sequence of environments for NUM_GEN cycles from index 1 to NUM_GEN; index 0 already done above for initializing
	for(int i = 0; i < NUM_GEN; i++){

		for(int j = 0; j < NUM_DEL_TYPES; j++){
		
			//special case for sDel0 (i.e. j = 0)
			if(j == 0){
				//Fluctuate = {0,1} means {remain the same, fluctuate}
				//Using 1-nsdf because nsdf denotes amount of autocorrelation/no change; we want prob of success for changing env = 1-nsdf
				double Fluctuate = gsl_ran_binomial (r, 1-nsd0f, 1);
				
				//if i == 0, initial env chosen at random
				//Env fluctuate when Fluctuate == 1, set N and sDel according to env
				//save currentEnv to result of fluctuation (possible to fluctuate to same env)
				if( Fluctuate == 1 || i == 0 ){
		
					//Simulate multinomial sampling by picking uniform RV (between 0 and 1)
					//bins correspond to size of each environment, env determined by which bin uniform RV falls into
					double UnifRV = gsl_ran_flat(r, 0.0, 1.0);	
	
					if( UnifRV <= bin[0] ){
						sequenceN[i] = nda * N;
						sequenceSDel[i][0] = sdda[0] * SD[0];
						currentEnv = 1;
					}
					else if ( UnifRV > bin[0] && UnifRV <= bin[1] ){
						sequenceN[i] = nda * N;
						sequenceSDel[i][0] = sddb[0] * SD[0];
						currentEnv = 2;
					}
					else if ( UnifRV > bin[1] && UnifRV <= bin[2] ){
						sequenceN[i] = ndb * N;
						sequenceSDel[i][0] = sdda[0] * SD[0];
						currentEnv = 3;
					}
					else{
						sequenceN[i] = ndb * N;
						sequenceSDel[i][0] = sddb[0] * SD[0];
						currentEnv= 4;
					}
				}
				//Fluctuate == 0, means no change from prev gen
				else{
					sequenceN[i] = sequenceN[i-1];
					sequenceSDel[i][0] = sequenceSDel[i-1][0];
				}
				//Save sequence currentEnv in envRecord, in case I need it
				envRecord_V1[i] = currentEnv;
			}
			//sDel type != 0 fluctuate as normal (i.e. j != 0)
			//fluctuate according to sdf[] (NOT nsd0f)
			else{
				//Binomial variable: Fluctuate = {0,1} means {remain the same, fluctuate}
				double Fluctuate = gsl_ran_binomial (r, 1-sdf[j], 1);
				
				//if i == 0, initial env chosen at random
				//Env fluctuate when Fluctuate == 1, set sequenceSDel according to env
				if( Fluctuate == 1 || i == 0){
		
					double UnifRV = gsl_ran_flat(r, 0.0, 1.0);	
			
					//if UnifRV <= sdpa (prob of env A), sequenceSDel in this generation is env A (modifier dan)
					//otherwise it is modifier dbn
					if( UnifRV <= sdpa[j] ){
						sequenceSDel[i][j] = sdda[j] * SD[j];
					}
					else{
						sequenceSDel[i][j] = sddb[j] * SD[j];
					}
				}
				//if Fluctuate = 0, currentEnv[j] not change
				else{
					sequenceSDel[i][j] = sequenceSDel[i-1][j];
				}
			}
		}
		//adjust seed using time
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
	}
	
	//free memory; prevent memory leaks
	gsl_rng_free(r);	
}






//multiply elements of Vec1 and Vec2 in sequential order, give product to ProductVec
void World::multiplyVector_Double(vector<double> &ProductVec, vector<double> &Vec1, vector<double> &Vec2){

	if(Vec1.size() == Vec2.size()){
		
		vector<double> tempProductVec;
		
		for(int i = 0; i < Vec1.size(); i++){
			tempProductVec.push_back(Vec1[i] * Vec2[i]);
		}
		
		ProductVec = tempProductVec;
	}
	else{
		cout << "ERROR: Vectors are of different size" << "\n";
	}

}







/*Old constructor, not used anymore

//Constructor to create one population and store in private variable Population
World::World(int ng, int n, const vector<double> &ud, const vector<double> &sd, const vector<double> &ub, const vector<double> &sb, int ndt, int nbt){

	NUM_GEN = ng;
	N = n;
	UD = ud;
	SD = sd;
	UB = ub;
	SB = sb;
	NUM_DEL_TYPES = ndt;
	NUM_BEN_TYPES = nbt;
	
	Popul newPopulation(UD, SD, UB, SB, NUM_DEL_TYPES, NUM_BEN_TYPES, N);

	Population = newPopulation;
	
	createStochFlucN(0.95,1,0.01,0.5);
}
*/
