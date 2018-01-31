
#include "IncludeFiles.h"
#include "World.h"


//input int vector pointer, prints by dereferencing
void PrintVec(vector<int> *vVec)                                       
{
    vector<int>::iterator itI;                                          

	if(vVec->empty()==1){
	cout << "Vector is empty" << "\n";
	}
	else{
		for(itI=vVec->begin() ; itI < vVec->end(); itI++){
			cout << *itI << " ";                                            
		}
		cout << "\n";
		//cout << "Size of vector is: "<< vVec->size() << endl;
    }
}

//input double vector pointer, prints by dereferencing
void PrintVecDouble(vector<double> *vVec)
{
    vector<double>::iterator itI;

	if(vVec->empty()==1){
	cout << "Vector is empty" << "\n";
	}
	else{
		for(itI=vVec->begin() ; itI < vVec->end(); itI++){
			cout << *itI << " ";
		}
		cout << "\n";
		//cout << "Size of vector is: "<< vVec->size() << endl;
    }
}



// Main function for the program
int main(int argc, char* argv[])
{

	//*****-----SET PARAMETERS FOR SIMULATIONS-----*****

	//Set file directory
	string FileDir = "/Users/eddieho/Desktop/AsexMR_Model/Simulation/AsexMR_Test/Results/";
	//Set file date
	string FileDate = "Dec24";

	//Set replicate per set and Generation time
	int numReps = 500;
	int numGen = 100;

	//Set list of parameters to test
	double listU[] = {0.5};
	double listPU0[] = {0.87};
	double listSDel0[] = {0.003};
	double listSDel1[] = {0.4};
	double listN[] = {1,2,4,10,20,100};
	
	//*****-----SET PARAMETERS FOR SIMULATIONS-----*****
	
	//Store list of parameters, 2-D vector parListAll will store all parameter sets
	//parameter set holds values for these parameters: {totalU, propUdel0, sDel[0], sDel[1], N}
	vector <double> parList(5);
	vector < vector <double> > parListAll;
	
	//count num of parameter sets
	int numParSets = 0;
	
	//*****Generate list of parameter set to run model
	for(int a = 0; a < sizeof(listU) / sizeof(listU[0]); a++){
	
		parList[0] = listU[a];
	
		for(int b = 0; b < sizeof(listPU0) / sizeof(listPU0[0]); b++){
			
			parList[1] = listPU0[b];
		
			for(int c = 0; c < sizeof(listSDel0) / sizeof(listSDel0[0]); c++){
				
				parList[2] = listSDel0[c];
				
				for(int d = 0; d < sizeof(listSDel1) / sizeof(listSDel1[0]); d++){
					
					parList[3] = listSDel1[d];
					
					for(int e = 0; e < sizeof(listN) / sizeof(listN[0]); e++){
						
						//create new parameter set
						parList[4] = listN[e];
								
						//store new parameter set				
						parListAll.push_back(parList);
						//print new parameter set
						PrintVecDouble(&parListAll[numParSets]);
						//increase count
						numParSets++;
					}
				}
			}
		}
	}
	
	//display number of parameter sets
	cout << numParSets << "\n";
	cout << parListAll.size() << "\n";
	//*****Generate list of parameter set to run model

	//*****Save parListAll into a .csv file
	string parListFileName(FileDir);
	parListFileName.append(FileDate);
	parListFileName.append("_parListAll.csv");
	
	ofstream parListFile;
	parListFile.open(parListFileName.c_str(), ios::out);
	
	//Print file headers
	parListFile		<< "FileDate,";
	parListFile		<< "parSet,";
	parListFile		<< "numReps,numGen,";
	parListFile		<< "totalUDel,propUdel0,";
	parListFile		<< "UD[0],UD[1],";
	parListFile		<< "SD[0],SD[1],";
	parListFile		<< "N";
	parListFile		<< "\n";
	
	//Print all parameter sets
	for(int n = 0; n < numParSets; n++){
	
		parListFile		<< FileDate << ",";
		parListFile		<< n << ",";
		parListFile		<< numReps << "," << numGen << ",";
		parListFile		<< parListAll[n][0] << "," << parListAll[n][1] << ",";
		parListFile		<< parListAll[n][0]*parListAll[n][1] << "," << parListAll[n][0]*(1-parListAll[n][1]) << ",";
		parListFile		<< parListAll[n][2] << "," << parListAll[n][3] << ",";
		parListFile		<< parListAll[n][4];
		parListFile		<< "\n";
	
	}
	//close file
	parListFile.close();
	//*****Save parListAll into a .csv file
	
	
	//Run simulation, looping through all parameter sets in parListAll[][]
  	for(int n = 0; n < numParSets; n++){	  
	  
		//Set number of types of del and ben mutations
		int numDelMutTypes = 2;
		int numBenMutTypes = 1;
	
		vector <double> UDel(numDelMutTypes);
		vector <double> sDel(numDelMutTypes);
		vector <double> UBen(numBenMutTypes);
		vector <double> sBen(numBenMutTypes);
	
		//Set mut, sel parameters for del and ben mutations
		double totalUDel = parListAll[n][0];
		double propUdel0 = parListAll[n][1];
	
		UDel[0] = propUdel0*totalUDel;	
		sDel[0] = parListAll[n][2];

		UDel[1] = (1-propUdel0)*totalUDel;
		sDel[1] = parListAll[n][3];
	
		//Set PopSize
		int PopSize = parListAll[n][4];
	
	
		//beneficial mutations
		UBen[0] = 0;
		sBen[0] = 0;
	
		//*****Parameters to alter fluctuating N and S
		double NF = 1;
		double NDeltaA = 1;
		double NDeltaB = 1;
		double NProbA = 1;
	
		vector<double> sDelF(numDelMutTypes);
		vector<double> sDelDeltaA(numDelMutTypes);
		vector<double> sDelDeltaB(numDelMutTypes);
		vector<double> sDelProbA(numDelMutTypes);
	
		sDelF[0] = 1;
		sDelF[1] = 1;
	
		sDelDeltaA[0] = 1;
		sDelDeltaA[1] = 1;

		sDelDeltaB[0] = 1;
		sDelDeltaB[1] = 1;
	
		sDelProbA[0] = 1;
		sDelProbA[1] = 1;
	
		//createWorld_V1 uses NSDel0F for autocorrelation of N and sDel0, not NF or sDelF[0]
		double NSDel0F = 0;
		double DNSDel0 = 0;
		//*****Parameters to alter fluctuating N and S
	
		
		//Construct file name as: FileDate_n_; replicate number will be added within runWorld->daveData		
		stringstream sstr1;
		sstr1 << n;

		string tempFileName(FileDate);
		tempFileName.append("_");
		tempFileName.append(sstr1.str());
		tempFileName.append("_");
	
		string FileName = tempFileName;
		
		
		World newWorld;
	
		newWorld.setParBasic(numGen, PopSize, UDel, sDel, UBen, sBen, numDelMutTypes, numBenMutTypes);
	
		newWorld.setParFlucN(NF, NDeltaA, NDeltaB, NProbA);
	
		newWorld.setParFlucSDel(sDelF, sDelDeltaA, sDelDeltaB, sDelProbA);
	
		newWorld.setParMisc(NSDel0F, DNSDel0);
	
		newWorld.setParSave(FileDir, FileName);
		
		for(int i = 0; i < numReps; i++){
			newWorld.createWorld_V0();
			newWorld.runWorld(i);
		}
	

	}
	
	
	
	
	
	/*
	vector<int> testVec;
	testVec.push_back(1);
	PrintVec(&testVec);
	
	testVec.clear();
	testVec.push_back(2);
	PrintVec(&testVec);
	
	
	vector< vector<double> > sequenceSd;
	
	for(int i = 0; i < 10; i++){
		vector<double> tempVec(5, 0);
		sequenceSd.push_back(tempVec);
		
		PrintVecDouble(&sequenceSd[i]);
	}
	
	
	
	vector<double> tempVec(5, 0);
	
	for(int i = 0; i < 10; i++){
		
		for(int j = 0; j < 5; j++){
		
			tempVec[j] = 1;
		
		}
		
		sequenceSd[i]=tempVec;
		
		PrintVecDouble(&sequenceSd[i]);
	}
	*/
	
	return 0;
}


