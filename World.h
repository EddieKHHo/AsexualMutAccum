#include "Popul.h"


class World{


	private:
		
		//Population that will be run in simulation
		Popul Population;
		
		int Generation;

		vector<double> sequenceN;
		
		vector< vector<double> > sequenceSDel;
		
		//parameters for running simulation
		int REPLICATE;
		
		int NUM_GEN;
		int N;
		vector<double> UD;
		vector<double> SD;
		vector<double> UB;
		vector<double> SB;
		int NUM_DEL_TYPES;
		int NUM_BEN_TYPES;
		
		//parameters for fluctuating N
		double N_F;		
		double N_DELTA_A, N_DELTA_B;
		double N_PROB_A;
		//parameters for fluctuating SDel
		vector<double> SD_F;		
		vector<double> SD_DELTA_A, SD_DELTA_B;
		vector<double> SD_PROB_A;
		
		//correlation in fluctuations between N and sDel mut type 0
		double N_SDEL0_F;
		double D_N_SDEL0;
		
		int SAVE_DATA_INTERVAL;
		
		string FILE_DIR;
		string FILE_NAME;
		
		
	
	public:
		
		World();
		
		void setParBasic(int ng, int n, const vector<double> &ud, const vector<double> &sd, const vector<double> &ub, const vector<double> &sb, int ndt, int nbt);
		
		void setParFlucN(double nf, double nda, double ndb, double npa);
		
		void setParFlucSDel(vector<double> sdf, vector<double> sdda, vector<double> sddb, vector<double> sdpa);
		
		void setParMisc(double nsd0f, double dnsd0);
		
		void setParSave(string fdir, string fn);
		
		void createWorld_V0();
		
		void createWorld_V1();
		
		void runWorld(int replicate);
		
		void displayData();
		
		void saveData();
		
		
		
		
		void createStochFlucN(double nf, double nda, double ndb, double npa);
		
		void createStochFlucSDel(vector<double> sdf, vector<double> sdda, vector<double> sddb, vector<double> sdpa);
		
		void createStochFlucNSDel_V1 (double nsdf, double nda, double ndb, double npa, vector<double> sdf, vector<double> sdda, vector<double> sddb, vector<double> sdpa , double dnsd0);
		
		
		void multiplyVector_Double(vector<double> &ProductVec, vector<double> &Vec1, vector<double> &Vec2);
		
		//old constructor, not used anymore
		//World(int ng, int n, const vector<double> &ud, const vector<double> &sd, const vector<double> &ub, const vector<double> &sb, int ndt, int nbt);

};