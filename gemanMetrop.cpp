// stuff for metropolis
#include <chrono>
#include <random>
#include <functional>
#include <array>
#include <math.h>
#include <cmath>
#include <time.h>

// file I/O stuff
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// hyperparameters
#define ITERS 10
#define SIGMA 20
#define TEMPCOEFF 1/10

// image dimensions
// Luckily they are the same. I think I have them confused throughout -G
#define IMG_WIDTH 512
#define IMG_LENGTH 512


//uses brace initialization, so must use the C++11 compiler. 
// call as g++ -std=c++11 file-name.pp -o file-handle

// functions for creating matrices
int **new_MatrixI(int row, int col){
		int ** matrix=new int *[row];
		int i,j;
		for(i=0; i<row; i++){
			matrix[i]=new int [col];
		}
		return (matrix);
	}
double **new_MatrixD(int row, int col){
	double ** matrix=new double *[row];
	int i,j;
	for(i=0; i<row; i++){
		matrix[i]=new double [col];
	}
	return (matrix);
}

// functions for deleting matrices (unused, I believe -G)
void del_MatrixI(int **mat, int row){
	int i;
	for(i=0; i<row; i++){
		delete[] mat[i];
	}
	delete[] mat;
}
void del_MatrixD(double **mat, int row){
	int i;
	for(i=0; i<row; i++){
		delete[] mat[i];
	}
	delete[] mat;
}

// functions for copying matrices
void set_eq_I(int* lhs, int* rhs){
	for (int i=0; i<DIM; i++){
		lhs[i]=rhs[i];
	}
}
void set_eq_D(double* lhs, double* rhs){
	for (int i=0; i<DIM; i++){
		lhs[i]=rhs[i];
	}
}

// Does sampling from distributions
class Sample
{
public:
	Sample():
	mt(std::random_device()()),
	G(0,SIGMA)
	{}
	int draw_bern(double);
	void randPerm(int*, int);
	double draw_gauss();
private:
  	std::mt19937 mt;
	std::normal_distribution<double> G;
}; //add constructor and destructor

// returns zero mean gaussian of variance SIGMA
double Sample::draw_gauss(){
	return G(mt);
}

// returns 1 with prob param
int Sample::draw_bern(double param){
	std::discrete_distribution<int> B {1-param, param};
	return B(mt);
}

// permutes a vector
void Sample::randPerm(int* vec, int len){
	// int index=len-1;
	srand(time(NULL));

	int temp, ind;

	for (int i=0; i<len; i++){
		ind=std::rand()%(len-i);
		temp=vec[len-i-1];
		vec[len-i-1]=vec[ind];
		vec[ind]=temp;
	}
}

void metSamp(int** observation, int** result){
	// define variables
	int i;
	int j;
	int val;

	// instantiate sampler
	Sample sampler;

	// randomization stuff (might need re-jigging)
	// int* s_I=new int[IMG_WIDTH];
	// int* s_J=new int[IMG_LENGTH];
	// for (int i=0; i<IMG_WIDTH; i++){
	// 	s_I[i]=i;
	// }
	// for (int i=0; i<IMG_LENGTH; i++){
	// 	s_J[i]=i;
	// }
	// sampler.randPerm(s_I, IMG_WIDTH);
	// sampler.randPerm(s_J, IMG_LENGTH);

	for(int sweep=2; sweep<ITERS; sweep++){
		for(int k=0; k<IMG_WIDTH*IMG_LENGTH; k++){
			// pick next site


			// calculate energies
			currentEnergy = localEnergy(observation, result, result[i][j], i, j);
			proposedEnergy = localEnergy(observation, result, val, i, j);

			if (proposedEnergy<=currentEnergy){
				// keep new pixel value
				// update energy for whole clique
			}
			else if(sampler.draw_bern(exp(-log(sweep)*TEMPCOEFF*(proposedEnergy-currentEnergy)))){
				// keep new pixel value
				// update energy for whole clique
			}
		}
	}
}

double localEnergy(int** observation, int** result, int val, int i, int j){
	return 0.5;
}

// Used for data I/O
class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str,line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }
    private:
        std::vector<std::string>    m_data;
};

// Used for data I/O
std::istream& operator>>(std::istream& str,CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   

// Read data from CSV into matrix
void readDataInt(int ** data, int numrows, int numcols){
	std::ifstream file("foo.csv");
	CSVRow row;
	int i=0;
	while(file >> row && i<numrows)
    {
    	for (int j=0; j<numcols; j++){
			
		    data[i][j] = atoi(row[j].c_str());
		}
		i++;
    }
}

// write matrix into a CSV
void writeDataInt(int ** data, int row, int col){
	std::ofstream dataFile;
	dataFile.open("gaussianObs.csv");
	for (int i=0; i<row; i++){
		for (int j=0; j<col-1; j++){
			dataFile<<data[i][j]<<",";
		}
		dataFile<<data[i][col-1]<<"\n";
	}
	dataFile.close();
}

int main(){
	// Declare variables
	int** observation=new_MatrixI(IMG_WIDTH, IMG_LENGTH);
	int** result=new_MatrixI(IMG_WIDTH, IMG_LENGTH);

	// load csv data into observation
	readDataInt(observation, IMG_WIDTH, IMG_LENGTH);
	// initialize result = observation
	set_eq_I(result, observation);

	// run metropolis sampler
	metSamp(observation, result);

	// write result to csv
	writeDataInt(result, IMG_WIDTH,IMG_LENGTH);
	return 0;
}
