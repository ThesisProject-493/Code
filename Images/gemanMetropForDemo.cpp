// stuff for metropolis
#include <algorithm>
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
#define ITERS 1000
#define SIGMA 40
#define TWOSIGMASQUARED 1600

// inverse temp
#define TEMPCOEFF 0.3

// cutoff for truncated parabola
#define ALPHA 7000

// higher lambda = less fidelity to observation data
#define LAMBDA 2

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

double localEnergy(int** observation, int** result, int val, int i, int j){
	double energy;
	energy = (val - observation[i][j])*(val - observation[i][j]);
	// energy = energy + (val-result[i][(j-1+IMG_LENGTH)%IMG_LENGTH])*(val-result[i][(j-1+IMG_LENGTH)%IMG_LENGTH]);
	// energy = energy + (val-result[i][(j+1)%IMG_LENGTH])*(val-result[i][(j+1)%IMG_LENGTH]);
	// energy = energy + (val-result[(i-1+IMG_WIDTH)%IMG_WIDTH][j])*(val-result[(i-1+IMG_WIDTH)%IMG_WIDTH][j]);
	// energy = energy + (val-result[(i+1)%IMG_WIDTH][j])*(val-result[(i+1)%IMG_WIDTH][j]);
	// min(squared distance,alpha)
	energy = energy + std::min(LAMBDA*(val-result[i][(j-1+IMG_LENGTH)%IMG_LENGTH])*(val-result[i][(j-1+IMG_LENGTH)%IMG_LENGTH]), ALPHA);
	energy = energy + std::min(LAMBDA*(val-result[i][(j+1)%IMG_LENGTH])*(val-result[i][(j+1)%IMG_LENGTH]), ALPHA);
	energy = energy + std::min(LAMBDA*(val-result[(i-1+IMG_WIDTH)%IMG_WIDTH][j])*(val-result[(i-1+IMG_WIDTH)%IMG_WIDTH][j]), ALPHA);
	energy = energy + std::min(LAMBDA*(val-result[(i+1)%IMG_WIDTH][j])*(val-result[(i+1)%IMG_WIDTH][j]), ALPHA);
	energy = energy / TWOSIGMASQUARED;
	// absolute distance
	// energy = energy + std::min(LAMBDA*abs(val-result[i][(j-1+IMG_LENGTH)%IMG_LENGTH]), ALPHA);
	// energy = energy + std::min(LAMBDA*abs(val-result[i][(j+1)%IMG_LENGTH]), ALPHA);
	// energy = energy + std::min(LAMBDA*abs(val-result[(i-1+IMG_WIDTH)%IMG_WIDTH][j]), ALPHA);
	// energy = energy + std::min(LAMBDA*abs(val-result[(i+1)%IMG_WIDTH][j]), ALPHA);
	return energy;
}

void metSamp(int** observation, int** result){
	// define variables
	int i, j, k, val;
	double proposedEnergy, currentEnergy, param;

	// vector for randomly visitng sites
	std::vector<int> indices;
	for (k=0; k<IMG_WIDTH*IMG_LENGTH; k++){
		indices.push_back(k);
	}
	// machinery for randomness
	std::random_device random_dev;
	std::mt19937       generator(random_dev());
	std::uniform_int_distribution<> dis(0, 255);

	for(int sweep=2; sweep<ITERS; sweep++){
		// permute indices
		std::shuffle(indices.begin(), indices.end(), generator);
		// perform sweep
		for(int k=0; k<IMG_WIDTH*IMG_LENGTH; k++){
			// pick next site
			j = indices[k] % IMG_LENGTH;
			i = (indices[k] - j) / IMG_WIDTH;

			// pick a proposal
			val = dis(generator);

			// calculate energies
			currentEnergy = localEnergy(observation, result, result[i][j], i, j);
			proposedEnergy = localEnergy(observation, result, val, i, j);

			// acceptance stage (after denial, anger, bargaining, and depression)
			if (proposedEnergy<=currentEnergy){
				// keep new pixel value
				result[i][j] = val;
			}
			else{
				param = exp(-log(sweep)*TEMPCOEFF*(proposedEnergy-currentEnergy));
				std::discrete_distribution<int> B {1-param, param};
			 	if(B(generator))
					result[i][j] = val;
			}
		}
	}
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
	std::ifstream file("noisyData40.csv");
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
	dataFile.open("restoredData40_13.csv");
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
	for (int i=0; i<IMG_WIDTH; i++){
		for (int j=0; j<IMG_LENGTH; j++){
			result[i][j] = observation[i][j];
		}
	}

	// run metropolis sampler
	metSamp(observation, result);

	// write result to csv
	writeDataInt(result, IMG_WIDTH,IMG_LENGTH);
	return 0;
}
