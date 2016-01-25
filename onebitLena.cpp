// Math stuff
#include <random>
#include <math.h>
#include <cmath>
// I don't know what this is for
#include <stdlib.h>
// For reading in csvs
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#define ARRAYDIMX 512
#define ARRAYDIMY 512

#define SIGMA 100

//uses brace initialization, so must use the C++11 compiler. 
// call as g++ -std=c++11 file-name.pp -o file-handle
class Sample
{
public:
	// void set_d(int dim){
	// 	d=dim;
	// }
	Sample():
	mt(std::random_device()()),
	G(0,SIGMA)
	{}
	// void draw_kern(int*, int*);
	int draw_bern(double);
	double draw_gauss();
	// Sample();
	// ~Sample();
private:
	// int d=DIM;//dimension of sample space

	// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	// std::default_random_engine generator; //(seed);

  	// std::random_device rd;
  	std::mt19937 mt;
	std::normal_distribution<double> G;
};//add constructor and destructor

// void Sample::draw_kern(int* x, int* y){
// 	int z=0;
// 	for (int i=0; i<DIM; i++){
// 		z=round(G(mt))+x[i];
// 		if (z>SPACE_SIZE){
// 			z=2*SPACE_SIZE-z;
// 		}else if(z<0){
// 			z=-z;
// 		}
// 		y[i]=z;
// 	}//this implements an uncorrelated multivar DIM gaussian 
// 	//with zero mean and sigma*I correlation matrix. This is comp.
// 	//bottleneck
// }

int Sample::draw_bern(double param){
	std::discrete_distribution<int> B {param, 1-param};

	return B(mt);
}

double Sample::draw_gauss(){
	return G(mt);
}

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

std::istream& operator>>(std::istream& str,CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   

void writeDataInt(int ** data, int row, int col){
	std::ofstream dataFile;
	dataFile.open("foo2.csv");
	for (int i=0; i<row; i++){
		for (int j=0; j<col-1; j++){
			if(data[i][j]==1){
				dataFile<<255<<",";
			}
			else{
				dataFile<<0<<",";
			}
		}
		if(data[i][col-1]==1){
			dataFile<<255<<"\n";
		}
		else{
			dataFile<<0<<"\n";
		}
	}
	dataFile.close();
}

void readDataInt(int ** data, int numrows, int numcols){
	std::ifstream file("foo.csv");
	CSVRow row;
	int i=0;
	int val;
	while(file >> row && i<numrows)
    {
    	for (int j=0; j<numcols; j++){
		    val = atoi(row[j].c_str());
		    if(val>127){
		    	data[i][j] = 1;
		    }
		    else{
		    	data[i][j] = -1;
			}
		}
		i++;
    }
}

void addNoise(int ** data, int row, int col, double p){
	Sample samp;
	for (int i=0; i<row; i++){
		for (int j=0; j<col; j++){
			//with prob p flip the bit
			if(!samp.draw_bern(p))
				data[i][j] = -1 * data[i][j];
		}
	}
}

int main(){
	double p = 0.2;
	int** data=new_MatrixI(ARRAYDIMX,ARRAYDIMY);
	readDataInt(data, ARRAYDIMX, ARRAYDIMY);
	addNoise(data, ARRAYDIMX, ARRAYDIMY,p);
	writeDataInt(data, ARRAYDIMX, ARRAYDIMY);
	return 0;
}
