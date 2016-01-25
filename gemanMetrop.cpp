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

#define END_POINT 5 //interval [-END_POINT, END_POINT]
#define DIM 2
#define DELTA 0.01
// #define SPACE_SIZE 2*END_POINT/DELTA+1
#define SPACE_SIZE 1001
#define ITERS 3000
#define SIGMA 100
#define CYCLES 1
#define PI 3.14159265
#define XSTART 1.75
#define YSTART 2.5
#define TEMPCOEFF 1/10

#define IMG_WIDTH 512
#define IMG_LENGTH 512
#define NHBD_SIZE 4
#define IsingBETA 1
#define Ising_h 1
#define Ising_p 0.1

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


int mod(int n, int m){
	if (n>=0) return n%m;
	else return n%m+m;
}

class Sample
{
public:
	Sample():
	mt(std::random_device()()),
	G(0,SIGMA)
	{}
	void draw_kern(int*, int*);
	int draw_bern(double);
	void randPerm(int*, int);
	// void vising_scheme(int* orderI, int* orderJ){
	// int maxI=IMG_LENGTH;
	// int maxJ=IMG_WIDTH;

	// 	// for(int i=0; i<maxI; i++){
	// 	// 	draw_bern(1/(maxI-i))
	// 	// }
	// }
private:

  	std::mt19937 mt;
	std::normal_distribution<double> G;
};//add constructor and destructor

void Sample::draw_kern(int* x, int* y){
	int z=0;
	for (int i=0; i<DIM; i++){
		z=round(G(mt))+x[i];
		if (z>SPACE_SIZE){
			z=2*SPACE_SIZE-z;
		}else if(z<0){
			z=-z;
		}
		y[i]=z;
	}//this implements an uncorrelated multivar DIM gaussian 
	//with zero mean and sigma*I correlation matrix. This is comp.
	//bottleneck
}

int Sample::draw_bern(double param){
	std::discrete_distribution<int> B {param, 1-param};
	return B(mt);
}

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

class Ising//modelling the image on a torus.
{
public:
	double H_loc(int**, int**, int, int, double p);
	void nhbd(int sI, int sJ, int* site_nhbdsI, int* site_nhbdsJ){
		//implements a NN neighborhood structure
		//on a torus.
		site_nhbdsI[0]=sI;
		site_nhbdsI[1]=mod(sI-1,IMG_LENGTH);
		site_nhbdsI[2]=sI;
		site_nhbdsI[3]=mod(sI+1,IMG_LENGTH);

		site_nhbdsJ[0]=mod(sJ-1,IMG_WIDTH);
		site_nhbdsJ[1]=sJ;
		site_nhbdsJ[2]=mod(sJ+1,IMG_WIDTH);
		site_nhbdsJ[3]=sJ;
	}
private:
	int* s_nhbdsI=new int[NHBD_SIZE];
	int* s_nhbdsJ=new int[NHBD_SIZE];

	// void nhbd(int sI, int sJ, int* site_nhbdsI, int* site_nhbdsJ){
	// 	//implements a NN neighborhood structure
	// 	//on a torus.
	// 	site_nhbdsI[0]=sI;
	// 	site_nhbdsI[1]=(sI-1)%IMG_WIDTH;
	// 	site_nhbdsI[2]=sI;
	// 	site_nhbdsI[3]=(sI+1)%IMG_LENGTH;

	// 	site_nhbdsJ[0]=(sJ-1)%IMG_WIDTH;
	// 	site_nhbdsJ[1]=sJ;
	// 	site_nhbdsJ[2]=(sJ+1)%IMG_LENGTH;
	// 	site_nhbdsJ[3]=sJ;
	// }
};

double Ising::H_loc(int** config, int** obs, int sI, int sJ, double p){
	nhbd(sI, sJ, s_nhbdsI, s_nhbdsJ);
	int nhbd_sum=0;
	int site=config[sI][sJ];
	int obs_site=obs[sI][sJ];

	for (int i=0; i<NHBD_SIZE; i++){
		for (int j=0; j<NHBD_SIZE; j++){
			nhbd_sum+=site*config[s_nhbdsI[i]][s_nhbdsJ[j]];
		}
	}
	return -IsingBETA*nhbd_sum-Ising_h*site-0.5*log(p/(1-p))*site*obs_site;
}

class Metropolis//all internal mechanics are done in integers
//this is data abstraction.
{
public:
	void set_initVal(double* init){
		x0=init;
		double_to_index(init,i0);
	}
	void double_to_index(double * x, int* ind){
		for (int i=0; i<DIM; i++){
			ind[i]=(x[i]+END_POINT)/DELTA;
		}
	}
	void index_to_double(int* ind, double *x){
		for (int i=0;i<DIM; i++){
			x[i]=ind[i]*DELTA-END_POINT;
		}
	}
	void metSamp(int**, int**);
private:
	int iters;
	int* i0=new int[DIM];
	double* x0=new double[DIM];
	// int **chain;
	Sample sampler;
	Ising IModel;

	// double objFunc(double* y){
	double objFunc(double* x){
		// double x=y-1;
		double out=0;
		// double x[1];
		// x[0]=y[0]-1;
		// out=(x[0])*(x[0])+x[1]*x[1]+sin(3*x[0])+sin(4*x[1]);

		//out=sin((2*PI/3.5)*x[0])*cos((2*PI/3.5)*x[1]);//Egg carton
		// out=sin((2*PI/3.5)*x[0])*cos((2*PI/3.5)*x[1])+0.1*(x[0]*x[0]+x[1]*x[1]);
		out=5*sin((2*PI/3.5)*x[0])*cos((2*PI/3.5)*x[1])+0.1*(x[0]*x[0]+x[1]*x[1]);
		// out=5*sin(2*PI/(0.75*3.5)*x[0])+5*sin(2*PI/(2*3.5)*x[0])+5*sin(2*PI/(4*3.5)*x[0])+0.1*x[0]*x[0];

		// for(int i=0; i<DIM; i++){
		// 	out+=x[i]*x[i];
		// }//x^Tx
		return out;
	}
};

void Metropolis::metSamp(int** obs, int** res){
	printf("entering metSamp\n");
	// set_eq_I(chain[0],i0);
	// chain[0]=i0;
	// int s0[2]={0,0};
	int* s_I=new int[IMG_WIDTH];
	int* s_J=new int[IMG_LENGTH];

	for (int i=0; i<IMG_WIDTH; i++){
		s_I[i]=i;
	}
	for (int i=0; i<IMG_LENGTH; i++){
		s_J[i]=i;
	}

	sampler.randPerm(s_I, IMG_WIDTH);
	sampler.randPerm(s_J, IMG_LENGTH);
	//Make vising scheme:

	// printf("vising scheme done\n");

	int s0[DIM]={s_I[0], s_J[0]};

	int *sNow=new int[DIM];
	// int *sNext=new int[DIM]; 
	// double *xNow=new double[DIM];
	// double *xNext=new double[DIM]; 
	// double *xTest=new double[DIM]; 
	int A;
	double hNow, hNext;
	// int** config=new_MatrixI(IMG_LENGTH, IMG_WIDTH);
	// int** propConfig=new_MatrixI(IMG_LENGTH, IMG_WIDTH);

	for(int i=0; i<IMG_LENGTH; i++){
		for(int j=0; j<IMG_WIDTH; j++){
			// set_eq_I(config[i][j], obs[i][j]);
			// config[i][j]=obs[i][j];
			res[i][j]=obs[i][j];
		}
	}
	set_eq_I(sNow, s0);
	// index_to_double(sNow,xNow);
	// fNow=objFunc(xNow);
	// printf("set to enter restoration loop\n");

	hNow=IModel.H_loc(res, obs, sNow[0],sNow[1], Ising_p);
	// printf("H.loc is good\n");

	// printf("hNow=, %lf\n", hNow);
	for(int i=1; i<ITERS; i++){
		for(int k=0; k<IMG_WIDTH; k++){
			for(int q=0; q<IMG_LENGTH; q++){
				// printf("k=%d, q=%d\n", k, q);
				// sampler.draw_kern(sNow,sNext);
				// sampler.draw_site(sNow, sNext);//pick new site

				// index_to_double(sNext,xNext);

				// fNext=objFunc(xNext);
				res[sNow[0]][sNow[1]]*=-1;//flip bit
				// printf("bit fip is good\n");
				hNext=IModel.H_loc(res, obs, sNow[0], sNow[1], Ising_p);//calculate potential

				// printf("nHext=%lf\n", hNext);

				if (hNext<=hNow){
					hNow=hNext;
					// printf("keeping\n");
					// set_eq_I(sNow, sNext);
				}else{
					A=sampler.draw_bern(exp(-(log(1+i)*TEMPCOEFF)*(hNow-hNext)));
					// printf("drawing\n");
					if(A==1){
						hNow=hNext;
						// printf("keeping\n");
						// set_eq_I(sNow, sNext);
					} else{
						res[sNow[0]][sNow[1]]*=-1;
						// printf("not keeping\n");
					}
				}
				// set_eq_I(chain[i],sNow);
				// printf("looping\n");

				sNow[0]=s_I[k];
				sNow[1]=s_J[q];

				// printf("new points: i=%d, j=%d\n\n", sNow[0], sNow[1]);
			}
		}
	}
	printf("done\n");
	// res=res;
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
	std::ifstream file("foo2.csv");
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

// write matrix into a CSV
void writeDataInt(int ** data, int row, int col){
	std::ofstream dataFile;
	dataFile.open("restoration.csv");
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


int main(){
	// Declare variables
	int** observation=new_MatrixI(IMG_WIDTH, IMG_LENGTH);
	int** result=new_MatrixI(IMG_WIDTH, IMG_LENGTH);

	// load csv data into observation
	readDataInt(observation, IMG_WIDTH, IMG_LENGTH);
	// initialize result = observation
	set_eq_I(result, observation);

	// instantiate metropolis sampler
	Metropolis met;
	// run metropolis sampler
	met.metSamp(observation, result);

	// write result to csv
	writeDataInt(result, IMG_WIDTH,IMG_LENGTH);
	return 0;
}
