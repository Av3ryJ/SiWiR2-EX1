#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "grid.hpp"

using namespace std;

double l2_norm(double *data, double *rhs, int gridsize, double stepsize) {
    double sum = 0;
    double alpha = 4.0/(stepsize*stepsize);
    double beta = 1.0/(stepsize*stepsize);
    for (int row = 1; row < gridsize-1; row++) {
        for (int col = 1; col < gridsize-1; col++) {
            double residual = rhs[row*gridsize + col] - alpha*data[row*gridsize + col]
                              +beta*data[(col-1) + row*gridsize]
                              +beta*data[(col+1) + row*gridsize]
                              +beta*data[col + (row-1)*gridsize]
                              +beta*data[col + (row+1)*gridsize];
            sum += residual*residual;
        }
    }
    int inner_points = gridsize-2;
    double factor = 1.0/(inner_points*inner_points); 
    return sqrt(factor*sum);
}
void calculate_residuum(double *res, double *data, double *rhs, int gridsize, int stepsize){
    double alpha = 4.0/(stepsize*stepsize);
    double beta = 1.0/(stepsize*stepsize);
    for (int row = 1; row < gridsize-1; row++) {
        for (int col = 1; col < gridsize-1; col++) {
            res[row*gridsize+col] = rhs[row*gridsize + col] - alpha*data[row*gridsize + col]
                              +beta*data[(col-1) + row*gridsize]
                              +beta*data[(col+1) + row*gridsize]
                              +beta*data[col + (row-1)*gridsize]
                              +beta*data[col + (row+1)*gridsize];
        }
    }
} 

void gauss_seidel_smoother(double *data,double *rhs, int gridsize, double stepsize){ 
    double alpha = 4.0/(stepsize*stepsize);
    double beta = 1.0/(stepsize*stepsize);
    for (int row = 1; row < gridsize-1; ++row) { // row is y col is x
        for (int col = 1; col < gridsize-1; ++col) {    // nicht ueber Rand iterieren
                data[col + row*gridsize] = (1/alpha)*
                        (rhs[col+row*gridsize]
                        +beta*data[(col-1) + row*gridsize]
                        +beta*data[(col+1) + row*gridsize]
                        +beta*data[col + (row-1)*gridsize]
                        +beta*data[col + (row+1)*gridsize]);
        }
    }
}

void addInterpolation(double *data,double *interpol, int gridsize){
     for (int row = 1; row < gridsize-1; ++row) { 
        for (int col = 1; col < gridsize-1; ++col) {
            data[col+row*gridsize] += interpol[col+row*gridsize]; 
        }
    }    
}

 void multigrid(double *u, double *rhs, int lvl){

    int gridsize = (int) pow(2,lvl) + 1;
    double stepsize = 1./(gridsize-1);

    // Presmoothing of Au = b -->uh = Gauss-Seidel(uh, f, h)
    gauss_seidel_smoother(u, rhs, gridsize, stepsize);

    //Compute Residual r = b-Auh 
     Grid *res = new Grid(lvl, true);
     calculate_residuum(res->data_, u, rhs, gridsize, stepsize);

    //restrict residual to level l-1  
    Grid *restricted = res->restrict();

    //leeres Gitter anlegen eps--> so groß wie restricted size //set u=0 on l-1 --> auch Randbedingungen werden 0
    Grid *eps = new Grid(lvl-1, false); //hier muss alles 0 sein!

    // Abbruchbedingung Rekursion: wenn kleinste Gittergröße erreicht ist
    //sonst: eps = V-Cycle(leeres Gitter eps, restricted Gitter, 2*h) --> Rekursion
    if(lvl==1) {
        eps->data_ = u; //nicht sicher 
    } else {
        multigrid(eps->data_, restricted->data_, lvl-1);
    } 
    // Ab hier wieder aufwärts in der Rekurison 
    //TODO Funktion zur Addition 
    Grid *eps2 = eps->interpolate();
    //correction
    addInterpolation(u, eps2->data_, gridsize);

    // postsmoothing(uh, f, h)
    gauss_seidel_smoother(u,rhs,gridsize, stepsize);

}



int main(int argc, char* argv[]) {
    //Usage
    if (argc < 3) {
        cout << "Usage: ./mgsolve <lvl> <iterations>" << endl;
        return -1;
    }
    //parse args
    int lvl = atoi(argv[1]);
    int iterations = atoi(argv[2]);

    int number_of_unknowns = (int) pow(pow(2, lvl)-1, 2);
    int grid_size = (int) pow(2,lvl) + 1;
    double step_size = 1./(grid_size-1);
    cout << "Level: " << lvl << "  Iterations: " << iterations << "  N: " << number_of_unknowns << "  h: " << step_size << endl;

    //initialize right-hand side values
    double *rhs = new double[grid_size*grid_size];
    for (int y = 0; y < grid_size; y++) {
        for (int x = 0; x < grid_size; x++) {
            rhs[x + y*grid_size] =2*M_PI*M_PI*cos(M_PI*x*step_size)*cos(M_PI*y*step_size);
        }
    }

    // ===== TESTS =====
    Grid *test = new Grid(lvl, true);
    test->printGrid();

    multigrid(test->data_, rhs, lvl);

    //Grid *small = test->restrict();
    //small->printGrid();

    //Grid *big = small->interpolate();
    //big->printGrid();
    delete test;
    delete []rhs;
    //delete small;
    //delete big;
    // ===== TESTS =====


   

    //double residual = residual(grid->data_, rhs, grid_size, step_size);

    // === Output ===
    // # x y u(x,y)
    //  0 0 0
    // .. .. ..

    //x und y Werte müssen in dem definierten Bereich liegen -> [0,1] und [0,1]
    std::ofstream fileO("solution.txt");
    fileO << "# x y u(x,y)"<< std::endl;
    //int grid_size = 17; //TODO: change to correct grid size (grid->getSize())
    //int step_size = 17; //TODO: change to correct step size (grid->getStepSize())
    for (int col = 0; col < grid_size; col++) {
        for (int row = 0; row < grid_size; row++) {
            fileO << col*step_size << " " << row*step_size << " " << endl; //TODO: << grid->data_[row*stepsize + col] << endl;
        }
    }
    fileO.close();

    return 0;
}