#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "grid.hpp"

using namespace std;

double residual(double *data, double *rhs, int gridsize, double stepsize) {
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
    double factor = 1.0/(inner_points*inner_points); //TODO: Wurzel aus 1/n ? 
    return sqrt(factor*sum);
}

void gauss_seidel_smoother(double *data,double *rhs, int gridsize,double stepsize){ 
    double alpha = 4.0/(stepsize*stepsize);
    double beta = 1.0/(stepsize*stepsize);
    for (int row = 1; row < ny; ++row) { // row is y col is x
        for (int col = 1; col < nx; ++col) {    // nicht ueber Rand iterieren
                data[col + row*gridsize] = (1/alpha)*
                        (rhs[col+row*gridsize]
                        +beta*data[(col-1) + row*gridsize]
                        +beta*data[(col+1) + row*gridsize]
                        +beta*data[col + (row-1)*gridsize]
                        +beta*data[col + (row+1)*gridsize]);
        }
    }
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
            //rechte Seite sollte f(x,y) sein und nicht g(x,y) oder? 
            rhs[x + y*grid_size] =2*M_PI*M_PI*cos(M_PI*x*step_size)*cos(M_PI*y*step_size);
        }
    }

    // ===== TESTS =====
    Grid *test = new Grid(lvl);
    test->printGrid();

    //Grid *small = test->restrict();
    //small->printGrid();

    //Grid *big = small->interpolate();
    //big->printGrid();
    delete test;
    delete []rhs;
    //delete small;
    //delete big;
    // ===== TESTS =====

    //TODO: Implement proper V-Cycles (aka implement gaus seidel again..)
    // Orientiert an Tafelübung und Wikipedia Pseudo-code
    // uh = V-Cycle(uh, f, h)


    // Presmoothing of Au = b -->uh = Gauss-Seidel(uh, f, h)
    gauss_seidel_smoother(test->data_, rhs, grid_size, step_size)

    //Compute Residual r = b-Auh 

    //restrict residual to level l-1  

    //leeres Gitter anlegen eps--> so groß wie restricted size //set u=0 on l-1 --> auch Randbedingungen werden 0

    // Abbruchbedingung Rekursion: wenn kleinste Gittergröße erreicht ist

    //sonst: eps = V-Cycle(leeres Gitter eps, restricted Gitter, 2*h) --> Rekursion

    // Ab hier wieder aufwärts in der Rekurison 
    // uh = uh + Interpolation(eps)
    // correct: uh = smoothing(uh, f, h)



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