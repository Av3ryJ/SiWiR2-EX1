#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "grid.hpp"
#include "Timer.h"

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
    return sqrt(sum);
}
void calculate_residuum(double *res, double *data, double *rhs, int gridsize, double stepsize) {
    //weights
    double alpha = 4.0/(stepsize*stepsize);
    double beta = 1.0/(stepsize*stepsize);
    for (int row = 1; row < gridsize-1; row++) {
        for (int col = 1; col < gridsize-1; col++) {
            // r = f - A*u
            res[row*gridsize+col] = rhs[row*gridsize + col] - alpha*data[row*gridsize + col]
                              +beta*data[(col-1) + row*gridsize]
                              +beta*data[(col+1) + row*gridsize]
                              +beta*data[col + (row-1)*gridsize]
                              +beta*data[col + (row+1)*gridsize];
        }
    }
} 

void gauss_seidel_smoother(double *data, double *rhs, int gridsize, double stepsize) {
    //weights
    double alpha = 4.0/(stepsize*stepsize);
    double beta = 1.0/(stepsize*stepsize);
    for (int row = 1; row < gridsize-1; ++row) { // row is y col is x
        for (int col = 1; col < gridsize-1; ++col) {    // nicht ueber Rand iterieren
                data[col + row*gridsize] = (1./alpha)*
                        (rhs[col+row*gridsize]
                        +beta*data[(col-1) + row*gridsize]
                        +beta*data[(col+1) + row*gridsize]
                        +beta*data[col + (row-1)*gridsize]
                        +beta*data[col + (row+1)*gridsize]);
        }
    }
}

void addInterpolation(double *data, double *interpol, int gridsize) {
     for (int row = 1; row < gridsize-1; ++row) { 
        for (int col = 1; col < gridsize-1; ++col) {
            // add two data grids, save in first function argument
            data[col+row*gridsize] += interpol[col+row*gridsize]; 
        }
    }    
}

void multigrid(double *u, double *rhs, int lvl) {

    //compute grid and stepsize from level lvl
    int gridsize = (int) pow(2,lvl) + 1;
    double stepsize = 1./(gridsize-1);

    // Presmoothing of Au = b -->uh = Gauss-Seidel(uh, f, h),
    gauss_seidel_smoother(u, rhs, gridsize, stepsize);
    gauss_seidel_smoother(u, rhs, gridsize, stepsize);

    //Compute Residual r = b-Auh 
    Grid *res = new Grid(lvl, false, false);
    calculate_residuum(res->data_, u, rhs, gridsize, stepsize);

    //restrict residual to level l-1  
    Grid *restricted = res->restrict();

    //empty grid eps--> same size as restricted
    Grid *eps = new Grid(lvl-1, false, true); //init with all zeros

    if (lvl == 2) {
        // base case: coarsest grid -> exact solution on 3x3 grid with 1 inner gridpoint
        eps->data_[4] = stepsize*stepsize*restricted->data_[4]/4.;
    }
    else {
        //recursion: call multigrid with eps, restricted and coarser level
        multigrid(eps->data_, restricted->data_, lvl-1);
    }

    // interpolate eps to finer level
    Grid *eps2 = eps->interpolate();
    // add interpolated eps to current approximation for correction
    addInterpolation(u, eps2->data_, gridsize);

    // postsmoothing(uh, f, h)
    gauss_seidel_smoother(u,rhs,gridsize, stepsize);
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

    //right-hand side vector
    double *rhs = new double[grid_size*grid_size];
    // Grid with solution
    Grid *solution = new Grid(lvl, false, false);
    //initialize both rhs and solution
    for (int y = 0; y < grid_size; y++) {
        for (int x = 0; x < grid_size; x++) {
            rhs[x + y*grid_size] =2*M_PI*M_PI*cos(M_PI*x*step_size)*cos(M_PI*y*step_size);
            solution->data_[x + y * grid_size] = cos(M_PI * x * step_size) * cos(M_PI * y * step_size);
        }
    }

    // ===== start multigrid algorithm =====

    Grid *grid = new Grid(lvl, true, true);

    double residual_of_last_iteration = l2_norm(grid->data_, rhs, grid_size, step_size);

    //Timing start
    double time = 100.0;
    siwir::Timer timer;
    //perform given number of iterations
    for (int i = 0; i < iterations; i++) {
        multigrid(grid->data_, rhs, lvl);
        double residual = l2_norm(grid->data_, rhs, grid_size, step_size);
        double discrete_residual = residual / sqrt(number_of_unknowns);
        double convergence_factor = residual / residual_of_last_iteration;
        residual_of_last_iteration = residual;
        cout << "iteration " << i+1 << ": discrete residual norm = " << discrete_residual << "; convergence factor = " << convergence_factor << endl;
    }
    //Timing stoppen & ausgeben
    time = std::min(time, timer.elapsed());

    // compute total error
    double disc_l2_norm_of_error = 0;
    for (int y = 1; y < grid_size; y++) {
        for (int x = 1; x < grid_size; x++) {
            double err = grid->data_[x + grid_size*y] - solution->data_[x + grid_size*y];
            disc_l2_norm_of_error += err * err;
        }
    }
    disc_l2_norm_of_error = sqrt(disc_l2_norm_of_error/(number_of_unknowns));

    cout << endl << "runtime: " << time << "; total error: " << disc_l2_norm_of_error << endl;

    // === Output ===
    // # x y u(x,y)
    //  0 0 0
    // .. .. ..

    //x und y Werte müssen in dem definierten Bereich liegen -> [0,1] und [0,1]
    std::ofstream fileO("solution.txt");
    fileO << "# x y u(x,y)"<< std::endl;
    for (int col = 0; col < grid_size; col++) {
        for (int row = 0; row < grid_size; row++) {
            fileO << col*step_size << " " << row*step_size << " " << grid->data_[row*grid_size + col] << endl;
        }
    }
    fileO.close();

    delete grid;
    delete []rhs;
    delete solution;

    return 0;
}