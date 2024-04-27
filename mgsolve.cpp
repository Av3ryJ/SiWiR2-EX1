#include <iostream>
#include <cmath>

#include "grid.hpp"

using namespace std;

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
    int starting_grid_size = (int) pow(2,lvl) + 1;
    cout << "Level: " << lvl << "  Iterations: " << iterations << "  N: " << number_of_unknowns << endl;

    Grid *test = new Grid(lvl);
    cout << "Grid size: " << test->getSize() << endl;
    for (int x = 1; x < 4; x++) {
        for (int y = 1; y < 4; y++) {
            int absolute_position = x + (test->getSize()*y);
            test->data_[absolute_position] = 1;
        }
    }
    test->printGrid();
    Grid *small = test->restrict();
    small->printGrid();

    Grid *big = small->interpolate();
    big->printGrid();
    delete test;
    delete small;
    delete big;

    //TODO: Implement proper V-Cycles (aka implement gaus seidel again..)
    //TODO: write to output file (gnuplot aka steal from last year)

    return 0;
}