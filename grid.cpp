#include "grid.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>

Grid::Grid(int lvl) {
    this->level_ = lvl;
    this->size_ = (int) pow(2,lvl) + 1;
    this->step_size_ = 1./(size_-1);
    this->data_ = new double[size_*size_];
    for (int x = 0; x < size_; x++) {
        for (int y = 0; y < size_; y++) {
            int absolute_position = x + (size_*y);
            double value = (x==0 || y==0 || x==size_-1 || y==size_-1) ? cos(M_PI*x*step_size_)*cos(M_PI*y*step_size_) : 0.00;
            this->data_[absolute_position] = value;
        }
    }
}

int Grid::getSize() const {
    return this->size_;
}

Grid::~Grid() {
    delete this->data_;
}

void Grid::printGrid() {
    for (int y = 0; y < size_; y++) {
        for (int x = 0; x < size_; x++) {
            int absolute_position = x + (size_*y);
            std::cout << this->data_[absolute_position] << std::setprecision(2) << " ";
        }
        std::cout << std::endl;
    }
}

Grid *Grid::restrict() {
    Grid *smaller = new Grid(this->level_-1);
    smaller->weightedRestriction(this, this->size_);
    return smaller;
}

void Grid::weightedRestriction(Grid *bigger, int size_big) {
    double *big_d = bigger->data_;
    //go over smaller (inner) grid and restrict from bigger one
    for (int x = 1; x < size_-1; x++) {
        for (int y = 1; y < size_-1; y++) {
            int absolute_position_small = x + (size_*y);
            int positionBig = 2 * x + (size_big * y * 2);
                                        //top left                   top right
            double weighted_sum = 0.125*(big_d[positionBig-size_big-1] + big_d[positionBig-size_big+1]
                                        //bottom left               bottom right
                                  + big_d[positionBig+size_big-1]+big_d[positionBig+size_big+1])
                                        //left                          right
                                  + 0.25*(big_d[positionBig - 1] + big_d[positionBig + 1]
                                        //bottom                        top
                                  + big_d[positionBig + size_big] + big_d[positionBig - size_big])
                                        //self
                                  + 0.5 * (big_d[positionBig]);
            this->data_[absolute_position_small] = weighted_sum;
        }
    }
}
