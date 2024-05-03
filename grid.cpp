#include "grid.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>

Grid::Grid(int lvl, bool initBoundary, bool initInnerPoints) {
    this->level_ = lvl;
    this->size_ = (int) pow(2,lvl) + 1;
    this->step_size_ = 1./(size_-1);
    this->data_ = new double[size_*size_];
    this->initializedBoundary_ = initBoundary;
    this->initializedInnerPoints_ = initInnerPoints;

    //init boundary: if initBoundary is True: fill with boundary values from exercise sheet; else: fill with zeros
    this->initBoundary(initBoundary);

    // init inner points only if initInnerPoints is True
    if (initInnerPoints) {
        for (int y = 1; y < size_-1; y++) {
            for (int x = 1; x < size_-1; x++) {
                this->data_[x + (size_*y)] = 0;
            }
        }
    }
}


void Grid::initBoundary(bool init) {
    if (init) {
        for (int x = 0; x < size_; x++) {
            this->data_[x] = cos(M_PI * x * step_size_); // lower
            this->data_[x * size_] = cos(M_PI * x * step_size_); // left
            this->data_[(x + 1) * size_ - 1] = -1 * cos(M_PI * x * step_size_); // right
            this->data_[(size_ - 1) * size_ + x] = -1 * cos(M_PI * x * step_size_); // upper
        }
    }
    else {
        for (int x = 0; x < size_; x++) {
            this->data_[x] = 0; // lower
            this->data_[x * size_] = 0; // left
            this->data_[(x + 1) * size_ - 1] = 0; // right
            this->data_[(size_ - 1) * size_ + x] = 0; // upper
        }
    }
}

int Grid::getSize() const {
    return this->size_;
}

Grid::~Grid() {
    delete []this->data_;
}

void Grid::printGrid() {
    for (int y = 0; y < size_; y++) {
        for (int x = 0; x < size_; x++) {
            std::cout << this->data_[x + (size_*y)] << std::setprecision(4) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

Grid *Grid::restrict() {
    // create new grid with level l-1 for coarser grid
    Grid *smaller = new Grid(this->level_-1, this->initializedBoundary_, false);
    // fill coarser grid with values
    smaller->weightedRestriction(this);
    return smaller;
}

void Grid::weightedRestriction(Grid *bigger) {
    int size_big = bigger->getSize();
    double *big_d = bigger->data_;
    //go over smaller (inner) grid and restrict from bigger one
    for (int y = 1; y < size_-1; y++) {
        for (int x = 1; x < size_-1; x++) {
            //find correct indices in both grids
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

Grid *Grid::interpolate() {
    //create new grid with level+1 for finer grid
    Grid *bigger = new Grid(this->level_+1, this->initializedBoundary_, false);
    //fill finer grid with values
    bigger->interpolation(this);
    return bigger;
}

void Grid::interpolation(Grid *smaller) {
    int size_small = smaller->getSize();
    double *small_d = smaller->data_;
    // go over smaller grid and copy/interpolate to big one
    bool yOdd = true;
    // Willst du hier x und y vertauschen, vertausche auch die x,yOdd änderung!!!!!
    for (int y = 1; y < size_-1; y++, yOdd = !yOdd) {
        bool xOdd = true;
        for (int x = 1; x < size_-1; x++, xOdd = !xOdd) {
            int positionBig = x + (size_ * y);
            double value;
            if (xOdd) {
                if (yOdd) { //both odd -> value between 4 points of coarser grid -> weighted sum of 4 points
                                    // top left                                     bottom left
                    value = 0.125 * (small_d[(x-1 + (size_small*(y-1)))/2] + small_d[(x+1 + (size_small*(y-1)))/2]
                                    // bottom right                                 top right
                            + small_d[(x-1 + (size_small*(y+1)))/2] + small_d[(x+1 + (size_small*(y+1)))/2]);
                }
                else { // xOdd, yEven -> horizontal between 2 -> weighted sum of two adjacent points
                    value = 0.25 * (small_d[((x-1) + (size_small * y))/2] + small_d[((x+1) + (size_small * y))/2]);
                }
            }
            else if (yOdd) { // xEven, yOdd -> vertical between 2 points -> weighted sum of two adjacent points
                value = 0.25 * (small_d[(x + (size_small * (y-1)))/2] + small_d[(x + (size_small * (y+1)))/2]);
            }
            else { //both even -> on point -> weight with value from self
                value = 0.5*small_d[x/2 + (size_small * y/2)];
            }
            this->data_[positionBig] = value;
        }
    }
}
