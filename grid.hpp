#include <vector>

#ifndef SIWIR2_EX1_GRID_HPP
#define SIWIR2_EX1_GRID_HPP
class Grid {
public:
    double *data_;

    explicit Grid(int lvl);
    int getSize() const;
    void printGrid();
    Grid *restrict();
    void weightedRestriction(Grid *bigger, int size_big);
    ~Grid();
private:
    int level_;
    int size_;
    double step_size_;
};
#endif //SIWIR2_EX1_GRID_HPP
