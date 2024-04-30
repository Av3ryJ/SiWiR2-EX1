#include <vector>

#ifndef SIWIR2_EX1_GRID_HPP
#define SIWIR2_EX1_GRID_HPP
class Grid {
public:
    double *data_;

    explicit Grid(int lvl, bool initBoundary, bool initInnerPoints);
    void initBoundary(bool init);

    int getSize() const;
    void printGrid();

    Grid *restrict();
    void weightedRestriction(Grid *bigger);

    Grid *interpolate();
    void interpolation(Grid *smaller);

    //TODO: Function returning convergence factor q

    ~Grid();
private:
    int level_;
    int size_;
    double step_size_;
    bool initializedBoundary_;
    bool initializedInnerPoints_;
};
#endif //SIWIR2_EX1_GRID_HPP
