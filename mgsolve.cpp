#include <iostream>
#include <cmath>

int main(int argc, char* argv[]) {
    //Usage
    if (argc < 3) {
        std::cout << "Usage: ./mgsolve <lvl> <iterations>" << std::endl;
        return -1;
    }
    //parse args
    int lvl = atoi(argv[1]);
    int iterations = atoi(argv[2]);

    int number_of_unknowns = (int) pow(pow(2, lvl)-1, 2);
    std::cout << "Level: " << lvl << "  Iterations: " << iterations << "  N: " << number_of_unknowns << std::endl;
}