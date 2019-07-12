
#include "cell.h"

cell::cell() {

    Neighbour = new int **[3];
    for (int j = 0; j < 3; j++) {
        Neighbour[j] = new int *[3];
        for (int k = 0; k < 3; k++) {
            Neighbour[j][k] = new int[3];
        }
    }
}
cell::~cell() {}
