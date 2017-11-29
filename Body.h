//
// Created by Clément Lefebvre on 29.04.17.
//

#ifndef NBODYSIMULATION_BODY_H
#define NBODYSIMULATION_BODY_H


#include <iostream>

struct Body {
    double x,y; //positions
    double vx,vy; //velocities
    double ax, ay; // acceleration
    double mass;
    bool inSpace;
};


std::string writePosition(const Body toWrite);


bool operator==(const Body &lhs, const Body &rhs);

bool operator!=(const Body &lhs, const Body &rhs);

#endif //NBODYSIMULATION_BODY_H
