//
// Created by Clément Lefebvre on 29.04.17.
//

#include <cmath>
#include <iostream>
#include <sstream>
#include "Body.h"









std::string writePosition(const Body toWrite) {
    std::ostringstream oss;
    oss << toWrite.x << " " << toWrite.y << " ";
    std::string var = oss.str();
    return var;
}









bool operator==(const Body &lhs, const Body &rhs){
    return lhs.x == rhs.x &&
           lhs.y == rhs.y &&
           lhs.vx == rhs.vx &&
           lhs.vy == rhs.vy &&
           lhs.ax == rhs.ax &&
           lhs.ay == rhs.ay &&
           lhs.mass == rhs.mass;
}

bool operator!=(const Body &lhs, const Body &rhs){
    return !(lhs == rhs);
}


