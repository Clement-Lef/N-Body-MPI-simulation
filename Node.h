//
// Created by Clément Lefebvre on 11.05.17.
//

#ifndef NBODYSIMULATION_NODE_H
#define NBODYSIMULATION_NODE_H


#include "Body.h"
#include <vector>

enum quad{qNW,qNE,qSE,qSW};


class Node {
private:
    double xmin;
    double xmax;
    double ymin;
    double ymax;

    double diag_length;
    Body* body;
    Node* NW;
    Node* NE;
    Node* SE;
    Node* SW;

public:
    double mass;
    double centerx;
    double centery;
    Node(double xmin, double xmax, double ymin, double ymax, Body *body);
    enum quad getQuadrant(double x, double y, double xmin, double xmax, double ymin, double ymax);

    void updateMass(Node* node, Body* nbody);

    void insertBody(Body* bodyInsert, Node* node);

    void computeForce(Node* node, Body* nbody, double threshold);

    std::vector<double> searchForBodyQuadrant(Node *node, Body *nbody);

    Body* searchForQuadrant(Node* node, double x, double y);

    void checkCollision(Node* node, Body* nbody);

    void freeTree(Node* node);



    void checkFarSpace(Node* root, Body* nbody,double threshold_farspace);

};


#endif //NBODYSIMULATION_NODE_H
