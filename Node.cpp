//
// Created by Clément Lefebvre on 11.05.17.
//

#include <cmath>
#include "Node.h"

//Constructor of a node
Node::Node(double xmin, double xmax, double ymin, double ymax, Body *body) {

    mass = body->mass;
    centerx = body->x;
    centery = body->y;
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    // diagonal length of the quadrant
    this->diag_length = sqrt(pow(xmax - xmin,2)+pow(ymax-ymin,2));

    this->body = body;
    NW = NULL;
    NE = NULL;
    SE = NULL;
    SW = NULL;
}





//Update the center of mass of the node
void Node::updateMass(Node *node, Body *nbody) {
    node->centerx = (node->centerx*node->mass + nbody->x*nbody->mass)/(node->mass + nbody->mass);
    node->centery = (node->centery*node->mass + nbody->y*nbody->mass)/(node->mass + nbody->mass);
    node->mass += nbody->mass;
}

// Search for the quadrant which contains the body nbody
std::vector<double> Node::searchForBodyQuadrant(Node *node, Body *nbody) {
    double x = nbody->x;
    double y = nbody->y;
    std::vector<double> Quad(6,0);
    if (node->body != nbody) {
        if (x < node->centerx && y < node->centery && node->SW != NULL) {
            return searchForBodyQuadrant(node->SW, nbody);
        } else if (x < node->centerx && y > node->centery && node->NW != NULL) {
            return searchForBodyQuadrant(node->NW, nbody);
        } else if (x > node->centerx && y > node->centery && node->NE != NULL) {
            return searchForBodyQuadrant(node->NE, nbody);
        } else if (x > node->centerx && y < node->centery && node-> SE != NULL) {
            return searchForBodyQuadrant(node->SE, nbody);
        }
    }
    Quad[0] = node->xmin; Quad[1] = node->xmax; Quad[2] = node->ymin; Quad[3] = node->ymax; Quad[4] = node->centerx;
    Quad[5] = node->centery;
    return Quad;
}

// Search for the body contained in a quadrant
Body* Node::searchForQuadrant(Node* node, double x, double y) {
        if (x < node->centerx && y < node->centery && node->SW != NULL) {
            return searchForQuadrant(node->SW, x,y);
        } else if (x < node->centerx && y > node->centery && node->NW != NULL) {
            return searchForQuadrant(node->NW, x,y);
        } else if (x > node->centerx && y > node->centery && node->NE != NULL) {
            return searchForQuadrant(node->NE, x,y);
        } else if (x > node->centerx && y < node->centery && node-> SE != NULL) {
            return searchForQuadrant(node->SE, x,y);
        }
    return node->body;
}

void Node::checkCollision(Node* node, Body* nbody) {
    double epsilon = 1e-6;

    std::vector<double> inQuad(6,0);
    std::vector<double> borderQuad(6,0);
    borderQuad[0] = node->xmin; borderQuad[1] = node->xmax; borderQuad[2] = node->ymin; borderQuad[3] = node->ymax; borderQuad[4] = node->centerx;
    borderQuad[5] = node->centery;
    double collision_threshold = 10;
    Body* Foundbody;
    //Now search for the 8 quadrants around the body -----------------------------
    //Do only 1 collision
    inQuad = searchForBodyQuadrant(node,nbody);
    //Diagonal North West Quadrant
    if (inQuad[0] > borderQuad[0] && inQuad[3] < borderQuad[3]) {
        Foundbody = searchForQuadrant(node, inQuad[0] - epsilon, inQuad[3] + epsilon);
        if (body != NULL) {
            double dx = Foundbody->x - nbody->x;
            double dy = Foundbody->y - nbody->y;
            double dist = sqrt(pow(dx, 2) + pow(dy, 2));
            if (dist < collision_threshold) {
                nbody->ax = -nbody->ax;
                nbody->ay = -nbody->ay;
                nbody->vx = -nbody->vx;
                nbody->vy = -nbody->vy;
                return;
            }
        }
    }
    // North Quadrant
    if (inQuad[3] < borderQuad[3]) {
        Foundbody = searchForQuadrant(node, inQuad[4], inQuad[3] + epsilon);
        if (body != NULL) {
            double dx = Foundbody->x - nbody->x;
            double dy = Foundbody->y - nbody->y;
            double dist = sqrt(pow(dx, 2) + pow(dy, 2));
            if (dist < collision_threshold) {
                nbody->ax = -nbody->ax;
                nbody->ay = -nbody->ay;
                nbody->vx = -nbody->vx;
                nbody->vy = -nbody->vy;
                return;
            }
        }
    }
    // North-East Quadrant
    if (inQuad[1] < borderQuad[1] && inQuad[3] < borderQuad[3]) {
        Foundbody = searchForQuadrant(node, inQuad[1] + epsilon, inQuad[3] + epsilon);
        if (body != NULL) {
            double dx = Foundbody->x - nbody->x;
            double dy = Foundbody->y - nbody->y;
            double dist = sqrt(pow(dx, 2) + pow(dy, 2));
            if (dist < collision_threshold) {
                nbody->ax = -nbody->ax;
                nbody->ay = -nbody->ay;
                nbody->vx = -nbody->vx;
                nbody->vy = -nbody->vy;
                return;
            }
        }
    }

    // East Quadrant
    if (inQuad[1] < borderQuad[1]) {
        Foundbody = searchForQuadrant(node, inQuad[1] + epsilon, inQuad[5]);
        if (body != NULL) {
            double dx = Foundbody->x - nbody->x;
            double dy = Foundbody->y - nbody->y;
            double dist = sqrt(pow(dx, 2) + pow(dy, 2));
            if (dist < collision_threshold) {
                nbody->ax = -nbody->ax;
                nbody->ay = -nbody->ay;
                nbody->vx = -nbody->vx;
                nbody->vy = -nbody->vy;
                return;
            }
        }
    }

    // South-East Quadrant
    if (inQuad[1] < borderQuad[1] && inQuad[2] > borderQuad[2]) {
        Foundbody = searchForQuadrant(node, inQuad[1] + epsilon, inQuad[2] - epsilon);
        if (body != NULL) {
            double dx = Foundbody->x - nbody->x;
            double dy = Foundbody->y - nbody->y;
            double dist = sqrt(pow(dx, 2) + pow(dy, 2));
            if (dist < collision_threshold) {
                nbody->ax = -nbody->ax;
                nbody->ay = -nbody->ay;
                nbody->vx = -nbody->vx;
                nbody->vy = -nbody->vy;
                return;
            }
        }
    }

    // South Quadrant
    if (inQuad[0] > borderQuad[0]) {
        Foundbody = searchForQuadrant(node, inQuad[4], inQuad[2] - epsilon);
        if (body != NULL) {
            double dx = Foundbody->x - nbody->x;
            double dy = Foundbody->y - nbody->y;
            double dist = sqrt(pow(dx, 2) + pow(dy, 2));
            if (dist < collision_threshold) {
                nbody->ax = -nbody->ax;
                nbody->ay = -nbody->ay;
                nbody->vx = -nbody->vx;
                nbody->vy = -nbody->vy;
                return;
            }
        }
    }

    // South-West Quadrant
    if (inQuad[0] > borderQuad[0] && inQuad[2] > borderQuad[2]) {
        Foundbody = searchForQuadrant(node, inQuad[0] - epsilon, inQuad[2] - epsilon);
        if (body != NULL) {
            double dx = Foundbody->x - nbody->x;
            double dy = Foundbody->y - nbody->y;
            double dist = sqrt(pow(dx, 2) + pow(dy, 2));
            if (dist < collision_threshold) {
                nbody->ax = -nbody->ax;
                nbody->ay = -nbody->ay;
                nbody->vx = -nbody->vx;
                nbody->vy = -nbody->vy;
                return;
            }
        }
    }

    // West Quadrant
    if (inQuad[0] > borderQuad[0]) {
        Foundbody = searchForQuadrant(node, inQuad[0] - epsilon, inQuad[5]);
        if (body != NULL) {
            double dx = Foundbody->x - nbody->x;
            double dy = Foundbody->y - nbody->y;
            double dist = sqrt(pow(dx, 2) + pow(dy, 2));
            if (dist < collision_threshold) {
                nbody->ax = -nbody->ax;
                nbody->ay = -nbody->ay;
                nbody->vx = -nbody->vx;
                nbody->vy = -nbody->vy;
                return;
            }
        }
    }
}


// Insert body into the tree
void Node::insertBody(Body *bodyInsert, Node *node) {
    enum quad oldquad, newquad;
    double xmid = node->xmin + 0.5*std::abs(node->xmax-node->xmin);
    double ymid = node->ymin + 0.5*std::abs(node->ymax-node->ymin);

    //Test if node is empty
    if (node->body != NULL) {
        //If not empty get the quadrant where the body is already inserted
        oldquad = getQuadrant(node->body->x, node->body->y, node->xmin,node->xmax,node->ymin,node->ymax);

        switch (oldquad) {
            case qNW:
                node->NW = new Node(node->xmin, xmid, ymid, node->ymax, node->body);
                break;
            case qNE:
                node->NE = new Node(xmid, node->xmax, ymid, node->ymax, node->body);
                break;
            case qSE:
                node->SE = new Node(xmid, node->xmax, node->ymin, node->ymax, node->body);
                break;
            case qSW:
                node->SW = new Node(node->xmin, xmid, node->ymin, ymid, node->body);
                break;
        }
        node->body = NULL;

    }
    //If the node is empty, get the position of the new quadrant to be created
    newquad = getQuadrant(bodyInsert->x, bodyInsert->y,node->xmin,node->xmax,node->ymin,node->ymax);

    //Update the mass of the node
    updateMass(node,bodyInsert);

    //Create the new quadrant and insert the body
    switch (newquad) {
        case qNW:
            if (node->NW == NULL) {
                node->NW = new Node(node->xmin, xmid, ymid, node->ymax, bodyInsert);
            } else {
                insertBody(bodyInsert, node->NW);
            }
            break;
        case qNE:
            if (node->NE == NULL) {
                node->NE = new Node(xmid, node->xmax, ymid, node->ymax, bodyInsert);
            } else {
                insertBody(bodyInsert, node->NE);
            }
            break;
        case qSE:
            if (node->SE == NULL) {
                node->SE = new Node(xmid, node->xmax, node->ymin, node->ymax, bodyInsert);
            } else {
                insertBody(bodyInsert, node->SE);
            }
            break;
        case qSW:
            if (node->SW == NULL) {
                node->SW = new Node(node->xmin, xmid, node->ymin, ymid, bodyInsert);
            } else {
                insertBody(bodyInsert, node->SW);
            }
            break;
    }




}


// Get the quadrant of a certain position
enum quad Node::getQuadrant(double x, double y, double xmin, double xmax, double ymin, double ymax) {
    double xmid = xmin + 0.5*std::abs(xmax-xmin);
    double ymid = ymin + 0.5*std::abs(ymax-ymin);

    if (y > ymid) {
        if (x > xmid) {
            return qNE;
        }
        else {
            return qNW;
        }
    }
    else {
        if (x > xmid) {
            return qSE;
        }
        else {
            return qSW;
        }
    }
}

// Compute the force of a body in a tree
void Node::computeForce(Node *node, Body *nbody, double threshold) {
    double dx = node->centerx - nbody->x;
    double dy = node->centery - nbody->y;

    double dist = sqrt(pow(dx,2) + pow(dy,2));
    double dist_square = pow(dist,2);

    double G = 6.673e-11;
    double eps = 1e-4; // To not divide by 0
    // Check if the quadrant is too far from the body
    if (((dist/node->diag_length > threshold) || (node->body)) && node->body != nbody) {
        double a = (G*node->mass)/(eps + dist);

        nbody->ax = (nbody->ax +  a*dx/dist_square);
        nbody->ay = (nbody->ay + a*dy/dist_square);
    }
    else {
        if(node->NW) {
            computeForce(node->NW,nbody,threshold);
        }
        if(node->NE) {
            computeForce(node->NE,nbody,threshold);
        }
        if(node->SE) {
            computeForce(node->SE,nbody,threshold);
        }
        if(node->SW) {
            computeForce(node->SW,nbody,threshold);
        }
    }

}


// Delete a tree
void Node::freeTree(Node *node) {

    if (node != NULL) {
        if (node->NW != NULL) {
            freeTree(node->NW);
        }
        if (node->NE != NULL) {
            freeTree(node->NE);
        }
        if (node->SE != NULL) {
            freeTree(node->SE);
        }
        if (node->SW != NULL) {
            freeTree(node->SW);
        }

        delete node;
    }
}



//Check if a body is too far from the center of mass of the system and update a boolean variable
void Node::checkFarSpace(Node *root, Body *nbody, double threshold_farspace) {
    double x_center_mass = ((root->mass-nbody->mass)*root->centerx + nbody->x*nbody->mass)/(root->mass);
    double y_center_mass = ((root->mass-nbody->mass)*root->centery + nbody->y*nbody->mass)/(root->mass);

    double dx = root->centerx - nbody->x;
    double dy = root->centery - nbody->y;


    double dist = sqrt(pow(dx,2) + pow(dy,2));

    if (dist > threshold_farspace) {
        nbody->inSpace = false;
    }
}
