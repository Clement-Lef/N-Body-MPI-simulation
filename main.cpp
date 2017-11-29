#include <iostream>
#include <math.h>
#include <fstream>
#include "Body.h"
#include "Node.h"
#include "mpi.h"
#include "Time.h"
#include "stddef.h"





void sendBodies(Body* bodies, int *number_bodies, int rank) {
    /* create a type for struct Body */
    const int nitems=8;
    int          blocklengths[8] = {1,1,1,1,1,1,1,1};
    MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE,MPI_CXX_BOOL};
    MPI_Datatype MPI_Body;
    MPI_Aint     offsets[8];

    offsets[0] = offsetof(Body, x);
    offsets[1] = offsetof(Body, y);
    offsets[2] = offsetof(Body, vx);
    offsets[3] = offsetof(Body, vy);
    offsets[4] = offsetof(Body, ax);
    offsets[5] = offsetof(Body, ay);
    offsets[6] = offsetof(Body, mass);
    offsets[7] = offsetof(Body, inSpace);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_Body);
    MPI_Type_commit(&MPI_Body);




    //Broadcast the number of bodies
    MPI_Bcast(number_bodies, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        //Allocate the memory for the bodies on each core
        bodies = (Body*)malloc((*number_bodies)*sizeof(Body));
    }

    for (int i = 0; i < *number_bodies; i++) {
        //Broadcast the bodies from rank 0 on each cores
        MPI_Bcast(&bodies[i],1,MPI_Body,0,MPI_COMM_WORLD);
    }

    MPI_Type_free(&MPI_Body);

}



// Update the velocity and position of the bodies
void updateBody(Body * bodies, double dt) {
        bodies->vx = (bodies->vx + bodies->ax*dt);
        bodies->vy = (bodies->vy + bodies->ay*dt);
        bodies->x = (bodies->x + bodies->vx*dt);
        bodies->y = (bodies->y + bodies->vy*dt);
}

// Get a random number between 0 and 1
double randomUnit() {
    return (double)rand()/(RAND_MAX);
}

Body* initBody(const int number_bodies,const int number_cluster, const double mass,const double m0, const double rotation_speed = 0.0) {

    Body* bodies;
    // Allocate memory for bodies on processor 0
    bodies = (Body*)malloc(number_bodies*sizeof(Body));

    int bodies_per_cluster = number_bodies/number_cluster;

    for (int j = 0; j < number_cluster; j++) {

        for (int i = 0; i < bodies_per_cluster; i++) {
            int index = j*bodies_per_cluster + i;
            double theta = 2*M_PI*randomUnit(); // Get random angle
            double rand = randomUnit(); // Get random number between 0 and 1

            bodies[index].mass = mass;
            bodies[index].vx = 0.0;
            bodies[index].vy = 0.0;
            bodies[index].inSpace = true;
            bodies[index].x = (2*j + rand*cos(theta));
            bodies[index].y = (2*j + rand*sin(theta));

            bodies[index].vx = (sin(theta)*(rotation_speed+rotation_speed*0.1*randomUnit())/(1+rand));
            bodies[index].vy = (-cos(theta)*(rotation_speed+rotation_speed*0.1*randomUnit())/(1+rand));



        }
        // Set  bodies of superior mass to orbit around
        bodies[j*bodies_per_cluster].mass = m0;
        bodies[j*bodies_per_cluster].x = 2*j;
        bodies[j*bodies_per_cluster].y = 2*j;
        bodies[j*bodies_per_cluster].vx = 0;
        bodies[j*bodies_per_cluster].vy = 0;
    }


    return bodies;
}


void deleteBodies(Body* bodies) {
    free(bodies);
}


void timeStep(Body* bodies, const int number_bodies, const double dt, const double threshold, const int rank,const int size) {

    Node *root;

    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;


    for (int i = 0; i < number_bodies; i++) {
        bodies[i].ax = 0.0;
        bodies[i].ay = 0.0;
        //Search for the size of the quadrant
        xmin = std::min(xmin, bodies[i].x);
        xmax = std::max(xmax, bodies[i].x);
        ymin = std::min(xmin, bodies[i].y);
        ymax = std::max(xmax, bodies[i].y);
    }
    // Create the root node of the tree
    root = new Node(xmin, xmax, ymin, ymax, &bodies[0]);


        // Insert each body into the quadtree
        for (int i = 1; i < number_bodies; i++) {
            if (bodies[i].inSpace) {
                root->insertBody(&bodies[i], root);
            }
        }
    int local_size = 0;

        // Compute the force on each body in parallel
        for (int i = rank; i < number_bodies; i += size) {
            if (bodies[i].inSpace) {
                root->computeForce(root, &bodies[i], threshold);
            }
            local_size++;
        }

        // Send the computed acceleration to each core
        for (int i = 0; i < number_bodies; i++) {
            MPI_Bcast(&(bodies[i].ax), 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
            MPI_Bcast(&(bodies[i].ay), 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
        }



        // Update the body position and velocity on each core
        for (int i = 0; i < number_bodies; i++) {
            updateBody(&bodies[i], dt);
        }



    // Check collision and if the body is in far space in parallel
    for (int i = rank; i < number_bodies; i+=size) {
        root->checkCollision(root, &bodies[i]);

        root->checkFarSpace(root,&bodies[i],5);
    }
    // Broadcast the boolean variable farSpace to each core
    for(int i = 0; i < number_bodies; i++)
    {
        MPI_Bcast(&(bodies[i].inSpace),1,MPI_CXX_BOOL,i%size,MPI_COMM_WORLD);
    }

    root->freeTree(root);



}


void simulate(Body* bodies, const int number_bodies, const double dt, const double tf, const double threshold,const int rank,const int size,const int output_number) {

    std::ofstream output;
    output.open("output.dat");

    int iter = 0;

    for (double t = 0; t <= tf; t = t+dt) {
        timeStep(bodies, number_bodies, dt, threshold, rank,size);

        if (rank == 0  and output_number != 0 and iter%output_number == 0) {
        for (int i = 0; i < number_bodies; i++) {
            output << writePosition(bodies[i]);
        }
        output << std::endl;
        }
        iter++;
    }
}

int main(int argc, char *argv[]) {



    int number_bodies = 200;
    int number_cluster = 1;
    Body* bodies;
    double mass = 2000;
    double dt = 0.01;
    double tf = 10;
    double m0 = 10000000;
    double threshold = 3;
    int output_number = 1;



    if (argc < 6 ) {
        std::cout << "Enter the argument in the following order : " << std::endl;
        std::cout << "Number_bodies Number_clusters dt tf threshold" << std::endl;
        return 1;
    }
    number_bodies = atoi(argv[1]);
    number_cluster = atoi(argv[2]);
    dt = atof(argv[3]);
    tf = atof(argv[4]);
    threshold = atof(argv[5]);
    output_number = atoi(argv[6]);
    int size, rank;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf("nbody: MPI init on %d process\n", size);
    }


    double starttime = 0;
    double ttimefinal = 0;

    if (rank == 0) {
        starttime = second();
    }

    //Initialize the bodies
    bodies = initBody(number_bodies,number_cluster,mass,m0,0.05);

    //Broadcast the bodies from rank 0 to each cores
    sendBodies(bodies,&number_bodies,rank);

    //Run the simulation
    simulate(bodies,number_bodies,dt,tf,threshold,rank,size,output_number);

    deleteBodies(bodies);

    MPI_Finalize();

    if (rank == 0) {
        ttimefinal = second() - starttime;
        printf("%d %f %f %f\n", number_bodies, dt, tf, ttimefinal);

    }


    return 0;
}


