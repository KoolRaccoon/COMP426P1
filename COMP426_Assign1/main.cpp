/*
 Issues:
 - The speeds in the left bottom corner galaxy are messed up
 - Need to find a better speed - distance from center ration
 - Find a better way to take care of particles going out of screen?
 -> Currently just pop back from the other side with 0 speed
 - For debug purpose, make two functions for galaxy generation
 -> 1 in center for debug and 2 for demo
 */

//#include "pch.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
//#include <D:\repos\GLEW\include\GL\glew.h>
//#include <D:\repos\GLFW\include\GLFW\glfw3.h>
#include <math.h>
#include <chrono>
#include <thread>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
#include <mutex>
#include "tbb/task_group.h"
#include "tbb/tbb.h"

using namespace tbb;
using namespace std;

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800
double MoveDistance = 0;
GLFWwindow* window;

//Structure used to store the information about planets
struct Point {
    
    double X;
    double Y;
    double Size = 0.002f;
    double Mass = (double(rand()) / double(RAND_MAX) * 1) + 1.0;
    
    vector<double> Force			= { 0.0, 0.0 };
    vector<double> InitialSpeed		= { 0.0, 0.0 };
    vector<double> FinalSpeed		= { 0.0, 0.0 };
    vector<double> DistanceToTravel = { 0.0, 0.0 };
};

//Structure used to create the nodes of the Tree
struct Node {
    vector<Point*> PointsInNodeQuadrant;
    vector<Node*> LeafNodes;
    Node * Parent;
    double QuadrantSize;
    double OriginCoordinates[2] = { 0.0 };
    double CenterOfMass[2]		= { 0.0 };
    int Num_Children			= 0;
    
    int PlanetCount = 0;
    double Mass		= 0.0;
    
    bool MassComputed			= false;
    bool NodeHasOnlyOnePoint	= false;
    bool NodeIsEmpty			= false;
    
};

int TotalPlanets = 500; //The total number of planets to be generated
vector<vector<double>> PlanetCoordinates(TotalPlanets, vector<double>(2));
vector<Point*> Pointss;
double PI = 3.14159265359;

void display(vector<Point*>);
void GenerateRandomPoints1Galaxy(int);
void GenerateRandomPoints2Galaxy(int);
//void display(double, double);
void Tree(Node*);
void ComputeMassDistribution(Node*);
void MassDistribution(Node*);
void CalculateForceOnPoint(Node*);
vector<double> CalculateResultingForce(Node*, Point*);
void CalculateMoveDistance(vector<Point*>&, Node*);
void ResetPointsForce(Point*);
void ResetPointsInitialVelocity(Point*);
void Cleanup(Node*);

const double G = 6.67398 * 0.00000000001;


//class ApplyFoo {
//	float *const my_a;
//public:
//	void operator()(const blocked_range<size_t>& r) const {
//		float *a = my_a;
//		for (size_t i = r.begin(); i != r.end(); ++i)
//			Foo(a[i]);
//	}
//	ApplyFoo(float a[]) :
//		my_a(a)
//	{}
//};

int main() {
    srand(time(NULL));
    //tbb::task_scheduler_init init(300);  // Limiting number of threads to 300
    
//    GenerateRandomPoints1Galaxy(TotalPlanets);
    GenerateRandomPoints2Galaxy(TotalPlanets);
    
    //Initialize the library
    if (!glfwInit()) {
        std::cout << "Error initializing GLFW" << std::endl;
        return -1;
    }
    //     Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Hello World", NULL, NULL);
    if (!window)
    {
        std::cout << "Error creating window" << std::endl;
        glfwTerminate();
        return -1;
    }
    
    //     Make the window's context current
    glfwMakeContextCurrent(window);
    
    //     Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        glClearColor(0.0, 0.0, 0.0, 0.5);
        
        //clear color and depth buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        //Initializing the Root node of the tree
        Node * Root = new Node;
        Root->QuadrantSize = 2.0f;
        Root->OriginCoordinates[0] = 0.0f;
        Root->OriginCoordinates[1] = 0.0f;
        
        //Populating the root node with the previously generated points
        for (int i = 0; i < Pointss.size(); i++) {
            Root->PointsInNodeQuadrant.push_back(Pointss[i]);
            Root->PlanetCount++;
        }
        
        display(Pointss);//Displaying the points
        Tree(Root);//Building the tree
        //ComputeMassDistribution(Root);
        MassDistribution(Root);
        CalculateForceOnPoint(Root);
        CalculateMoveDistance(Pointss, Root);
        
        
        //ResetPointsForce(Pointss);
        Cleanup(Root);//clear up the memory used by the tree
        delete Root;
        //std::cout << "Finished iterating through game loop" << std::endl;
        
        //Swap front and back buffers
        glfwSwapBuffers(window);
        
        //Poll for and process events
        glfwPollEvents();
    }
    glfwTerminate();
    return 0;
}

//Randomly generating points on the screen
void GenerateRandomPoints2Galaxy(int TotalPlanets) {
    
    double R1 = 0.2;
    double R2 = 0.5;
    
    Point * BlackHole1 = new Point;
    BlackHole1->Mass = 20000000000;
    //BlackHole->Mass = (double(rand()) / double(RAND_MAX) * 10000000) + 1.0;
    BlackHole1->X = 0.5;
    BlackHole1->Y = 0.5;
    BlackHole1->InitialSpeed[0] = -0.4;
    BlackHole1->Size = 0.005f;
    //BlackHole1->InitialSpeed[1] = -0.2;
    
    
    Point * BlackHole2 = new Point;
    BlackHole2->Mass = 20000000000;
    //BlackHole->Mass = (double(rand()) / double(RAND_MAX) * 10000000) + 1.0;
    BlackHole2->X = -0.5;
    BlackHole2->Y = -0.5;
    BlackHole2->InitialSpeed[0] = 0.4;
    BlackHole2->Size = 0.005f;
    
    //BlackHole2->InitialSpeed[1] = 0.2;
    
    for (int i = 0; i< TotalPlanets; i++) {
        
        if (i == 0) {// Inserting Blackhole in galaxy 1
            Pointss.push_back(BlackHole1);
        }
        else if (i == TotalPlanets / 2) { // Inserting Black Hole in Galaxy 2
            Pointss.push_back(BlackHole2);
        }
        else {
            
            // Generating coordinates in such a way that the points appear in a circle on the screen
            Point * P = new Point;
            
            //Generating coordinates for galaxy 1
            double a1 = ((double(rand()) / double(RAND_MAX)) * 1.0) * 2 * PI;
            double r1 = R1 * sqrt((double(rand()) / double(RAND_MAX)) * 1.0);
            
            // Generating coordinates for galaxy 2;
            double a2 = ((double(rand()) / double(RAND_MAX)) * 1.0) * 2 * PI;
            double r2 = R2 * sqrt((double(rand()) / double(RAND_MAX)) * 1.0);
            
            //int randGalaxy = int(rand()) % 2;
            double a3 = 0.0;
            if (i<TotalPlanets / 2) {
                //Calculating Initial Velocity of the Planet
                a3 = a1 - (PI / 2);
//                P->InitialSpeed[0] = -0.4 + (1.5 / R1) * cos(a3) * (R1 - r1) / R1;
//                P->InitialSpeed[1] = (1.5 / R1) * sin(a3) * (R1 - r1) / R1;
            
                float orbital_velocity = sqrt((G*P->Mass*BlackHole1->Mass)/r1);
                //cout << "orbital_velocity = " << orbital_velocity << endl;
                P->InitialSpeed[0] =0.6 * orbital_velocity * cos(a3) - 0.4;
                P->InitialSpeed[1] =0.6 * orbital_velocity * sin(a3);
                
                
                PlanetCoordinates[i][0] = BlackHole1->X + r1 * cos(a1);
                PlanetCoordinates[i][1] = BlackHole1->Y + r1 * sin(a1);
            }
            else {
                //Calculating Initial Velocity of the Planet
                a3 = a2 - (PI / 2);
//                P->InitialSpeed[0] = 0.4 + (1.5 / R2) * cos(a3) * (R2 - r2) / R2;
//                P->InitialSpeed[1] = (1.5 / R2) * sin(a3) * (R2 - r2) / R2;

                float orbital_velocity = sqrt((G*P->Mass*BlackHole2->Mass)/r2);
                //cout << "orbital_velocity = " << orbital_velocity << endl;
                P->InitialSpeed[0] = 0.6 * orbital_velocity * cos(a3) + 0.4;
                P->InitialSpeed[1] = 0.6 * orbital_velocity * sin(a3);
                
                PlanetCoordinates[i][0] = BlackHole2->X + r2 * cos(a2);
                PlanetCoordinates[i][1] = BlackHole2->Y + r2 * sin(a2);
                
            }
            bool DuplicatePoint = false; //Used to make sure that there are no duplicate points generated
            
            P->X = PlanetCoordinates[i][0];
            P->Y = PlanetCoordinates[i][1];
            
            for (int j = 1; j < Pointss.size(); j++) {
                if (Pointss[j]->X == P->X && Pointss[j]->Y == P->Y) {
                    i--;
                    DuplicatePoint = true;
                }
            }
            if (DuplicatePoint == false)
                Pointss.push_back(P); //Storing all the points in a vector of pointers
        }
    }
}

//===========================================================================================//
//Randomly generating points on the screen
void GenerateRandomPoints1Galaxy(int TotalPlanets) {
    
    double R1 = 0.7;
    
    Point * BlackHole1 = new Point;
    BlackHole1->Mass = 100000000000;
    //BlackHole->Mass = (double(rand()) / double(RAND_MAX) * 10000000) + 1.0;
    BlackHole1->X = 0.0;
    BlackHole1->Y = 0.0;
    BlackHole1->Size = 0.005f;
    //    BlackHole1->InitialSpeed[1] = -0.2;
    
    for (int i = 0; i< TotalPlanets; i++) {
        
        if (i == 0) {// Inserting Blackhole in galaxy 1
            Pointss.push_back(BlackHole1);
        }
        else {
            
            // Generating coordinates in such a way that the points appear in a circle on the screen
            Point * P = new Point;
            
            //Generating coordinates for galaxy 1
            double a1 = ((double(rand()) / double(RAND_MAX)) * 1.0) * 2 * PI;
            double r1 = R1 * sqrt((double(rand()) / double(RAND_MAX)) * 1.0);
            
            //Calculating Initial Velocity of the Planet
            double a3 = a1 - (PI / 2);
//            if (r1 > (R1 - R1/2)){
//                P->InitialSpeed[0] = (1.5 / R1) * cos(a3) * r1;
//                P->InitialSpeed[1] = (1.5 / R1) * sin(a3) * r1;
//            }
//            else{
//                P->InitialSpeed[0] = (1.5 / R1) * cos(a3) * (R1 - r1) / R1;
//                P->InitialSpeed[1] = (1.5 / R1) * sin(a3) * (R1 - r1) / R1;
//                
//            }
            float orbital_velocity = sqrt((G*P->Mass*BlackHole1->Mass)/r1);
            //cout << "orbital_velocity = " << orbital_velocity << endl;
            P->InitialSpeed[0] =0.6 * orbital_velocity * cos(a3);
            P->InitialSpeed[1] =0.6 * orbital_velocity * sin(a3);
            
            //int randGalaxy = int(rand()) % 2;
            
            PlanetCoordinates[i][0] = BlackHole1->X + r1 * cos(a1);
            PlanetCoordinates[i][1] = BlackHole1->Y + r1 * sin(a1);
            
            bool DuplicatePoint = false; //Used to make sure that there are no duplicate points generated
            
            P->X = PlanetCoordinates[i][0];
            P->Y = PlanetCoordinates[i][1];
            
            for (int j = 1; j < Pointss.size(); j++) {
                if (Pointss[j]->X == P->X && Pointss[j]->Y == P->Y) {
                    i--;
                    DuplicatePoint = true;
                }
            }
            if (DuplicatePoint == false)
                Pointss.push_back(P); //Storing all the points in a vector of pointers
        }
    }
}

//Function used to display the points on the screen
void display(vector<Point*> Points) {
    
    for (int j = 0; j<Points.size(); j++) {
        double P_Size = Points[j]->Size;
        if (Points[j]->X + P_Size < 1.0 || Points[j]->X - P_Size > -1.0 || Points[j]->Y + P_Size < 1.0 || Points[j]->Y - P_Size > -1.0) {
            if (j < Points.size() / 2) {
                //Gives the color of the points that will appear on the screen for the first galaxy
                glColor3f(1, 0.5, 0.5);
            }
            else {
                //Gives the color of the points that will appear on the screen for the second galaxy
                glColor3f(0.5, 0.5, 1);
            }
            glBegin(GL_POLYGON);
            
            glVertex2f((Points[j]->X - P_Size), (Points[j]->Y + P_Size));
            glVertex2f((Points[j]->X - P_Size), (Points[j]->Y - P_Size));
            glVertex2f((Points[j]->X + P_Size), (Points[j]->Y - P_Size));
            glVertex2f((Points[j]->X + P_Size), (Points[j]->Y + P_Size));
            
            glEnd();
            glPopMatrix();
        }
    }
}


void Tree(Node * Current_Node) {
    //cout << "Calling Tree Function" << endl;
    if (Current_Node->PointsInNodeQuadrant.size() <= 1) {
        if (Current_Node->PointsInNodeQuadrant.size() == 0)
            Current_Node->NodeIsEmpty = true;
        else
            Current_Node->NodeHasOnlyOnePoint = true;
    }
    else {
        
        Current_Node->NodeIsEmpty = false;
        Current_Node->NodeHasOnlyOnePoint = false;
        
        //Each leaf will represent a different quadrant. Leaf1 is Quadrant 1, Leaf2 is Quadrant 2, ... and so forth.
        Node *Leaf1 = new Node;
        Node *Leaf2 = new Node;
        Node *Leaf3 = new Node;
        Node *Leaf4 = new Node;
        
        Leaf1->QuadrantSize = Current_Node->QuadrantSize / 2;
        Leaf2->QuadrantSize = Current_Node->QuadrantSize / 2;
        Leaf3->QuadrantSize = Current_Node->QuadrantSize / 2;
        Leaf4->QuadrantSize = Current_Node->QuadrantSize / 2;
        
        //Finding New Origin Coordinates
        Leaf1->OriginCoordinates[0] = Current_Node->OriginCoordinates[0] + Current_Node->QuadrantSize / 2;
        Leaf1->OriginCoordinates[1] = Current_Node->OriginCoordinates[1] + Current_Node->QuadrantSize / 2;
        
        Leaf2->OriginCoordinates[0] = Current_Node->OriginCoordinates[0] - Current_Node->QuadrantSize / 2;
        Leaf2->OriginCoordinates[1] = Current_Node->OriginCoordinates[1] + Current_Node->QuadrantSize / 2;
        
        Leaf3->OriginCoordinates[0] = Current_Node->OriginCoordinates[0] - Current_Node->QuadrantSize / 2;
        Leaf3->OriginCoordinates[1] = Current_Node->OriginCoordinates[1] - Current_Node->QuadrantSize / 2;
        
        Leaf4->OriginCoordinates[0] = Current_Node->OriginCoordinates[0] + Current_Node->QuadrantSize / 2;
        Leaf4->OriginCoordinates[1] = Current_Node->OriginCoordinates[1] - Current_Node->QuadrantSize / 2;
        
        //Find the number of points in the quadrant.
        for (int j = 0; j < Current_Node->PointsInNodeQuadrant.size(); j++) {
            //Populating each node's planet list with planets that belong in it's quadrant.
            
            if (Current_Node->PointsInNodeQuadrant[j]->X >= Current_Node->OriginCoordinates[0] && Current_Node->PointsInNodeQuadrant[j]->Y >= Current_Node->OriginCoordinates[1]) {
                Leaf1->PointsInNodeQuadrant.push_back(Current_Node->PointsInNodeQuadrant[j]);
                Leaf1->PlanetCount++;
            }
            else if (Current_Node->PointsInNodeQuadrant[j]->X < Current_Node->OriginCoordinates[0] && Current_Node->PointsInNodeQuadrant[j]->Y >= Current_Node->OriginCoordinates[1]) {
                Leaf2->PointsInNodeQuadrant.push_back(Current_Node->PointsInNodeQuadrant[j]);
                Leaf2->PlanetCount++;
            }
            else if (Current_Node->PointsInNodeQuadrant[j]->X < Current_Node->OriginCoordinates[0] && Current_Node->PointsInNodeQuadrant[j]->Y < Current_Node->OriginCoordinates[1]) {
                Leaf3->PointsInNodeQuadrant.push_back(Current_Node->PointsInNodeQuadrant[j]);
                Leaf3->PlanetCount++;
            }
            else if (Current_Node->PointsInNodeQuadrant[j]->X >= Current_Node->OriginCoordinates[0] && Current_Node->PointsInNodeQuadrant[j]->Y < Current_Node->OriginCoordinates[1]) {
                Leaf4->PointsInNodeQuadrant.push_back(Current_Node->PointsInNodeQuadrant[j]);
                Leaf4->PlanetCount++;
            }
        }
        
        if (Leaf1->PointsInNodeQuadrant.size() >= 1) {
            Current_Node->Num_Children++;
            Leaf1->Parent = Current_Node;
            Current_Node->LeafNodes.push_back(Leaf1);
            //thread ExpandLeaf1(Tree, Leaf1);
            //ExpandLeaf1.join();
            Tree(Leaf1);
            //(Leaf1);
        }
        if (Leaf2->PointsInNodeQuadrant.size() >= 1) {
            Current_Node->Num_Children++;
            Leaf2->Parent = Current_Node;
            Current_Node->LeafNodes.push_back(Leaf2);
            //thread ExpandLeaf2(Tree, Leaf2);
            //ExpandLeaf2.join();
            //Tree(Leaf2);
            Tree(Leaf2);
        }
        if (Leaf3->PointsInNodeQuadrant.size() >= 1) {
            Current_Node->Num_Children++;
            Leaf3->Parent = Current_Node;
            Current_Node->LeafNodes.push_back(Leaf3);
            //thread ExpandLeaf3(Tree, Leaf3);
            //ExpandLeaf3.join();
            //Tree(Leaf3);
            Tree(Leaf3);
        }
        if (Leaf4->PointsInNodeQuadrant.size() >= 1) {
            Current_Node->Num_Children++;
            Leaf4->Parent = Current_Node;
            Current_Node->LeafNodes.push_back(Leaf4);
            //thread ExpandLeaf4(Tree, Leaf4);
            //ExpandLeaf4.join();
            //Tree(Leaf4);
            Tree(Leaf4);
        }
    }
}

//Computing the Mass of each Node and it's center of mass
void ComputeMassDistribution(Node *Parent) {
    if (Parent->NodeHasOnlyOnePoint == true) {
        Parent->CenterOfMass[0] = Parent->PointsInNodeQuadrant[0]->X;
        Parent->CenterOfMass[1] = Parent->PointsInNodeQuadrant[0]->Y;
        Parent->Mass = Parent->PointsInNodeQuadrant[0]->Mass;
    }
    else {
        for (int i = 0; i < Parent->LeafNodes.size(); i++) {
            ComputeMassDistribution(Parent->LeafNodes[i]);
            Parent->Mass += Parent->LeafNodes[i]->Mass;
            Parent->CenterOfMass[0] += (Parent->LeafNodes[i]->Mass * Parent->LeafNodes[i]->CenterOfMass[0]);
            Parent->CenterOfMass[1] += (Parent->LeafNodes[i]->Mass * Parent->LeafNodes[i]->CenterOfMass[1]);
        }
        Parent->CenterOfMass[0] /= Parent->Mass;
        Parent->CenterOfMass[1] /= Parent->Mass;
    }
}



//Non-Recursive Mass Distribution
void MassDistribution(Node* Root) {
    Node* Current_Node;
    bool TreeTraversed = false;
    Current_Node = Root;
    
    while (Root->MassComputed == false) {
        if (Current_Node->NodeHasOnlyOnePoint == true) {
            Current_Node->CenterOfMass[0] = Current_Node->PointsInNodeQuadrant[0]->X;
            Current_Node->CenterOfMass[1] = Current_Node->PointsInNodeQuadrant[0]->Y;
            Current_Node->Mass = Current_Node->PointsInNodeQuadrant[0]->Mass;
            
            Current_Node->MassComputed = true;
            Current_Node = Current_Node->Parent;
        }
        else {
            int i = 0;
            //bool LeafNodeVisited = false;
            for (; i < Current_Node->LeafNodes.size(); i++) {
                //Visit the leaf nodes first
                if (Current_Node->LeafNodes[i]->MassComputed == false) { // && LeafNodeVisited == false) {
                    Current_Node = Current_Node->LeafNodes[i];
                    //LeafNodeVisited = true;
                    i = 0;
//                    break;
                }
                //Compute mass for Quadrant
                else if ( i == Current_Node->Num_Children - 1){
                    for (int j = 0; j < Current_Node->Num_Children; j++) {
                        Current_Node->Mass += Current_Node->LeafNodes[j]->Mass;
                        Current_Node->CenterOfMass[0] += (Current_Node->LeafNodes[j]->Mass * Current_Node->LeafNodes[j]->CenterOfMass[0]);
                        Current_Node->CenterOfMass[1] += (Current_Node->LeafNodes[j]->Mass * Current_Node->LeafNodes[j]->CenterOfMass[1]);
                    }
                    Current_Node->CenterOfMass[0] /= Current_Node->Mass;
                    Current_Node->CenterOfMass[1] /= Current_Node->Mass;

                    Current_Node->MassComputed = true;
                    i = Current_Node->Num_Children;
                    if (Current_Node == Root)
                        Current_Node = Root;
                    else
                        Current_Node = Current_Node->Parent;
                }
            }
        }
    }
}

//Computes the Total Forces acting on each Planet.
void CalculateForceOnPoint(Node * Root) {
    for ( int i = 0; i < Pointss.size(); i++){
        vector<double> Force = { 0.0, 0.0 };
        Force = CalculateResultingForce(Root, Pointss[i]);
        Pointss[i]->Force[0] += Force[0];
        Pointss[i]->Force[1] += Force[1];
    }
    
    //parallel_for(size_t(0), Pointss.size(), [&](size_t i) {
    //	vector<double> Force = { 0.0, 0.0 };
    //	Force = CalculateResultingForce(Root, Pointss[i]);
    //	Pointss[i]->Force[0] += Force[0];
    //	Pointss[i]->Force[1] += Force[1];
    
    //});
}

// Computes the Total Force of the Tree that is acting upon a Certain Planet.
vector<double> CalculateResultingForce(Node *Parent, Point *TargetPlanet) {
    vector<double> SumOfForces = { 0.0, 0.0 };
    
    if (Parent->PlanetCount == 1) {
        //Compute Force between two planets
        double dis = 0;
        double Force = 0;
        double Theta = 0;
        double Xcomponent = 0;
        double Ycomponent = 0;
        Xcomponent = Parent->PointsInNodeQuadrant[0]->X - TargetPlanet->X;
        Ycomponent = Parent->PointsInNodeQuadrant[0]->Y - TargetPlanet->Y;
        Theta = atan2(Ycomponent, Xcomponent);
        dis = sqrt(pow((Xcomponent), 2.0) + pow((Ycomponent), 2.0));
        
        //If it's the same point or the points are overlapping, the resulting force will be 0.
        if (dis == 0)
            Force = 0;
        else
            Force = (G*TargetPlanet->Mass* Parent->Mass) / (dis*dis);
        
        SumOfForces[0] = Force * cos(Theta);
        SumOfForces[1] = Force * sin(Theta);
        
    }
    else {//Determine if the ratio between the target planet and the node is smaller than the fixed value of 1.
        double r = 0;
        double d = 0;
        double Ratio = 0;
        r = sqrt(pow((Parent->CenterOfMass[0] - TargetPlanet->X), 2.0) + pow((Parent->CenterOfMass[1] - TargetPlanet->Y), 2.0));
        d = Parent->QuadrantSize;
        if (r == 0)
            Ratio = 1;
        else
            Ratio = d / r;
        if (Ratio < 1) {
            //Compute Force between Target Point and Node
            double Force = 0;
            double Theta = 0;
            double Xcomponent = 0;
            double Ycomponent = 0;
            double dis = 0;
            Xcomponent = Parent->CenterOfMass[0] - TargetPlanet->X;
            Ycomponent = Parent->CenterOfMass[1] - TargetPlanet->Y;
            Theta = atan2(Ycomponent, Xcomponent);
            //                Theta = tan((Ycomponent)/(Xcomponent));
            dis = sqrt(pow((Xcomponent), 2.0) + pow((Ycomponent), 2.0));
            
            //If the point sits on top of the Center of Mass of the Node then the resulting force will be 0.
            if (dis == 0)
                Force = 0;
            else
                Force = (G*TargetPlanet->Mass * Parent->Mass) / (dis*dis);
            
            SumOfForces[0] = Force * cos(Theta);
            SumOfForces[1] = Force * sin(Theta);
            Theta = 0;
            dis = 0;
            Force = 0;
        }
        else {
            for (int i = 0; i < Parent->LeafNodes.size(); i++) {
                vector<double> IndividualForces = { 0.0, 0.0 };
                IndividualForces = CalculateResultingForce(Parent->LeafNodes[i], TargetPlanet);
                SumOfForces[0] += IndividualForces[0];
                SumOfForces[1] += IndividualForces[1];
            }
        }
    }
    return SumOfForces;
}


void CalculateMoveDistance(vector<Point*> &Points, Node * Root) {
    
    for (int i = 0; i < Points.size(); i++) {
        //Initialize the distance vectors of the points
        double Distance[2] = { 0.0 };
        //double DistanceLeftToTravel[2]  =  {0.0};
        double Time = 0.002; //time frame is 1ms
							 //double Theta                    =  0.0;
							 //double r                        =  0.0;
        const double MAX_ACC =  20000000.0;
        const double MIN_ACC = -20000000.0;
        const double MAX_SPEED =  8.15;
        const double MIN_SPEED = -8.15;
        
        //Determine the acceleration acting on each point
        double Accelaration[2] = { 0.0 };
        Accelaration[0] = Points[i]->Force[0] / Points[i]->Mass;
        Accelaration[1] = Points[i]->Force[1] / Points[i]->Mass;
        
        // Limiting the Acceleration if it is too high
        if (Accelaration[0] > MAX_ACC)
            Accelaration[0] = MAX_ACC;
        
        else if (Accelaration[0] < MIN_ACC)
            Accelaration[0] = MIN_ACC;
        
        
        //Initializing speed on the Y-axis
        if (Accelaration[1] > MAX_ACC)
            Accelaration[1] = MAX_ACC;
        
        else if (Accelaration[1] < MIN_ACC)
            Accelaration[1] = MIN_ACC;
        
        //Calculating the distance each point will travel on the screen
//        Distance[0] = Points[i]->InitialSpeed[0] * Time + Accelaration[0] * pow(Time, 2) / 2;
//        Distance[1] = Points[i]->InitialSpeed[1] * Time + Accelaration[1] * pow(Time, 2) / 2;
        Distance[0] = Points[i]->InitialSpeed[0] * Time;
        Distance[1] = Points[i]->InitialSpeed[1] * Time;
        
        //Calculating the Final Speed of each point after a time frame
        Points[i]->FinalSpeed[0] = Points[i]->InitialSpeed[0] + Time * Accelaration[0];
        Points[i]->FinalSpeed[1] = Points[i]->InitialSpeed[1] + Time * Accelaration[1];
        
        
        if (Points[i]->InitialSpeed[0] > MAX_SPEED)
            Points[i]->InitialSpeed[0] = MAX_SPEED;
        
        else if (Points[i]->InitialSpeed[0] < MIN_SPEED)
            Points[i]->InitialSpeed[0] = MIN_SPEED;
        else
            Points[i]->InitialSpeed[0] = Points[i]->FinalSpeed[0];
        
        
        //Initializing speed on the Y-axis
        if (Points[i]->InitialSpeed[1] > MAX_SPEED)
            Points[i]->InitialSpeed[1] = MAX_SPEED;
        
        else if (Points[i]->InitialSpeed[1] < MIN_SPEED)
            Points[i]->InitialSpeed[1] = MIN_SPEED;
        else
            Points[i]->InitialSpeed[1] = Points[i]->FinalSpeed[1];
        
        
        //if (i == 5) {
        //	cout << "Points 5  " <<
        //		"Acceleration : X:  " << Accelaration[0] << "  Y:  " << Accelaration[1] << endl;
        //}
        
        //Handling the cases where the point will move off screen and send it coming from the other side instead
        if (Points[i]->X + Distance[0] > 1.0) {//Handling the X-axis
            //DistanceLeftToTravel[0] = Points[i]->X + Distance[0] - 1.0;
            //Points[i]->X = -1.0 + DistanceLeftToTravel[0];
            //Points[i]->X = -1.0 + Points[i]->X;
            //			Points[i]->X = -1.0;
            ResetPointsForce(Points[i]);
            //            ResetPointsInitialVelocity(Points[i]);
        }
        else if (Points[i]->X + Distance[0] < -1.0) {
            //DistanceLeftToTravel[0] = Points[i]->X + Distance[0] + 1.0;
            //Points[i]->X = 1.0 + DistanceLeftToTravel[0];
            //Points[i]->X = 1.0 + Points[i]->X;
            //			Points[i]->X = 1.0;
            ResetPointsForce(Points[i]);
            //            ResetPointsInitialVelocity(Points[i]);
        }
        else {
            Points[i]->X += Distance[0];
            ResetPointsForce(Points[i]);
        }
        
        if (Points[i]->Y + Distance[1] > 1.0) {//Handling the Y-axis
            //DistanceLeftToTravel[1] = Points[i]->Y + Distance[1] - 1.0;
            //Points[i]->Y = -1.0 + DistanceLeftToTravel[1];
            //Points[i]->Y = -1.0 + Points[i]->Y;
            //			Points[i]->Y = -1.0;
            ResetPointsForce(Points[i]);
            //            ResetPointsInitialVelocity(Points[i]);
        }
        else if (Points[i]->Y + Distance[1] < -1.0) {
            //DistanceLeftToTravel[1] = Points[i]->Y + Distance[1] + 1.0;
            //Points[i]->Y = 1.0 + DistanceLeftToTravel[1];
            //Points[i]->Y = 1.0 + Points[i]->Y;
            //			Points[i]->Y = 1.0;
            ResetPointsForce(Points[i]);
            //            ResetPointsInitialVelocity(Points[i]);
        }
        else {
            Points[i]->Y += Distance[1];
            ResetPointsForce(Points[i]);
        }
    }
}

//Reset the Forces acting on each point back to zero
void ResetPointsForce(Point* Point) {
    Point->Force = { 0.0, 0.0 };
}

//Reset the Forces acting on each point back to zero
void ResetPointsInitialVelocity(Point* Point) {
    Point->InitialSpeed = { 0.0, 0.0 };
}

//Cleaning up the Pointers to the Root and all its children nodes
void Cleanup(Node* Parent) {
    if (Parent->LeafNodes.size() > 0) {
        for (int i = 0; i < Parent->LeafNodes.size(); i++) {
            Cleanup(Parent->LeafNodes[i]);
        }
    }
    else {
        delete Parent;
    }
}