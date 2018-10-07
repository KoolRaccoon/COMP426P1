#include "pch.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <math.h>

#include <chrono>
#include <thread>

//#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>

using namespace std;

#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 480
double MoveDistance = 0;
GLFWwindow* window;

//Structure used to store the information about planets
struct Point {
    
    double X;
    double Y;
    double Size = 0.002f;
    double Mass = (double(rand()) / double(RAND_MAX)* 1000000000) + 1.0;
	//double Mass = 1000000000.0;

    vector<double> Force			= {0.0, 0.0};
	vector<double> InitialSpeed		= {0.0, 0.0};
	vector<double> FinalSpeed		= {0.0, 0.0};
	vector<double> DistanceToTravel	= {0.0, 0.0};
};

//Structure used to create the nodes of the Tree
struct Node {
    vector<Point*> PointsInNodeQuadrant;
    vector<Node*> LeafNodes;
	double QuadrantSize;
	double OriginCoordinates[2] = {0.0};

    bool NodeHasOnlyOnePoint    = false;
    bool NodeIsEmpty            = false;

    int PlanetCount       = 0;
	double Mass            = 0.0;
	double CenterOfMass[2] = {0.0};
    //vector<double> Force   = {0.0f, 0.0f};
};

int TotalPlanets = 100; //The total number of planets to be generated
vector<vector<double>> PlanetCoordinates(TotalPlanets, vector<double>(2));
//vector<Point> Points;
vector<Point*> Pointss;

//double Size = 0.005f;

//void drawCircle( double x, double y, double z, double radius, GLint numberOfSides);


void display(vector<Point*>);
void GenerateRandomPoints(int);
//void display(double, double);
void Tree(Node*);
void ComputeMassDistribution(Node*);
void CalculateForceOnPoint(Node*);
vector<double> CalculateResultingForce(Node*, Point*);
void CalculateMoveDistance(vector<Point*>);
void ResetPointsForce(vector<Point*>);
void Cleanup(Node*);

const double G = 6.67398 * 0.00000000001;

int main() {
    srand(time(NULL));

    GenerateRandomPoints(TotalPlanets);
  
    //Initialize the library
    if (!glfwInit()){
        std::cout << "Error initializing GLFW" << std::endl;
        return -1;
    }
//     Create a windowed mode window and its OpenGL context 
    window = glfwCreateWindow(800, 1000, "Hello World", NULL, NULL);
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
			//CalculateForceOnPoint(Root);
			//CalculateMoveDistance(Pointss);
			//ResetPointsForce(Pointss);
			Cleanup(Root);//clear up the memory used by the tree
			delete Root;

			std::cout << "Finished iterating through game loop" << std::endl;

        //Swap front and back buffers
        glfwSwapBuffers(window);
        
        //Poll for and process events
        glfwPollEvents();
    }
    glfwTerminate();
    return 0;
    }

//Randomly generating points on the screen
void GenerateRandomPoints(int TotalPlanets){
	
    for (int i=0; i< TotalPlanets; i++){
		// Generating coordinates in such a way that the points appear in a circle on the screen
		Point * P = new Point;
		//Generating coordinates for galaxy 1
		double R1 = 0.2;
		double a1 = ((double(rand()) / double(RAND_MAX)) * 1.0) * 2 * 3.14159265359;
		double r1 = R1 * sqrt((double(rand()) / double(RAND_MAX)) * 1.0);
		
		// Generating coordinates for galaxy 2;
		double R2 = 0.4;
		double a2 = ((double(rand()) / double(RAND_MAX)) * 1.0) * 2 * 3.14159265359; 
		double r2 = R2 * sqrt((double(rand()) / double(RAND_MAX)) * 1.0);

		//int randGalaxy = int(rand()) % 2;

		if (i<TotalPlanets/2) {

			PlanetCoordinates[i][0] = 0.5 + r1 * cos(a1);
			PlanetCoordinates[i][1] = 0.5 + r1 * sin(a1);
		}
		else {
			PlanetCoordinates[i][0] = -0.5 + r2 * cos(a2);
			PlanetCoordinates[i][1] = -0.5 + r2 * sin(a2);

		}
		bool DuplicatePoint = false; //Used to make sure that there are no duplicate points generated

        //PlanetCoordinates[i][0] = ((double(rand()) / double(RAND_MAX)) * 0.5) + -0.25;
        //PlanetCoordinates[i][1] = ((double(rand()) / double(RAND_MAX)) * 0.5) + -0.25;
        //cout << "displaying X: " << PlanetCoordinates[i][0] << "Displaying Y: " << PlanetCoordinates[i][1] << endl;
        P->X = PlanetCoordinates[i][0];
        P->Y = PlanetCoordinates[i][1];

		for (int j = 0; j < Pointss.size(); j++) {
			if (Pointss[j]->X == P->X && Pointss[j]->Y == P->Y) {
				i--;
				DuplicatePoint = true;
			}
		}
		if (DuplicatePoint == false)
	        Pointss.push_back(P); //Storing all the points in a vector of pointers
    }
}

//Function used to display the points on the screen
void display(vector<Point*> Points) {

    //Square
//        glVertex2f(1.0f, 1.0f);
//        glVertex2f(1.0f, 0.0f);
//        glVertex2f(0.0f, 0.0f);
//        glVertex2f(0.0f, 1.0f);

//    glVertex2f( (P.X - P.Size) - MoveDistance, (P.Y + P.Size));
//    glVertex2f( (P.X - P.Size) - MoveDistance, (P.Y - P.Size));
//    glVertex2f( (P.X + P.Size) - MoveDistance, (P.Y - P.Size));
//    glVertex2f( (P.X + P.Size) - MoveDistance, (P.Y + P.Size));
	
	//Gives the color of the points that will appear on the screen
	glColor3f(1, 0.5, 0.5);
    
	for (int j=0; j<Points.size(); j++){
		if (j < Points.size() / 2) {
			//Gives the color of the points that will appear on the screen
			glColor3f(1, 0.5, 0.5);
		}
		else {
			//Gives the color of the points that will appear on the screen
			glColor3f(0.5, 0.5, 1);
		}
        glBegin(GL_POLYGON);
        
        /*glVertex2f((Points[j]->X - Points[j]->Size) + Pointss[j]->DistanceToTravel[0], (Points[j]->Y + Points[j]->Size) + Pointss[j]->DistanceToTravel[1]);
        glVertex2f((Points[j]->X - Points[j]->Size) + Pointss[j]->DistanceToTravel[0], (Points[j]->Y - Points[j]->Size) + Pointss[j]->DistanceToTravel[1]);
        glVertex2f((Points[j]->X + Points[j]->Size) + Pointss[j]->DistanceToTravel[0], (Points[j]->Y - Points[j]->Size) + Pointss[j]->DistanceToTravel[1]);
        glVertex2f((Points[j]->X + Points[j]->Size) + Pointss[j]->DistanceToTravel[0], (Points[j]->Y + Points[j]->Size) + Pointss[j]->DistanceToTravel[1]);
*/
		glVertex2f((Points[j]->X - Points[j]->Size), (Points[j]->Y + Points[j]->Size));
		glVertex2f((Points[j]->X - Points[j]->Size), (Points[j]->Y - Points[j]->Size));
		glVertex2f((Points[j]->X + Points[j]->Size), (Points[j]->Y - Points[j]->Size));
		glVertex2f((Points[j]->X + Points[j]->Size), (Points[j]->Y + Points[j]->Size));
		//Pointss[j]->DistanceToTravel[0] += 0.01;
		//Pointss[j]->DistanceToTravel[1] += 0.01;
        glEnd();
        glPopMatrix();
    }
    //this_thread::sleep_for(chrono::milliseconds(33));

    //MoveDistance = MoveDistance + 0.01f;
}


void Tree(Node * Parent){
	cout << "Calling Tree Function" << endl;
    if (Parent->PointsInNodeQuadrant.size() <= 1){
        if (Parent->PointsInNodeQuadrant.size() == 0)
            Parent->NodeIsEmpty = true;
        else
            Parent->NodeHasOnlyOnePoint = true;
    }
    else{
        Parent->NodeIsEmpty = false;
        Parent->NodeHasOnlyOnePoint = false;
        
        //Each leaf will represent a different quadrant. Leaf1 is Quadrant 1, Leaf2 is Quadrant 2, ... and so forth.
		Node *Leaf1 = new Node;
		Node *Leaf2 = new Node;
		Node *Leaf3 = new Node;
		Node *Leaf4 = new Node;
        //bool TreeCompleted = false;
        //while (TreeCompleted == false){
            //Calculating the new Quadrant Size of the leaf nodes

        Leaf1->QuadrantSize = Parent->QuadrantSize/2;
        Leaf2->QuadrantSize = Parent->QuadrantSize/2;
        Leaf3->QuadrantSize = Parent->QuadrantSize/2;
        Leaf4->QuadrantSize = Parent->QuadrantSize/2;
        
        //Finding New Origin Coordinates
        Leaf1->OriginCoordinates[0] = Parent->OriginCoordinates[0] + Parent->QuadrantSize/2;
        Leaf1->OriginCoordinates[1] = Parent->OriginCoordinates[1] + Parent->QuadrantSize/2;
        
        Leaf2->OriginCoordinates[0] = Parent->OriginCoordinates[0] - Parent->QuadrantSize/2;
        Leaf2->OriginCoordinates[1] = Parent->OriginCoordinates[1] + Parent->QuadrantSize/2;
        
        Leaf3->OriginCoordinates[0] = Parent->OriginCoordinates[0] - Parent->QuadrantSize/2;
        Leaf3->OriginCoordinates[1] = Parent->OriginCoordinates[1] - Parent->QuadrantSize/2;
        
        Leaf4->OriginCoordinates[0] = Parent->OriginCoordinates[0] + Parent->QuadrantSize/2;
        Leaf4->OriginCoordinates[1] = Parent->OriginCoordinates[1] - Parent->QuadrantSize/2;
        
        //Find the number of points in the quadrant.
        for (int j = 0; j < Parent->PointsInNodeQuadrant.size(); j++){
            //Populating each node's planet list with planets that belong in it's quadrant.
            
            if (Parent->PointsInNodeQuadrant[j]->X >= Parent->OriginCoordinates[0] && Parent->PointsInNodeQuadrant[j]->Y >= Parent->OriginCoordinates[1]){
                Leaf1->PointsInNodeQuadrant.push_back(Parent->PointsInNodeQuadrant[j]);
                Leaf1->PlanetCount++;
            }
            else if (Parent->PointsInNodeQuadrant[j]->X < Parent->OriginCoordinates[0] && Parent->PointsInNodeQuadrant[j]->Y >= Parent->OriginCoordinates[1]){
                Leaf2->PointsInNodeQuadrant.push_back(Parent->PointsInNodeQuadrant[j]);
                Leaf2->PlanetCount++;
            }
            else if (Parent->PointsInNodeQuadrant[j]->X < Parent->OriginCoordinates[0] && Parent->PointsInNodeQuadrant[j]->Y < Parent->OriginCoordinates[1]){
                Leaf3->PointsInNodeQuadrant.push_back(Parent->PointsInNodeQuadrant[j]);
                Leaf3->PlanetCount++;
            }
            else if (Parent->PointsInNodeQuadrant[j]->X >= Parent->OriginCoordinates[0] && Parent->PointsInNodeQuadrant[j]->Y < Parent->OriginCoordinates[1]){
                Leaf4->PointsInNodeQuadrant.push_back(Parent->PointsInNodeQuadrant[j]);
                Leaf4->PlanetCount++;
            }
        }
        
        if (Leaf1->PointsInNodeQuadrant.size() >= 1){
            Parent->LeafNodes.push_back(Leaf1);
			thread ExpandLeaf1(Tree, Leaf1);
            //Tree(Leaf1);
        }
        if (Leaf2->PointsInNodeQuadrant.size() >= 1){
            Parent->LeafNodes.push_back(Leaf2);
			thread ExpandLeaf2(Tree, Leaf2);
            //Tree(Leaf2);
        }
        if (Leaf3->PointsInNodeQuadrant.size() >= 1){
            Parent->LeafNodes.push_back(Leaf3);
			thread ExpandLeaf3(Tree, Leaf3);
            //Tree(Leaf3);
        }
        if (Leaf4->PointsInNodeQuadrant.size() >= 1){
            Parent->LeafNodes.push_back(Leaf4);
			thread ExpandLeaf4(Tree, Leaf4);
            //Tree(Leaf4);
        }
    }
}

//Computing the Mass of each Node and it's center of mass
void ComputeMassDistribution(Node *Parent){
    if (Parent->NodeHasOnlyOnePoint == true){
        Parent->CenterOfMass[0] = Parent->PointsInNodeQuadrant[0]->X;
        Parent->CenterOfMass[1] = Parent->PointsInNodeQuadrant[0]->Y;
        Parent->Mass = Parent->PointsInNodeQuadrant[0]->Mass;
    }
    else {
        for (int i = 0; i < Parent->LeafNodes.size(); i ++){
            ComputeMassDistribution(Parent->LeafNodes[i]);
            Parent->Mass += Parent->LeafNodes[i]->Mass;
            Parent->CenterOfMass[0] += (Parent->LeafNodes[i]->Mass * Parent->LeafNodes[i]->CenterOfMass[0]);
            Parent->CenterOfMass[1] += (Parent->LeafNodes[i]->Mass * Parent->LeafNodes[i]->CenterOfMass[1]);        
		}
		Parent->CenterOfMass[0] /= Parent->Mass;
		Parent->CenterOfMass[1] /= Parent->Mass;
	}
}

//Computes the Total Forces acting on each Planet.
void CalculateForceOnPoint(Node * Root){
    for ( int i = 0; i < Pointss.size(); i++){
        Pointss[i]->Force = CalculateResultingForce(Root, Pointss[i]);
    }
}

// Computes the Total Force of the Tree that is acting upon a Certain Planet.
vector<double> CalculateResultingForce(Node *Parent, Point *TargetPlanet){
    //TargetPlanet->Force[0] = 0;
    //TargetPlanet->Force[1] = 0;
    vector<double> SumOfForces = {0.0, 0.0};
	double ScaleDownFactor = 0.0;
    
    //for (int i = 0; i < Parent->LeafNodes.size(); i ++){
        //if ((Parent->PlanetCount == 1) && !((TargetPlanet->X == Parent->PointsInNodeQuadrant[0]->X) &&(TargetPlanet->Y == Parent->PointsInNodeQuadrant[0]->Y))){
	if (Parent->PlanetCount == 1){
			//Compute Force between two planets
			double dis = 0;
			double Force = 0;
			double Theta = 0;
			double Xcomponent = 0;
			double Ycomponent = 0;
			Xcomponent = Parent->PointsInNodeQuadrant[0]->X - TargetPlanet->X;
			Ycomponent = Parent->PointsInNodeQuadrant[0]->Y - TargetPlanet->Y;
            Theta = atan2(Ycomponent,Xcomponent);
            dis   = sqrt(pow((Xcomponent), 2.0) + pow((Ycomponent),2.0));
			
			//If it's the same point or the points are overlapping, the resulting force will be 0.
			if (dis == 0)
				Force = 0;
			else
	            Force = (G*TargetPlanet->Mass* Parent->Mass)/(dis*dis);
            
			SumOfForces[0] = Force * cos(Theta);
            SumOfForces[1] = Force * sin(Theta);
			//TargetPlanet->Force[0] += Force * cos(Theta);
			//TargetPlanet->Force[1] += Force * sin(Theta);
            Theta = 0;
            dis   = 0;
            Force = 0;
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
		    if (Ratio < 1){
                //Compute Force between Target Point and Node
				double Force = 0;
				double Theta = 0;
				double Xcomponent = 0;
				double Ycomponent = 0;
				double dis = 0;
				Xcomponent = Parent->CenterOfMass[0] - TargetPlanet->X;
				Ycomponent = Parent->CenterOfMass[1] - TargetPlanet->Y;
                Theta = tan((Ycomponent)/(Xcomponent));
                dis   = sqrt(pow((Xcomponent), 2.0) + pow((Ycomponent), 2.0));

				//If the point sits on top of the Center of Mass of the Node then the resulting force will be 0.	
				if (dis == 0)
					Force = 0;
				else
	                Force = (G*TargetPlanet->Mass * Parent->Mass)/(dis*dis);
				
				SumOfForces[0] = Force * cos(Theta);
				SumOfForces[1] = Force * sin(Theta);
                Theta = 0;
                dis   = 0;
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
		//TargetPlanet->Force[0] += SumOfForces[0];
		//TargetPlanet->Force[1] += SumOfForces[1];
    //}
    return SumOfForces;
}


void CalculateMoveDistance(vector<Point*> Points) {
	
	for (int i = 0; i < Points.size(); i++) {
		//Initialize the distance vectors of the points
		double Distance[2] = {0.0};
		double DistanceLeftToTravel[2] = {0.0};
		double TimeValue = 0.0030303; //30 frames per second

		//Determine the acceleration acting on each point
		double Accelaration[2] = {0.0};
		Accelaration[0] = Points[i]->Force[0] / Points[i]->Mass;
		Accelaration[1] = Points[i]->Force[1] / Points[i]->Mass;

		//Finding the final speed of each point after moving for certain time frame
		//Points[i]->FinalSpeed[0] = Accelaration[0] * (1.5) + Points[i]->InitialSpeed[0];
		//Points[i]->FinalSpeed[1] = Accelaration[1] * (1.5) + Points[i]->InitialSpeed[1];

		//Points[i]->DistanceToTravel[0] += Points[i]->FinalSpeed[0] - Accelaration[0] * pow((1 / 30), 2.0);
		//Points[i]->DistanceToTravel[1] += Points[i]->FinalSpeed[1] - Accelaration[1] * pow((1 / 30), 2.0);

		//Calculating the distance each point will travel on the screen
		Distance[0] = Accelaration[0] * pow(TimeValue, 2.0);
		Distance[1] = Accelaration[1] * pow(TimeValue, 2.0);

		//Reinitializing the initial speed of each point after moving
		Points[i]->InitialSpeed = Points[i]->FinalSpeed;
		
		//Handling the cases where the point will move off screen and send it coming from the other side instead
		if (Points[i]->X + Distance[0] > 1.0) {//Handling the X-axis
			DistanceLeftToTravel[0] = Points[i]->X + Distance[0] - 1.0;
			//Points[i]->X = -1.0 + DistanceLeftToTravel[0];
			Points[i]->X = -1.0 + Points[i]->X;
		}
		else if (Points[i]->X + Distance[0] < -1.0){
			DistanceLeftToTravel[0] = Points[i]->X + Distance[0] + 1.0;
			//Points[i]->X = 1.0 + DistanceLeftToTravel[0];
			Points[i]->X = 1.0 + Points[i]->X;
		}
		else {
			Points[i]->X += Distance[0];
		}

		if (Points[i]->Y + Distance[1] > 1.0) {//Handling the Y-axis
			DistanceLeftToTravel[1] = Points[i]->Y + Distance[1] - 1.0;
			//Points[i]->Y = -1.0 + DistanceLeftToTravel[1];
			Points[i]->Y = -1.0 + Points[i]->Y;
		}
		else if (Points[i]->Y + Distance[1] < -1.0) {
			DistanceLeftToTravel[1] = Points[i]->Y + Distance[1] + 1.0;
			//Points[i]->Y = 1.0 + DistanceLeftToTravel[1];
			Points[i]->Y = 1.0 + Points[i]->Y;
		}
		else {
			Points[i]->Y += Distance[1];
		}
	}
}

//Reset the Forces acting on each point back to zero
void ResetPointsForce(vector<Point*> Points) {
	for (int i = 0; i < Points.size(); i++) {
		Points[i]->Force = { 0.0, 0.0 };
	}
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
