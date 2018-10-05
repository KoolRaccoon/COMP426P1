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
GLfloat MoveDistance = 0;
GLFWwindow* window;

struct Point {
    
    float X;
    float Y;
    float Size = 0.002f;
    float Mass = (float(rand()) / float(RAND_MAX)* 1000000) + 1.0;
    vector<float> Force				= {0.0, 0.0};
	vector<float> InitialSpeed		= {0.0, 0.0};
	vector<float> FinalSpeed		= {0.0, 0.0};
	vector<float> DistanceToTravel	= {0.0, 0.0};

};


struct Node {
    vector<Point*> PointsInNodeQuadrant;
    vector<Node*> LeafNodes;
	float QuadrantSize;
	float OriginCoordinates[2] = {0.0};

    bool NodeHasOnlyOnePoint    = false;
    bool NodeIsEmpty            = false;

    int PlanetCount       = 0;
	float Mass            = 0.0;
	float CenterOfMass[2] = {0.0};
    //vector<GLfloat> Force   = {0.0f, 0.0f};
};

int TotalPlanets = 20;
vector<vector<float>> PlanetCoordinates(TotalPlanets, vector<float>(2));
vector<Point> Points;
vector<Point*> Pointss;

GLfloat Size = 0.005f;

//void drawCircle( GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides);


void display(vector<Point*>);
void GenerateRandomPoints(int);
//void display(float, float);
void Tree(Node*);
void ComputeMassDistribution(Node*);
void CalculateForceOnPoint(Node*);
vector<float> CalculateResultingForce(Node*, Point*);
void CalculateMoveDistance(vector<Point*>);
void ResetPointsForce(vector<Point*>);
void Cleanup(Node*);

const float G = 6.67398 * 0.00000000001;

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
        glClearColor(0.2, 0.3, 0.4, 0.5);
        
        //clear color and depth buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        
        glColor3f(1, 0.5, 0.5);
        

		Node * Root = new Node;
		Root->QuadrantSize = 2.0f;
		Root->OriginCoordinates[0] = 0.0f;
		Root->OriginCoordinates[1] = 0.0f;

		for (int i = 0; i < Pointss.size(); i++) {
			Root->PointsInNodeQuadrant.push_back(Pointss[i]);
			Root->PlanetCount++;
		}


            display(Pointss);
			Tree(Root);
			ComputeMassDistribution(Root);
			CalculateForceOnPoint(Root);
			CalculateMoveDistance(Pointss);
			ResetPointsForce(Pointss);
			Cleanup(Root);
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


void GenerateRandomPoints(int TotalPlanets){
	

    for (int i=0; i< TotalPlanets; i++){
		bool DuplicatePoint = false;
        Point * P = new Point;
        PlanetCoordinates[i][0] = ((float(rand()) / float(RAND_MAX)) * 0.5) + -0.25;
        PlanetCoordinates[i][1] = ((float(rand()) / float(RAND_MAX)) * 0.5) + -0.25;
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
	        Pointss.push_back(P);
        
    }
    
 
}

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

    for (int j=0; j<Points.size(); j++){
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
            Tree(Leaf1);
        }
        if (Leaf2->PointsInNodeQuadrant.size() >= 1){
            Parent->LeafNodes.push_back(Leaf2);
            Tree(Leaf2);
        }
        if (Leaf3->PointsInNodeQuadrant.size() >= 1){
            Parent->LeafNodes.push_back(Leaf3);
            Tree(Leaf3);
        }
        if (Leaf4->PointsInNodeQuadrant.size() >= 1){
            Parent->LeafNodes.push_back(Leaf4);
            Tree(Leaf4);
        }
        //}
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
vector<float> CalculateResultingForce(Node *Parent, Point *TargetPlanet){
    //TargetPlanet->Force[0] = 0;
    //TargetPlanet->Force[1] = 0;
    vector<float> SumOfForces = {0, 0};
    float dis = 0;
	float Force = 0;
	float Theta = 0;
    
    for (int i = 0; i < Parent->LeafNodes.size(); i ++){
        if (Parent->LeafNodes[i]->PlanetCount == 1){
            //Compute Force between two planets
			float Xcomponent = 0;
			float Ycomponent = 0;
			Xcomponent = TargetPlanet->X - Parent->LeafNodes[i]->PointsInNodeQuadrant[0]->X;
			Ycomponent = TargetPlanet->Y - Parent->LeafNodes[i]->PointsInNodeQuadrant[0]->Y;
            Theta = atan2(Ycomponent,Xcomponent);
            dis   = sqrt(pow(2.0, (Xcomponent)) + pow(2.0, (Ycomponent)));
            Force = (G*TargetPlanet->Mass* Parent->LeafNodes[i]->Mass)/(dis*dis);
            SumOfForces[0] += Force * cos(Theta);
            SumOfForces[1] += Force * sin(Theta);
            Theta = 0;
            dis   = 0;
            Force = 0;
        }
        else {
			float r = 0;
			float d = 0;
			float Ratio = 0;
			r = sqrt(pow(2.0, (TargetPlanet->X - Parent->LeafNodes[i]->CenterOfMass[0])) + pow(2.0, (TargetPlanet->Y - Parent->LeafNodes[i]->CenterOfMass[1])));
            d = Parent->LeafNodes[i]->QuadrantSize;
			if (r == 0)
				Ratio = 1;
		    if (d/r < 1){
                //Compute Force between Target Planet and Node
				float Force = 0;
				float Theta = 0;
				float Xcomponent = 0;
				float Ycomponent = 0;
				Xcomponent = TargetPlanet->X - Parent->LeafNodes[i]->CenterOfMass[0];
				Ycomponent = TargetPlanet->Y - Parent->LeafNodes[i]->CenterOfMass[1];
                Theta = tan((Ycomponent)/(Xcomponent));
                dis   = sqrt(pow(2.0, (Xcomponent)) + pow(2.0, (Ycomponent)));
                Force = (G*TargetPlanet->Mass* Parent->LeafNodes[i]->Mass)/(dis*dis);
                SumOfForces[0] += Force * cos(Theta);
                SumOfForces[1] += Force * sin(Theta);
                Theta = 0;
                dis   = 0;
                Force = 0;

            }
            else {
                SumOfForces = CalculateResultingForce(Parent->LeafNodes[i], TargetPlanet);
            }
        }
		TargetPlanet->Force[0] += SumOfForces[0];
		TargetPlanet->Force[1] += SumOfForces[1];
    }
    return TargetPlanet->Force;
}


void CalculateMoveDistance(vector<Point*> Points) {
	
	for (int i = 0; i < Points.size(); i++) {
		float Accelaration[2] = { 0.0 };
		Accelaration[0] = Points[i]->Force[0] / Points[i]->Mass;
		Accelaration[1] = Points[i]->Force[1] / Points[i]->Mass;

		Points[i]->FinalSpeed[0] = Accelaration[0] * (1 / 30) + Points[i]->InitialSpeed[0];
		Points[i]->FinalSpeed[1] = Accelaration[1] * (1 / 30) + Points[i]->InitialSpeed[1];

		//Points[i]->DistanceToTravel[0] += Points[i]->FinalSpeed[0] - Accelaration[0] * pow(2.0, (1 / 30));
		//Points[i]->DistanceToTravel[1] += Points[i]->FinalSpeed[1] - Accelaration[1] * pow(2.0, (1 / 30));
		Points[i]->X += Points[i]->FinalSpeed[0] - Accelaration[0] * pow(2.0, (1 / 30));
		Points[i]->Y += Points[i]->FinalSpeed[1] - Accelaration[1] * pow(2.0, (1 / 30));
	}
}

void ResetPointsForce(vector<Point*> Points) {
	for (int i = 0; i < Points.size(); i++) {
		Points[i]->Force[0] = 0.0;
		Points[i]->Force[1] = 0.0;

	}
}

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

/**********************************************************/


void display(GLfloat X, GLfloat Y) {
    
    
    
    
    //Square
    //        glVertex2f(1.0f, 1.0f);
    //        glVertex2f(1.0f, 0.0f);
    //        glVertex2f(0.0f, 0.0f);
    //        glVertex2f(0.0f, 1.0f);
    
    //    glVertex2f( (P.X - P.Size) - MoveDistance, (P.Y + P.Size));
    //    glVertex2f( (P.X - P.Size) - MoveDistance, (P.Y - P.Size));
    //    glVertex2f( (P.X + P.Size) - MoveDistance, (P.Y - P.Size));
    //    glVertex2f( (P.X + P.Size) - MoveDistance, (P.Y + P.Size));
    
    
    glVertex2f( (X - Size) - MoveDistance, (Y + Size));
    glVertex2f( (X - Size) - MoveDistance, (Y - Size));
    glVertex2f( (X + Size) - MoveDistance, (Y - Size));
    glVertex2f( (X + Size) - MoveDistance, (Y + Size));
    
    
    
    
    //this_thread::sleep_for(chrono::milliseconds(5));
    
    //MoveDistance = MoveDistance + 0.01f;
}

/*
void drawCircle( GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides )
{
    int numberOfVertices = numberOfSides + 2;
    
    GLfloat twicePi = 2.0f * M_PI;
    
    GLfloat circleVerticesX[numberOfVertices];
    GLfloat circleVerticesY[numberOfVertices];
    GLfloat circleVerticesZ[numberOfVertices];
    
    circleVerticesX[0] = x;
    circleVerticesY[0] = y;
    circleVerticesZ[0] = z;
    
    for ( int i = 1; i < numberOfVertices; i++ )
    {
        circleVerticesX[i] = x + ( radius * cos( i *  twicePi / numberOfSides ) );
        circleVerticesY[i] = y + ( radius * sin( i * twicePi / numberOfSides ) );
        circleVerticesZ[i] = z;
    }
    
    GLfloat allCircleVertices[( numberOfVertices ) * 3];
    
    for ( int i = 0; i < numberOfVertices; i++ )
    {
        allCircleVertices[i * 3] = circleVerticesX[i];
        allCircleVertices[( i * 3 ) + 1] = circleVerticesY[i];
        allCircleVertices[( i * 3 ) + 2] = circleVerticesZ[i];
    }
    
    glEnableClientState( GL_VERTEX_ARRAY );
    glVertexPointer( 3, GL_FLOAT, 0, allCircleVertices );
    glDrawArrays( GL_TRIANGLE_FAN, 0, numberOfVertices);
    glDisableClientState( GL_VERTEX_ARRAY );
}
*/