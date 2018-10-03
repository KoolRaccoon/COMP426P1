#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <math.h>

#include <chrono>
#include <thread>

#include <unistd.h>
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
    
    GLfloat X;
    GLfloat Y;
    GLfloat Size = 0.005f;
    float Mass = 2.0f;
    vector<float> Force = {0, 0};
};


struct Node {
    vector<Point> PointsInNodeQuadrant;
    vector<Node> LeafNodes;
    GLfloat QuadrantSize;
    GLfloat OriginCoordinates[2];

    bool NodeHasOnlyOnePoint    = false;
    bool NodeIsEmpty            = false;

    int PlanetCount             = 0;
    float Mass                  = 0;
    float CenterOfMass[2];
    vector<float> Force         = {0, 0};
};

vector<vector<GLfloat>> PlanetCoordinates(10, vector<GLfloat>(2));
vector<Point> Points;
vector<Point*> Pointss;

GLfloat Size = 0.005f;

void drawCircle( GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides);


void display(vector<Point*>);
void GenerateRandomPoints();
void display(GLfloat, GLfloat);
void Tree(Node);
void ComputeMassDistribution(Node);
void CalculateForceOnPoint(Node);
vector<float> CalculateResultingForce(Node, Point);
const float G = 6.67398 * 0.00000000001;

int main() {
    srand(time(NULL));
    Point P;
    //cout << "Displaying Point P's X:" << P.X << endl;
    //cout << "Displaying Point P's Y:" << P.Y << endl;
   
    GenerateRandomPoints();
    Node Root;
    Root.QuadrantSize = 2.0f;
    Root.OriginCoordinates[0] = 0.0f;
    Root.OriginCoordinates[1] = 0.0f;
    
    for (int i = 0; i < Points.size(); i++){
        Root.PointsInNodeQuadrant.push_back(Points[i]);
        Root.PlanetCount++;
    }
    
    
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
        Tree(Root);
        glClearColor(0.2, 0.3, 0.4, 0.5);
        
        //clear color and depth buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        
        glColor3f(1, 1, 1);
        

//      Render here
        //for (int j=0; j<10; j++){
           // glBegin(GL_POLYGON);
            
            //P.X = PlanetCoordinates[j][0];
            //P.Y = PlanetCoordinates[j][1];
            
            //display(PlanetCoordinates[j][0], PlanetCoordinates[j][1]);
            display(Pointss);
            //glEnd();
            //glPopMatrix();
            
//        }

        //Swap front and back buffers
        
        glfwSwapBuffers(window);
        
        //Poll for and process events
        glfwPollEvents();
    }
    
    glfwTerminate();
    return 0;
    }


void GenerateRandomPoints(){

    for (int i=0; i<10; i++){
        Point * P = new Point;
        PlanetCoordinates[i][0] = ((float(rand()) / float(RAND_MAX)) * 1) + -0.5;
        PlanetCoordinates[i][1] = ((float(rand()) / float(RAND_MAX)) * 1) + -0.5;
        //cout << "displaying X: " << PlanetCoordinates[i][0] << "Displaying Y: " << PlanetCoordinates[i][1] << endl;
        P->X = PlanetCoordinates[i][0];
        P->Y = PlanetCoordinates[i][1];
        Pointss.push_back(P);
        //display(P);
    }
    
 
}

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
    
    
   
    
    this_thread::sleep_for(chrono::milliseconds(5));
    
    //MoveDistance = MoveDistance + 0.01f;
}


void display(vector<Point*> Points) {


    //glClearColor(0.2, 0.3, 0.4, 0.5);
    
    //clear color and depth buffer
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    
    //glColor3f(1, 1, 1);

    //glBegin(GL_POLYGON);
    

    //Square
//        glVertex2f(1.0f, 1.0f);
//        glVertex2f(1.0f, 0.0f);
//        glVertex2f(0.0f, 0.0f);
//        glVertex2f(0.0f, 1.0f);

//    glVertex2f( (P.X - P.Size) - MoveDistance, (P.Y + P.Size));
//    glVertex2f( (P.X - P.Size) - MoveDistance, (P.Y - P.Size));
//    glVertex2f( (P.X + P.Size) - MoveDistance, (P.Y - P.Size));
//    glVertex2f( (P.X + P.Size) - MoveDistance, (P.Y + P.Size));

    for (int j=0; j<10; j++){
        glBegin(GL_POLYGON);
        
        glVertex2f((Points[j]->X - Points[j]->Size), (Points[j]->Y + Points[j]->Size));
        glVertex2f((Points[j]->X - Points[j]->Size), (Points[j]->Y - Points[j]->Size));
        glVertex2f((Points[j]->X + Points[j]->Size), (Points[j]->Y - Points[j]->Size));
        glVertex2f((Points[j]->X + Points[j]->Size), (Points[j]->Y + Points[j]->Size));
        glEnd();
        glPopMatrix();
    }
    
    //glEnd();
    //glPopMatrix();

    //this_thread::sleep_for(chrono::milliseconds(33));
    
    MoveDistance = MoveDistance + 0.01f;
}


void Tree(Node Parent){
    
    if (Parent.PointsInNodeQuadrant.size() <= 1){
        if (Parent.PointsInNodeQuadrant.size() == 0)
            Parent.NodeIsEmpty = true;
        else
            Parent.NodeHasOnlyOnePoint = true;
    }
    else{
        
        Parent.NodeIsEmpty = false;
        Parent.NodeHasOnlyOnePoint = false;
        
        //Each leaf will represent a different quadrant. Leaf1 is Quadrant 1, Leaf2 is Quadrant 2, ... and so forth.
        Node Leaf1, Leaf2, Leaf3, Leaf4;
        //bool TreeCompleted = false;
        //while (TreeCompleted == false){
            //Calculating the new Quadrant Size of the leaf nodes

        Leaf1.QuadrantSize = Parent.QuadrantSize/2;
        Leaf2.QuadrantSize = Parent.QuadrantSize/2;
        Leaf3.QuadrantSize = Parent.QuadrantSize/2;
        Leaf4.QuadrantSize = Parent.QuadrantSize/2;
        
        //Finding New Origin Coordinates
        Leaf1.OriginCoordinates[0] = Parent.OriginCoordinates[0] + Parent.QuadrantSize/2;
        Leaf1.OriginCoordinates[1] = Parent.OriginCoordinates[1] + Parent.QuadrantSize/2;
        
        Leaf2.OriginCoordinates[0] = Parent.OriginCoordinates[0] - Parent.QuadrantSize/2;
        Leaf2.OriginCoordinates[1] = Parent.OriginCoordinates[1] + Parent.QuadrantSize/2;
        
        Leaf3.OriginCoordinates[0] = Parent.OriginCoordinates[0] - Parent.QuadrantSize/2;
        Leaf3.OriginCoordinates[1] = Parent.OriginCoordinates[1] - Parent.QuadrantSize/2;
        
        Leaf4.OriginCoordinates[0] = Parent.OriginCoordinates[0] + Parent.QuadrantSize/2;
        Leaf4.OriginCoordinates[1] = Parent.OriginCoordinates[1] - Parent.QuadrantSize/2;
        
        //Find the number of points in the quadrant.
        for (int j = 0; j < Parent.PointsInNodeQuadrant.size(); j++){
            //Populating each node's planet list with planets that belong in it's quadrant.
            
            if (Parent.PointsInNodeQuadrant[j].X >=Leaf1.OriginCoordinates[0] && Parent.PointsInNodeQuadrant[j].Y >= Leaf1.OriginCoordinates[1]){
                Leaf1.PointsInNodeQuadrant.push_back(Parent.PointsInNodeQuadrant[j]);
                Leaf1.PlanetCount++;
            }
            else if (Parent.PointsInNodeQuadrant[j].X < Leaf2.OriginCoordinates[0] && Parent.PointsInNodeQuadrant[j].Y >= Leaf2.OriginCoordinates[1]){
                Leaf2.PointsInNodeQuadrant.push_back(Parent.PointsInNodeQuadrant[j]);
                Leaf2.PlanetCount++;
            }
            else if (Parent.PointsInNodeQuadrant[j].X < Leaf3.OriginCoordinates[0] && Parent.PointsInNodeQuadrant[j].Y < Leaf3.OriginCoordinates[1]){
                Leaf3.PointsInNodeQuadrant.push_back(Parent.PointsInNodeQuadrant[j]);
                Leaf3.PlanetCount++;
            }
            else if (Parent.PointsInNodeQuadrant[j].X >= Leaf4.OriginCoordinates[0] && Parent.PointsInNodeQuadrant[j].Y < Leaf4.OriginCoordinates[1]){
                Leaf4.PointsInNodeQuadrant.push_back(Parent.PointsInNodeQuadrant[j]);
                Leaf4.PlanetCount++;
            }
        }
        
        if (Leaf1.PointsInNodeQuadrant.size() >= 1){
            Parent.LeafNodes.push_back(Leaf1);
            Tree(Leaf1);
        }
        if (Leaf2.PointsInNodeQuadrant.size() >= 1){
            Parent.LeafNodes.push_back(Leaf1);
            Tree(Leaf2);
        }
        if (Leaf3.PointsInNodeQuadrant.size() >= 1){
            Parent.LeafNodes.push_back(Leaf1);
            Tree(Leaf3);
        }
        if (Leaf4.PointsInNodeQuadrant.size() >= 1){
            Parent.LeafNodes.push_back(Leaf1);
            Tree(Leaf4);
        }
        //}
    }
}

//Computing the Mass of each Node and it's center of mass
void ComputeMassDistribution(Node Parent){
    
    if (Parent.NodeHasOnlyOnePoint == true){
        Parent.CenterOfMass[0] = Parent.PointsInNodeQuadrant[0].X;
        Parent.CenterOfMass[1] = Parent.PointsInNodeQuadrant[0].Y;
        Parent.Mass = Parent.PointsInNodeQuadrant[0].Mass;
    }
    else {
        for (int i = 0; i < Parent.LeafNodes.size(); i ++){
            ComputeMassDistribution(Parent.LeafNodes[i]);
            Parent.Mass += Parent.LeafNodes[i].Mass;
            Parent.CenterOfMass[0] += (Parent.LeafNodes[i].Mass * Parent.LeafNodes[i].CenterOfMass[0]);
            Parent.CenterOfMass[1] += (Parent.LeafNodes[i].Mass * Parent.LeafNodes[i].CenterOfMass[1]);        }
    }
    Parent.CenterOfMass[0] /= Parent.Mass;
    Parent.CenterOfMass[1] /= Parent.Mass;
}

//Computes the Total Forces acting on each Planet.
void CalculateForceOnPoint(Node Root){
    for ( int i = 0; i < Points.size(); i++){
        Points[i].Force = CalculateResultingForce(Root, Points[i]);
    }
    
}

// Computes the Total Force of the Tree that is acting upon a Certain Planet.
vector<float> CalculateResultingForce(Node Parent, Point TargetPlanet){
    TargetPlanet.Force[0] = 0;
    TargetPlanet.Force[1] = 0;
    vector<float> SumOfForces = {0, 0};
    int dis = 0;
    
    for (int i = 0; i < Parent.LeafNodes.size(); i ++){
        if (Parent.LeafNodes[i].PlanetCount == 1){
            //Compute Force between two planets
            float Force = 0;
            float Theta = 0;
            Theta = tan((TargetPlanet.X - Parent.LeafNodes[i].PointsInNodeQuadrant[0].X)/(TargetPlanet.Y - Parent.LeafNodes[i].PointsInNodeQuadrant[0].Y));
            dis   = sqrt(pow(2.0, (TargetPlanet.X - Parent.LeafNodes[i].PointsInNodeQuadrant[0].X)) + pow(2.0, (TargetPlanet.Y - Parent.LeafNodes[i].PointsInNodeQuadrant[0].Y)));
            Force = (G*TargetPlanet.Mass* Parent.LeafNodes[i].Mass)/(dis*dis);
            SumOfForces[0] += Force * cos(Theta);
            SumOfForces[1] += Force * sin(Theta);
            Theta = 0;
            dis   = 0;
            Force = 0;
        }
        else {
            int r = 0;
            int d = 0;
            r = sqrt(pow(2.0, (TargetPlanet.X - Parent.LeafNodes[i].CenterOfMass[0])) + pow(2.0, (TargetPlanet.Y - Parent.LeafNodes[i].CenterOfMass[1])));
            d = Parent.LeafNodes[i].QuadrantSize;
            if (d/r < 1){
                //Compute Force between Target Planet and Node
                float Force = 0;
                float Theta = 0;
                Theta = tan((TargetPlanet.X - Parent.LeafNodes[i].CenterOfMass[0])/(TargetPlanet.Y - Parent.LeafNodes[i].CenterOfMass[1]));
                dis   = sqrt(pow(2.0, (TargetPlanet.X - Parent.LeafNodes[i].CenterOfMass[0])) + pow(2.0, (TargetPlanet.Y - Parent.LeafNodes[i].CenterOfMass[1])));
                Force = (G*TargetPlanet.Mass* Parent.LeafNodes[i].Mass)/(dis*dis);
                SumOfForces[0] += Force * cos(Theta);
                SumOfForces[1] += Force * sin(Theta);
                Theta = 0;
                dis   = 0;
                Force = 0;

            }
            else {
                SumOfForces = CalculateResultingForce(Parent.LeafNodes[i], TargetPlanet);
                TargetPlanet.Force[0] += SumOfForces[0];
                TargetPlanet.Force[1] += SumOfForces[1];
            }
        }
    }
    return TargetPlanet.Force;
}




/*
 
 Function MainApp::CalcForce
    for all particles
        force = RootNode.CalculateForceFromTree(particle)
    end for
end
 
    Function force = TreeNode::CalculateForce(targetParticle)
    force = 0
 
     if number of particle equals 1
        force = Gravitational force between targetParticle and particle
     else
        r = distance from nodes center of mass to targetParticle
        d = height of the node
     if (d/r < θ)
        force = Gravitational force between targetParticle and node
     else
        for all child nodes n
        force += n.CalculateForce(particle)
     end for
 end if
 end
 end
 
 */






/**********************************************************/

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










