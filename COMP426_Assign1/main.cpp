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

using namespace std;

#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 480
GLfloat MoveDistance = 0;
GLFWwindow* window;

struct Point {
    
    GLfloat X;
    GLfloat Y;
    GLfloat Size = 0.005f;
    
};


struct Node {
    vector<Point> PointsInNodeQuadrant;
    GLfloat QuadrantSize;
    bool NodeHasOnlyOnePoint = false;
    vector<Node> LeafNodes;
    GLfloat OriginCoordinates[2];
    int PlanetCount = 0;
};

vector<vector<GLfloat>> PlanetCoordinates(10, vector<GLfloat>(2));
vector<Point> Points;

GLfloat Size= 0.02f;

void drawCircle( GLfloat x, GLfloat y, GLfloat z, GLfloat radius, GLint numberOfSides );


void display(Point);
void GenerateRandomPoints();
void display(GLfloat, GLfloat);
void Tree(Node);


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
        for (int j=0; j<10; j++){
            glBegin(GL_POLYGON);
            
            //P.X = PlanetCoordinates[j][0];
            //P.Y = PlanetCoordinates[j][1];
            
            display(PlanetCoordinates[j][0], PlanetCoordinates[j][1]);
            glEnd();
            glPopMatrix();
            
        }

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
        Point P;
        PlanetCoordinates[i][0] = ((float(rand()) / float(RAND_MAX)) * 2) + -1;
        PlanetCoordinates[i][1] = ((float(rand()) / float(RAND_MAX)) * 2) + -1;
        //cout << "displaying X: " << PlanetCoordinates[i][0] << "Displaying Y: " << PlanetCoordinates[i][1] << endl;
        P.X = PlanetCoordinates[i][0];
        P.Y = PlanetCoordinates[i][1];
        Points.push_back(P);
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


void display(Point P) {


    glClearColor(0.2, 0.3, 0.4, 0.5);
    
    //clear color and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    
    glColor3f(1, 1, 1);

    glBegin(GL_POLYGON);
    

    //Square
//        glVertex2f(1.0f, 1.0f);
//        glVertex2f(1.0f, 0.0f);
//        glVertex2f(0.0f, 0.0f);
//        glVertex2f(0.0f, 1.0f);

//    glVertex2f( (P.X - P.Size) - MoveDistance, (P.Y + P.Size));
//    glVertex2f( (P.X - P.Size) - MoveDistance, (P.Y - P.Size));
//    glVertex2f( (P.X + P.Size) - MoveDistance, (P.Y - P.Size));
//    glVertex2f( (P.X + P.Size) - MoveDistance, (P.Y + P.Size));

    
    glVertex2f( (P.X - P.Size), (P.Y + P.Size));
    glVertex2f( (P.X - P.Size), (P.Y - P.Size));
    glVertex2f( (P.X + P.Size), (P.Y - P.Size));
    glVertex2f( (P.X + P.Size), (P.Y + P.Size));
    
    
    glEnd();
    glPopMatrix();

    this_thread::sleep_for(chrono::milliseconds(33));
    
    MoveDistance = MoveDistance + 0.01f;
}


void Tree(Node Parent){
    

    for (int i = 0; i < Points.size(); i++)
        Parent.PointsInNodeQuadrant.push_back(Points[i]);
    
    if (Parent.PointsInNodeQuadrant.size() < 1){
        Parent.NodeHasOnlyOnePoint = true;
    }
    else{
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
            
            if (Parent.PointsInNodeQuadrant[j].X >=Leaf1.OriginCoordinates[0] && Parent.PointsInNodeQuadrant[j].Y >= Leaf1.OriginCoordinates[1])
                Leaf1.PointsInNodeQuadrant.push_back(Parent.PointsInNodeQuadrant[j]);
            else if (Parent.PointsInNodeQuadrant[j].X < Leaf2.OriginCoordinates[0] && Parent.PointsInNodeQuadrant[j].Y >= Leaf2.OriginCoordinates[1])
                Leaf2.PointsInNodeQuadrant.push_back(Parent.PointsInNodeQuadrant[j]);
            else if (Parent.PointsInNodeQuadrant[j].X < Leaf3.OriginCoordinates[0] && Parent.PointsInNodeQuadrant[j].Y < Leaf3.OriginCoordinates[1])
                Leaf3.PointsInNodeQuadrant.push_back(Parent.PointsInNodeQuadrant[j]);
            else if (Parent.PointsInNodeQuadrant[j].X >= Leaf4.OriginCoordinates[0] && Parent.PointsInNodeQuadrant[j].Y < Leaf4.OriginCoordinates[1])
                Leaf4.PointsInNodeQuadrant.push_back(Parent.PointsInNodeQuadrant[j]);
        }
        
        if (Leaf1.PointsInNodeQuadrant.size() > 1)
            Tree(Leaf1);
        if (Leaf2.PointsInNodeQuadrant.size() > 1)
            Tree(Leaf2);
        if (Leaf3.PointsInNodeQuadrant.size() > 1)
            Tree(Leaf3);
        if (Leaf4.PointsInNodeQuadrant.size() > 1)
            Tree(Leaf4);
            
        //}
        
    }
    
    
}












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










