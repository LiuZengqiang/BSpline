// coding like poem

#include "GL/freeglut.h"
#include "GL/glu.h"
#include "./src/BSpline.h"
#include <vector>

using namespace std;
using namespace Eigen;
unsigned int SCR_WIDTH = 400;
unsigned int SCR_HEIGHT = 400;

vector<vector<Point>> show_points;

void init(void);

void display();

int main(int argc, char **argv) {


    BSpline bSpline;

    bSpline.init();

//    bSpline.calculateControlPointsInter();
    bSpline.calculateControlPointsAppro();

    show_points = bSpline.getShowPoints();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(SCR_WIDTH, SCR_HEIGHT);
    glutCreateWindow("B-Spline");

    init();

    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}

void init(void) {
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glShadeModel(GL_FLAT);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    for (int i = 0; i < show_points.size(); i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < show_points[i].size(); ++j) {
            glVertex3d(show_points[i][j].x, show_points[i][j].y, show_points[i][j].z);
        }
        glEnd();
    }

    for (int j = 0; j < show_points[0].size(); j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < show_points.size(); i++) {
            glVertex3d(show_points[i][j].x, show_points[i][j].y, show_points[i][j].z);
        }
        glEnd();
    }
    glutSwapBuffers();
}