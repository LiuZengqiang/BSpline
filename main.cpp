// coding like poem

#include "GL/freeglut.h"
#include "GL/glu.h"
#include "./src/BSpline.h"
#include <vector>

using namespace std;
using namespace Eigen;

// windows size
unsigned int SCR_WIDTH = 400;
unsigned int SCR_HEIGHT = 400;

float camera_pos_x = 0.0f;
float camera_pos_y = 0.0f;
float camera_pos_z = 1.0f;
float theta = 0.0f;
float alpha = 0.0f;

// points for rendering
vector<vector<Point>> show_points;
vector<vector<Point>> data_points;
vector<vector<Point>> control_points;

void init(void);

void display();

void keyBoards(unsigned char key, int x, int y);

int main(int argc, char **argv) {

    int m = 10, n = 10;
    Mode mode = INTERPOLATION;
    BSpline bSpline(m, n);
    bSpline.init();
//    bSpline.calculateControlPointsInter();
    bSpline.calculateControlPointsAppro();

    show_points = bSpline.getShowPoints();
    data_points = bSpline.getDataPoints();
    control_points = bSpline.getControlPoints();
    float error = bSpline.getError();

    cout << "error is:" << error << endl;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(SCR_WIDTH, SCR_HEIGHT);
    glutCreateWindow("B-Spline");
    init();
    glutKeyboardFunc(keyBoards);
    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}

void init(void) {
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
}

// todo::add show data points, control points
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluLookAt(camera_pos_x, camera_pos_y, camera_pos_z, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    glLoadIdentity();
    glPushMatrix();
    glRotatef(theta, 1.0f, 0.0f, 0.0f);
    glRotatef(alpha, 0.0f, 1.0f, 0.0f);
    // glLineWidth(10);

    // render data point
    glColor3f(255.0f, 0.0f, 0.0f);
    glPointSize(5.0f);
    for (int i = 0; i < data_points.size(); i++) {
        for (int j = 0; j < data_points[i].size(); j++) {
            glBegin(GL_POINTS);
            glVertex3d(data_points[i][j].x, data_points[i][j].y, data_points[i][j].z);
            glEnd();
        }
    }
    // render control points
    glColor3f(0.0f, 255.0f, 0.0f);
    for (int i = 0; i < control_points.size(); i++) {
        for (int j = 0; j < control_points[i].size(); j++) {
            glBegin(GL_POINTS);
            glVertex3d(control_points[i][j].x, control_points[i][j].y, control_points[i][j].z);
            glEnd();
        }
    }

//    for (int i = 0; i < control_points.size(); i++) {
//        glBegin(GL_LINE_STRIP);
//        for (int j = 0; j < control_points[i].size(); j++) {
//            glVertex3d(control_points[i][j].x, control_points[i][j].y, control_points[i][j].z);
//        }
//        glEnd();
//    }
//    for (int j = 0; j < control_points[0].size(); j++) {
//        glBegin(GL_LINE_STRIP);
//        for (int i = 0; i < control_points.size(); ++i) {
//            glVertex3d(control_points[i][j].x, control_points[i][j].y, control_points[i][j].z);
//        }
//        glEnd();
//    }

    glColor3f(255.0f, 255.0f, 255.0f);
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
    glPopMatrix();
}

void keyBoards(unsigned char key, int x, int y) {
    float speed = 1.0f;
    switch (key) {
        case 'w':
            theta += speed;
            glutPostRedisplay();
            break;
        case 's':
            theta -= speed;
            glutPostRedisplay();
            break;
        case 'a':
            alpha += speed;
            glutPostRedisplay();
            break;
        case 'd':
            alpha -= speed;
            glutPostRedisplay();
            break;
        default:
            break;
    }
}