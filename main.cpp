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

//
RenderMode renderMode = ALL;

// points for rendering
vector<vector<Point>> show_points;
vector<vector<Point>> data_points;
vector<vector<Point>> control_points;

void init(void);

void display();

void keyBoards(unsigned char key, int x, int y);

int main(int argc, char **argv) {
    int m = 10, n = 10;
    int num = 20;
    Mode mode = INTERPOLATION;

    // process parameters
    for (int i = 1; i < argc; i++) {
        string str(argv[i]);
        if (str.length() <= 1 || str[0] != '-') {
            cout << "parameter error: the " << i << " \'" << str
                 << "\' parameter is error, please use -m for m, -n for n, -M for fit model(approximation/interpolation), -g for grid size."
                 << endl;
        }
        if (str[1] == 'm') {
            if (i + 1 >= argc) {
                cout << "parameter error: the \'-m\' without value(unsigned int), please check the parameter list."
                     << endl;
            }
            string str_p = argv[++i];
            m = atoi(str_p.c_str());
            if (m <= 3) {
                cout << "parameter m should greater than 3, now the m used default value 4";
                m = 4;
            }
        } else if (str[1] == 'n') {
            if (i + 1 >= argc) {
                cout << "parameter error: the \'-n\' without value(unsigned int), please check the parameter list."
                     << endl;
            }
            string str_p = argv[++i];
            n = atoi(str_p.c_str());
            if (n <= 3) {
                cout << "parameter n should greater than 3, now the n used default value 4";
                n = 4;
            }
        } else if (str[1] == 'M') {
            if (i + 1 >= argc) {
                cout << "parameter error: the \'-M\' without value(bool), please check the parameter list." << endl;
            }
            string str_p = argv[++i];
            if (str_p == "i" || str_p == "I" || str_p == "interpolation" || str_p == "Interpolation" ||
                str_p == "INTERPOLATION") {
                mode = INTERPOLATION;
            } else {
                mode = APPROXIMATION;
            }
        } else if (str[1] == 'g') {
            if (i + 1 >= argc) {
                cout << "parameter error: the \'-g\' without value(unsigned int), please check the parameter list."
                     << endl;
            }
            string str_p = argv[++i];
            num = atoi(str_p.c_str());
        }
    }

    cout << "parameter list:" << endl;
    cout << "\tm:" << m << " n:" << n << " M:" << ((mode == INTERPOLATION) ? "interpolation" : "approximation") << " g:"
         << num << endl;

    BSpline bSpline(m, n, num);
    bSpline.init();
    if (mode == INTERPOLATION) {
        bSpline.calculateControlPointsInter();
    } else {
        bSpline.calculateControlPointsAppro();
    }
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
    glEnable(GL_DEPTH_TEST);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluLookAt(camera_pos_x, camera_pos_y, camera_pos_z, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
    glLoadIdentity();
    glPushMatrix();
    glRotatef(theta, 1.0f, 0.0f, 0.0f);
    glRotatef(alpha, 0.0f, 1.0f, 0.0f);
    // render data point
    if (renderMode == ALL || renderMode == DATAPOINT) {
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(5.0f);
        for (int i = 0; i < data_points.size(); i++) {
            for (int j = 0; j < data_points[i].size(); j++) {
                glBegin(GL_POINTS);
                glVertex3d(data_points[i][j].x, data_points[i][j].y, data_points[i][j].z);
                glEnd();
            }
        }
    }
    // render control points
    if (renderMode == ALL || renderMode == CONTROLPOINT) {
        glColor3f(0.0f, 1.0f, 0.0f);
        glPointSize(5.0f);
        for (int i = 0; i < control_points.size(); i++) {
            for (int j = 0; j < control_points[i].size(); j++) {
                glBegin(GL_POINTS);
                glVertex3d(control_points[i][j].x, control_points[i][j].y, control_points[i][j].z);
                glEnd();
            }
        }
    }

    if (renderMode == ALL || renderMode == SURFACE) {
        glColor3f(1.0f, 1.0f, 1.0f);
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

        for (int i = 0; i < show_points.size(); i++) {
            for (int j = 0; j < show_points[i].size(); j++) {
                if (i + 1 >= show_points.size() || j + 1 >= show_points[i].size()) {
                    break;
                }
                // four point render a rectangle
                Point p0 = show_points[i][j];
                Point p1 = show_points[i][j + 1];
                Point p2 = show_points[i + 1][j + 1];
                Point p3 = show_points[i + 1][j];
                Vector3f v0(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
                Vector3f v1(p3.x - p0.x, p3.y - p0.y, p3.z - p0.z);
                Vector3f nor = v0.cross(v1);
                nor = nor.normalized();

                glColor3f(nor(0) * 0.5f + 0.5f, nor(1) * 0.5f + 0.5f, nor(2) * 0.5f + 0.5f);
                glBegin(GL_TRIANGLE_FAN);
                glVertex3d(p0.x, p0.y, p0.z);
                glVertex3d(p1.x, p1.y, p1.z);
                glVertex3d(p2.x, p2.y, p2.z);
                glVertex3d(p3.x, p3.y, p3.z);
                glEnd();
            }
        }
        for (int j = 0; j < show_points[0].size(); j++) {
            glBegin(GL_LINE_STRIP);
            glVertex3d(0.0, -0.5, 0.0);
            glVertex3d(show_points.front()[j].x, show_points.front()[j].y, show_points.front()[j].z);
            glEnd();
            glBegin(GL_LINE_STRIP);
            glVertex3d(0.0, 0.5, 0.0);
            glVertex3d(show_points.back()[j].x, show_points.back()[j].y, show_points.back()[j].z);
            glEnd();
        }
    }
    glutSwapBuffers();
    glPopMatrix();
}

void keyBoards(unsigned char key, int x, int y) {
    float speed = 1.0f;

    if (key == 'W' || key == 'w') {
        theta += speed;
        glutPostRedisplay();
    } else if (key == 'S' || key == 's') {
        theta -= speed;
        glutPostRedisplay();
    } else if (key == 'A' || key == 'a') {
        alpha += speed;
        glutPostRedisplay();
    } else if (key == 'D' || key == 'd') {
        alpha -= speed;
        glutPostRedisplay();
    } else if (key == ' ') {
        switch (renderMode) {
            case ALL:
                renderMode = SURFACE;
                cout << "Surface." << endl;
                break;
            case SURFACE:
                renderMode = DATAPOINT;
                cout << "Data points." << endl;
                break;
            case DATAPOINT:
                renderMode = CONTROLPOINT;
                cout << "Control points." << endl;
                break;
            default:
                renderMode = ALL;
                cout << "All." << endl;
        }
        glutPostRedisplay();
    }

}