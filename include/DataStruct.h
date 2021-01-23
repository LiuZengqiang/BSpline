//
// Created by sth on 2021/1/22.
//

#ifndef BSPLINE_DATASTRUCT_H
#define BSPLINE_DATASTRUCT_H
enum Mode {
    INTERPOLATION,
    APPROXIMATION
};


struct Point {
    int id = -1;
    float x;
    float y;
    float z;

    Point() {
        x = y = z = 0.0f;
    }

    Point(float x, float y, float z) : x(x), y(y), z(z) {}
};

#endif //BSPLINE_DATASTRUCT_H
