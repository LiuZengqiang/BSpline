//
// Created by sth on 2021/1/20.
//

#ifndef BSPLINE_GLOBALFUNCTION_H
#define BSPLINE_GLOBALFUNCTION_H

#include <stdlib.h>
#include "DataStruct.h"
#include <cmath>

#define Pi 3.1415926536
#define Epsilon 1e-6

namespace globalFunction {
    float angle2Radian(const float &angle) {
        return angle / 180.0f * Pi;
    }

    bool floatEqual(float a, float b) {
        return std::abs(a - b) <= Epsilon;
    }

    bool floatGreatEqual(float a, float b) {
        if (floatEqual(a, b) || a > b) {
            return true;
        } else {
            return false;
        }
    }

    bool floatLessEqual(float a, float b) {
        if (floatEqual(a, b) || a < b) {
            return true;
        } else {
            return false;
        }
    }

    float getDis(Point &a, Point &b) {
        float ret = std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z));
        return floatEqual(ret, 0.0f) ? 0.0f : ret;
    }

}
#endif //BSPLINE_GLOBALFUNCTION_H
