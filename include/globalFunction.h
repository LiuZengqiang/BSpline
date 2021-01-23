//
// Created by sth on 2021/1/20.
//

#ifndef BSPLINE_GLOBALFUNCTION_H
#define BSPLINE_GLOBALFUNCTION_H

#include <stdlib.h>
#include "DataStruct.h"
#include <cmath>

#define Epsilon 1e-6

namespace globalFunction {

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

        return std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z));
    }

}
#endif //BSPLINE_GLOBALFUNCTION_H
