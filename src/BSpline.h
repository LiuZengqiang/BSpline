//
// Created by sth on 2021/1/21.
//

#ifndef BSPLINE_BSPLINE_H
#define BSPLINE_BSPLINE_H

#include "globalFunction.h"
#include "DataStruct.h"
#include "Eigen/Dense"
#include <iostream>
#include <algorithm>
#include <float.h>
#include <vector>

using namespace std;
using namespace Eigen;

class BSpline {
public:
    BSpline(int m, int n) {
        m_ = m;
        n_ = m;
    };

    ~BSpline() {};

    void init() {
        // todo:: init sphere data
//        data_points_.resize(m_);
//        float theta = 0.0f;
//        float alpha = 0.0f;
//        for (int i = 0; i <= m_; i++) {
//            data_points_[i].resize(n_ + 1);
//
//            for (int j = 0; j <= n_; j++) {
//
//
//
//
//            }
//        }

        data_points_.resize(m_ + 1);
        for (int i = 0; i < m_ + 1; i++) {
            data_points_[i].resize(n_ + 1);
            for (int j = 0; j < n_ + 1; j++) {
                float alpha = 180.0f * (i + 1) / (m_ + 2) - 90.0f;
                float theta = (360.0f) * j / (n_ + 1);
                float y = sin(globalFunction::angle2Radian(alpha));
                float temp = 1.0f - y * y;
                if (globalFunction::floatLessEqual(temp, 0.0f)) {
                    temp = 0.0f;
                }
                temp = sqrt(temp);
                float x = temp * sin(globalFunction::angle2Radian(theta));
                float z = temp * cos(globalFunction::angle2Radian(theta));
                data_points_[i][j].x = globalFunction::floatEqual(x, 0.0f) ? 0.0f : x;
                data_points_[i][j].x *= 100.0f;
                data_points_[i][j].y = globalFunction::floatEqual(y, 0.0f) ? 0.0f : y;
                data_points_[i][j].y *= 100.0f;
                data_points_[i][j].z = globalFunction::floatEqual(z, 0.0f) ? 0.0f : z;
                data_points_[i][j].z *= 100.0f;
            }
        }

        debugPoints(data_points_);
//        data_points_[0][0].x = -1.5;
//        data_points_[0][0].y = 1.5;
//        data_points_[3][0].x = -1.5;
//        data_points_[3][0].y = 1.5;

    }

    void calculateControlPointsInter() {

        // init s_, t_
        initST();
        // init u_, v_
        initUVInter();
        // init Nsi, Njt
        initBaseFunction();

        vector<MatrixXf> data_mat(3, MatrixXf(m_ + 1, n_ + 1));
        // init data_x_y_z;
        for (int i = 0; i < m_ + 1; i++) {
            for (int j = 0; j < n_ + 1; j++) {
                data_mat[0](i, j) = data_points_[i][j].x;
                data_mat[1](i, j) = data_points_[i][j].y;
                data_mat[2](i, j) = data_points_[i][j].z;
            }
        }
        //
        for (int i = 0; i < 3; i++) {
            MatrixXf Q = base_function_s_i_.colPivHouseholderQr().solve(data_mat[i]);
//            MatrixXf Q = base_function_s_i_.lu().solve(data_mat[i]);
            MatrixXf P = base_function_j_t_.transpose().colPivHouseholderQr().solve(Q.transpose());
//            MatrixXf P = base_function_j_t_.transpose().lu().solve(Q.transpose());
            data_mat[i] = P.transpose();
        }

        vector<vector<Point>> temp_vec_con(m_ + 1, vector<Point>(n_ + 1));
        control_points_.assign(temp_vec_con.begin(), temp_vec_con.end());

        for (int i = 0; i < m_ + 1; i++) {
            for (int j = 0; j < n_ + 1; j++) {
                control_points_[i][j].x = data_mat[0](i, j);
                control_points_[i][j].y = data_mat[1](i, j);
                control_points_[i][j].z = data_mat[2](i, j);
            }
        }

        calculateShowPoints();

        calError();
//        normalPoints(show_points_);
        normalPoints();
    };

    void calculateControlPointsAppro() {
        e_ = m_ - 1;
        f_ = n_ - 1;

        p_ = 3;
        q_ = 3;

//        p_ = 3;
//        q_ = 3;
        initST();

        initUVAppro();

        vector<vector<Point>> temp_vec_con(e_ + 1, vector<Point>(f_ + 1));
        // e_+1*f_+1
        control_points_.assign(temp_vec_con.begin(), temp_vec_con.end());

        for (int channel = 0; channel < 3; channel++) {
            MatrixXf D(data_points_.size(), data_points_[0].size());
            // e+1, n+1
            MatrixXf temp_D(e_ + 1, n_ + 1);

            for (int i = 0; i < data_points_.size(); i++) {
                for (int j = 0; j < data_points_[0].size(); j++) {
                    if (channel == 0) {
                        D(i, j) = data_points_[i][j].x;
                    } else if (channel == 1) {
                        D(i, j) = data_points_[i][j].y;
                    } else {
                        D(i, j) = data_points_[i][j].z;
                    }
                }
            }

            for (int col = 0; col < n_ + 1; col++) {
                // h_ -> e_;
                // P0 = D0
                float P0 = D(0, col);
                // Pe = Ph
                float Pe = D(m_, col);

                VectorXf P(e_ - 1);
                VectorXf QK(m_ - 1);
                VectorXf Q(e_ - 1);

                MatrixXf N(m_ - 1, e_ - 1);

                // init N(n-1*h-1)
                for (int i = 0; i < m_ - 1; i++) {
                    for (int j = 0; j < e_ - 1; j++) {
                        N(i, j) = getN(j + 1, p_, s_[i + 1], u_);
                    }
                }
                // init QK(h-1)
                for (int i = 0; i < m_ - 1; i++) {
                    int k = i + 1;
                    QK(i) = D(k, col) - getN(0, p_, s_[k], u_) * P0 - getN(e_, p_, s_[k], u_) * Pe;
                }
                // inti Q(h-1)
                for (int i = 0; i < e_ - 1; i++) {
                    float temp = 0.0f;
                    for (int j = 0; j < m_ - 1; j++) {
                        temp += QK(j) * getN(i + 1, p_, s_[j + 1], u_);
                    }
                    Q(i) = temp;
                }

//                Q = (QK.transpose() * N).transpose();
                // NT*N*P = Q
                MatrixXf NTN = N.transpose() * N;
                P = NTN.lu().solve(Q);
//                P = NTN.colPivHouseholderQr().solve(Q);

                temp_D(0, col) = P0;
                temp_D(e_, col) = Pe;

                for (int i = 0; i < e_ - 1; i++) {
                    temp_D(i + 1, col) = P(i);
                }
            }
            // temp_D(e_1 * n_+1)
            D = temp_D;

            // after processing column ,the row size is decreased
            for (int row = 0; row < e_ + 1; row++) {
                float P0 = D(row, 0);
                float Pf = D(row, n_);

                VectorXf P(f_ - 1);
                VectorXf QK(n_ - 1);
                VectorXf Q(f_ - 1);

                MatrixXf N(n_ - 1, f_ - 1);

                for (int i = 0; i < n_ - 1; i++) {
                    for (int j = 0; j < f_ - 1; j++) {
                        N(i, j) = getN(j + 1, q_, t_[i + 1], v_);
                    }
                }
                // init QK
                for (int i = 0; i < n_ - 1; i++) {
                    int k = i + 1;
                    QK(i) = D(row, k) - getN(0, q_, t_[k], v_) * P0 - getN(f_, q_, t_[k], v_) * Pf;
                }

                // inti Q(h-1)
                for (int i = 0; i < f_ - 1; i++) {
                    float temp = 0.0f;
                    for (int j = 0; j < n_ - 1; j++) {
                        temp += QK(j) * getN(i + 1, q_, t_[j + 1], v_);
                    }
                    Q(i) = temp;
                }
//                Q = (QK.transpose() * N).transpose();
                MatrixXf NTN = N.transpose() * N;
                // P()
                P = NTN.lu().solve(Q);
//                P = NTN.colPivHouseholderQr().solve(Q);
                if (channel == 0) {
                    control_points_[row][0].x = P0;
                    control_points_[row][f_].x = Pf;
                    for (int i = 0; i < f_ - 1; i++) {
                        control_points_[row][i + 1].x = P(i);
                    }
                } else if (channel == 1) {
                    control_points_[row][0].y = P0;
                    control_points_[row][f_].y = Pf;
                    for (int i = 0; i < f_ - 1; i++) {
                        control_points_[row][i + 1].y = P(i);
                    }
                } else {
                    control_points_[row][0].z = P0;
                    control_points_[row][f_].z = Pf;
                    for (int i = 0; i < f_ - 1; i++) {
                        control_points_[row][i + 1].z = P(i);
                    }
                }
            }
        }
        calculateShowPoints();
//        normalPoints(show_points_);
        calError();
        normalPoints();
    }

    vector<vector<Point>> &getShowPoints() {
        return show_points_;
    }

    vector<vector<Point>> &getDataPoints() {
        return data_points_;
    }

    vector<vector<Point>> &getControlPoints() {
        return control_points_;
    }

    void debugMat(MatrixXf &m) {
        cout << "Mat:" << endl;
        cout << m << endl;
    }

    void debugPoints(vector<vector<Point>> &points) {
        for (int i = 0; i < points.size(); i++) {
            for (int j = 0; j < points[i].size(); j++) {
                cout << "(" << points[i][j].x << "," << points[i][j].y << "," << points[i][j].z << ")" << " ";
            }
            cout << endl;
        }
    }

    float getError() {

        return error_dis_;
    }


private:
    /**
     * calculate surface points.
     */
    void calculateShowPoints() {

        show_points_.resize(num_points_on_line_ + 1);

        for (int i = 0; i <= num_points_on_line_; i++) {
            for (int j = 0; j <= num_points_on_line_; j++) {
                float s = ((float) i) / num_points_on_line_;
                float t = ((float) j) / num_points_on_line_;

                if (globalFunction::floatEqual(s, 1.0f)) {
                    s = 1.0f - Epsilon * 2;
                }
                if (globalFunction::floatEqual(t, 1.0f)) {
                    t = 1.0f - Epsilon * 2;
                }

                // s,t
                Point p(0.0f, 0.0f, 0.0f);
                for (int x = 0; x < control_points_.size(); x++) {
                    float n_i_p = getN(x, p_, s, u_);
                    for (int y = 0; y < control_points_[x].size(); y++) {
                        float n_j_q = getN(y, q_, t, v_);
                        p.x += n_i_p * n_j_q * control_points_[x][y].x;
                        p.y += n_i_p * n_j_q * control_points_[x][y].y;
                        p.z += n_i_p * n_j_q * control_points_[x][y].z;
                    }
                }
                show_points_[i].push_back(p);
            }
        }
    }

    /**
     * calculate base function N_i_p(t),N_j_q(t)...
     * @param i:i/j
     * @param p:p/q
     * @param t:s_[i]/t_[i]
     * @param u:knot vector
     * @return:N_i_p(t),N_j_q(t)
     */
    float getN(const int i, const int p, const float t, vector<float> &u) {
        if (p == 0) {
            return (t >= u[i] && t < u[i + 1]) ? 1.0f : 0.0f;
        }

        float left = u[i + p] - u[i];
        float right = u[i + p + 1] - u[i + 1];

        if (globalFunction::floatEqual(left, 0.0f)) {
//            left = 0.0f;
            left = (t - u[i]) * getN(i, p - 1, t, u);
        } else {
            left = (t - u[i]) / left * getN(i, p - 1, t, u);
        }
        if (globalFunction::floatEqual(right, 0.0f)) {
//            right = 0.0f;
            right = (u[i + p + 1] - t) * getN(i + 1, p - 1, t, u);
        } else {
            right = (u[i + p + 1] - t) / right * getN(i + 1, p - 1, t, u);
        }
        return left + right;
    }

    float calError() {
        float error = 0.0f;
        for (int i = 0; i < data_points_.size(); i++) {
            for (int j = 0; j < data_points_[i].size(); j++) {
                Point pre = data_points_[i][j];
                Point p(0.0f, 0.0f, 0.0f);
                float s = s_[i];
                float t = t_[j];
                if (globalFunction::floatEqual(s, 1.0f)) {
                    s = 1.0f - Epsilon * 2;
                }
                if (globalFunction::floatEqual(t, 1.0f)) {
                    t = 1.0f - Epsilon * 2;
                }
                for (int x = 0; x < control_points_.size(); x++) {
                    float n_i_p = getN(x, p_, s, u_);
                    for (int y = 0; y < control_points_[x].size(); y++) {
                        float n_j_q = getN(y, q_, t, v_);
                        p.x += n_i_p * n_j_q * control_points_[x][y].x;
                        p.y += n_i_p * n_j_q * control_points_[x][y].y;
                        p.z += n_i_p * n_j_q * control_points_[x][y].z;
                    }
                }
                error += globalFunction::getDis(pre, p);
            }
        }

        error_dis_ = error / (data_points_.size() * data_points_[0].size());
        error_dis_ = globalFunction::floatEqual(error_dis_, 0.0f) ? 0.0f : error_dis_;
    }

    /**
     * init s_,t_
     */
    void initST() {
        s_.resize(m_ + 1);
        vector<vector<float>> u(m_ + 1, vector<float>(n_ + 1, 0.0f));
        for (int j = 0; j < n_ + 1; j++) {
            float l = 0.0f;
            for (int i = 1; i < m_ + 1; i++) {
                float temp_l = globalFunction::getDis(data_points_[i - 1][j], data_points_[i][j]);
                l += temp_l;
                u[i][j] = l;
            }
            for (int i = 1; i < m_ + 1; i++) {
                u[i][j] /= l;
            }
            u[0][j] = 0.0f + Epsilon;
            u[m_][j] = 1.0f - Epsilon;
        }
        for (int i = 0; i < m_ + 1; i++) {
            float temp_s = 0.0f;
            for (int j = 0; j < n_ + 1; j++) {
                temp_s += u[i][j];
            }
            temp_s /= (n_ + 1);
            s_[i] = temp_s;
        }

        // init t
        t_.resize(n_ + 1);
        vector<vector<float>> v(m_ + 1, vector<float>(n_ + 1, 0.0f));
        for (int i = 0; i < m_ + 1; i++) {
            float l = 0.0f;
            for (int j = 1; j < n_ + 1; j++) {
                float temp_l = globalFunction::getDis(data_points_[i][j - 1], data_points_[i][j]);
                l += temp_l;
                v[i][j] = l;
            }
            for (int j = 1; j < n_ + 1; j++) {
                v[i][j] /= l;
            }
            v[i][0] = 0.0f + Epsilon;
            v[i][n_] = 1.0f - Epsilon;
        }

        for (int j = 0; j < n_ + 1; j++) {
            float temp_t = 0.0f;
            for (int i = 0; i < m_ + 1; i++) {
                temp_t += v[i][j];
            }
            temp_t /= m_ + 1;
            t_[j] = temp_t;
        }
        cout << "s_ t_:" << endl;
        for (int i = 0; i < m_ + 1; i++) {
            cout << s_[i] << " ";
        }
        cout << endl;
        for (int i = 0; i < n_ + 1; i++) {
            cout << t_[i] << " ";
        }
        cout << endl;
    }

    // init u,v
    void initUVAppro() {
        u_.resize(e_ + p_ + 1 + 1);
        v_.resize(f_ + q_ + 1 + 1);
        for (int i = 0; i < u_.size(); i++) {
            if (i <= p_) {
                u_[i] = 0.0f;
            } else if (i <= e_) {
                u_[i] = u_[i - 1] + 1.0f / (u_.size() - 2 * p_ - 2 + 1);
            } else {
                u_[i] = 1.0f;
            }
        }

        for (int i = 0; i < v_.size(); i++) {
            if (i <= q_) {
                v_[i] = 0.0f;
            } else if (i <= f_) {
                v_[i] = v_[i - 1] + 1.0f / (v_.size() - 2 * q_ - 2 + 1);
            } else {
                v_[i] = 1.0f;
            }
        }

        cout << "u_ v_:" << endl;
        for (int i = 0; i < u_.size(); i++) {
            cout << u_[i] << " ";
        }
        cout << endl;
        for (int i = 0; i < v_.size(); i++) {
            cout << v_[i] << " ";
        }
        cout << endl;

    }

    //init u_, v_
    void initUVInter() {
        u_.resize(m_ + p_ + 1 + 1);
        v_.resize(n_ + q_ + 1 + 1);
        for (int i = 0; i < u_.size(); i++) {
            if (i <= p_) {
                u_[i] = 0.0f;
            } else if (i <= m_) {
                float temp = 0.0f;
                for (int j = i - p_; j < i; j++) {
                    temp += s_[j];
                }
                u_[i] = 1.0f / p_ * temp;
            } else {
                u_[i] = 1.0f;
            }
        }

        for (int i = 0; i < v_.size(); i++) {
            if (i <= q_) {
                v_[i] = 0.0f;
            } else if (i <= n_) {
                float temp = 0.0f;
                for (int j = i - q_; j < i; j++) {
                    temp += t_[j];
                }
                v_[i] = 1.0f / q_ * temp;
            } else {
                v_[i] = 1.0f;
            }
        }
        cout << "u_ v_:" << endl;
        for (int i = 0; i < u_.size(); i++) {
            cout << u_[i] << " ";
        }
        cout << endl;
        for (int i = 0; i < v_.size(); i++) {
            cout << v_[i] << " ";
        }
        cout << endl;
    }

    // init base function
    void initBaseFunction() {
        base_function_s_i_.resize(m_ + 1, m_ + 1);
        for (int s = 0; s < m_ + 1; s++) {
            for (int i = 0; i < m_ + 1; i++) {
                float temp = getN(i, p_, s_[s], u_);
                temp = globalFunction::floatEqual(temp, 0.0f) ? 0.0f : temp;
                base_function_s_i_(s, i) = temp;
            }
        }

        base_function_j_t_.resize(n_ + 1, n_ + 1);
        for (int j = 0; j < n_ + 1; j++) {
            for (int t = 0; t < n_ + 1; t++) {
                base_function_j_t_(j, t) = getN(j, q_, t_[t], v_);
            }
        }
    }

    // translate data points,control points,show points to middle of windows
    void normalPoints() {
        normalPoints(data_points_);
        normalPoints(control_points_);
        normalPoints(show_points_);
    }

    void normalPoints(vector<vector<Point>> &points) {
        float max_x = -FLT_MAX;
        float min_x = FLT_MAX;
        float max_y = -FLT_MAX;
        float min_y = FLT_MAX;
        float max_z = -FLT_MAX;
        float min_z = FLT_MAX;
        for (int i = 0; i < points.size(); i++) {
            for (int j = 0; j < points[i].size(); j++) {
                max_x = max(max_x, points[i][j].x);
                max_y = max(max_y, points[i][j].y);
                max_z = max(max_z, points[i][j].z);

                min_x = min(min_x, points[i][j].x);
                min_y = min(min_y, points[i][j].y);
                min_z = min(min_z, points[i][j].z);
            }
        }

        for (int i = 0; i < points.size(); i++) {
            for (int j = 0; j < points[i].size(); j++) {
                if (!globalFunction::floatEqual(max_x, min_x)) {
                    points[i][j].x = (points[i][j].x - min_x) / (max_x - min_x) * 1.0f + (-0.5f);
                } else {
                    points[i][j].x = 0.0f;
                }
                if (!globalFunction::floatEqual(max_y, min_y)) {
                    points[i][j].y = (points[i][j].y - min_y) / (max_y - min_y) * 1.0f + (-0.5f);
                } else {
                    points[i][j].y = 0.0f;
                }
                if (!globalFunction::floatEqual(max_z, min_z)) {
                    points[i][j].z = (points[i][j].z - min_z) / (max_z - min_z) * 1.0f + (-0.5f);
                } else {
                    points[i][j].z = 0.0f;
                }
            }
        }
    }

    // for interpolation
    int n_ = 3;
    int m_ = 3;

    int p_ = 3;
    int q_ = 3;

    // for approximation
    // int m>e>=p>=1
    // int n>f>=q>=1
    // we define e = m-1, f=n-;
    // int m_ = m_;
    // int n_ = n_;
    int e_;
    int f_;
    float error_dis_ = 0.0f;

    // data points
    vector<vector<Point>> data_points_;
    // control points;
    vector<vector<Point>> control_points_;
    // show points
    vector<vector<Point>> show_points_;

    // base function
    MatrixXf base_function_s_i_;
    MatrixXf base_function_j_t_;

    // parameter vector
    vector<float> s_;
    vector<float> t_;
    // knots vector
    vector<float> u_;
    vector<float> v_;

    // show points size
    unsigned int num_points_on_line_ = 10;
};

#endif //BSPLINE_BSPLINE_H
