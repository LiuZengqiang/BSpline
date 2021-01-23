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
    BSpline() {};

    ~BSpline() {};

    void initInter() {
        data_points_.resize(4);
        for (int i = 0; i < 4; i++) {
            data_points_[i].resize(4);
            for (int j = 0; j < 4; j++) {
                data_points_[i][j].x = (float) i;
                data_points_[i][j].y = (float) j;
                data_points_[i][j].z = 0.0f;
            }
        }
        data_points_[1][1].x = 1.5;
    }

    void calculateControlPointsInter() {

        m_ = data_points_.size() - 1;
        n_ = data_points_[0].size() - 1;

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

        show_points_.resize(num_points_on_line_);

        for (int i = 0; i < num_points_on_line_; i++) {

            for (int j = 0; j < num_points_on_line_; j++) {
                float s = ((float) i) / num_points_on_line_;
                float t = ((float) j) / num_points_on_line_;

                if (globalFunction::floatEqual(s, 1.0f)) {
                    s = 1.0f;
                }
                if (globalFunction::floatEqual(t, 1.0f)) {
                    t = 1.0f;
                }

                // s,t

                Point p(0.0f, 0.0f, 0.0f);

                for (int x = 0; x < m_ + 1; x++) {
                    float n_i_p = getN(x, p_, s, u_);

                    for (int y = 0; y < n_ + 1; y++) {
                        float n_j_q = getN(y, q_, t, v_);
                        p.x += n_i_p * n_j_q * control_points_[x][y].x;
                        p.y += n_i_p * n_j_q * control_points_[x][y].y;
                        p.z += n_i_p * n_j_q * control_points_[x][y].z;
                    }
                }
                show_points_[i].push_back(p);
            }
        }
        cout << "init show_points." << endl;

        normalShowPoints();
    };

    // todo::
    void calculateControlPointsAppro() {
        m_ = data_points_.size() - 1;
        n_ = data_points_[0].size() - 1;

        e_ = m_ - 1;
        f_ = n_ - 1;
        p_ = e_;
        q_ = f_;

        initST();

        initUVAppro();

//        initBaseFunction();
        // init control points' size

        vector<vector<Point>> temp_vec_con(e_ + 1, vector<Point>(f_ + 1));
        control_points_.assign(temp_vec_con.begin(), temp_vec_con.end());


        for (int channel = 0; channel < 3; channel++) {

            for (int col = 0; col < n_ + 1; col++) {
                // h_ -> e_;
                // P0 = D0
                float P0 = data_points_[0][col].x;
                float Pe = data_points_[m_][col].x;
                VectorXf P(e_ - 1);
                MatrixXf QK(e_ - 1);
                VectorXf Q(e_ - 1);
                MatrixXf N(e_ - 1, e_ - 1);

                // init N(h-1*h-1)
                for (int i = 0; i < e_ - 1; i++) {
                    for (int j = 0; j < e_ - 1; j++) {
                        N(i, j) = getN(j + 1, p_, s_[i + 1], u_);
                    }
                }
                // init QK
                for(int i=0; i<){

                }
            }
        }


    }

    vector<vector<Point>> &getShowPoints() {
        return show_points_;
    }

    void debugMat(MatrixXf &m) {
        cout << "Mat:" << endl;
        cout << m << endl;
    }

private:
    float getN(const int i, const int p, const float t, vector<float> &u) {
        if (p == 0) {
            return (t >= u[i] && t < u[i + 1]) ? 1.0f : 0.0f;
        }

        float left = u[i + p] - u[i];
        float right = u[i + p + 1] - u[i + 1];

        if (globalFunction::floatEqual(left, 0.0f)) {
            left = 0.0f;
            //left = (t - u[i]) * getN(i, p - 1, t, u);
        } else {
            left = (t - u[i]) / left * getN(i, p - 1, t, u);
        }
        if (globalFunction::floatEqual(right, 0.0f)) {
            right = 0.0f;
            //right = (u[i + p + 1] - t) * getN(i + 1, p - 1, t, u);
        } else {
            right = (u[i + p + 1] - t) / right * getN(i + 1, p - 1, t, u);
        }
        return left + right;
    }

    // init s,t
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

    void normalShowPoints() {
        float max_x = -FLT_MAX;
        float min_x = FLT_MAX;
        float max_y = -FLT_MAX;
        float min_y = FLT_MAX;
        float max_z = -FLT_MAX;
        float min_z = FLT_MAX;
        for (int i = 0; i < show_points_.size(); i++) {
            for (int j = 0; j < show_points_[i].size(); j++) {
                max_x = max(max_x, show_points_[i][j].x);
                max_y = max(max_y, show_points_[i][j].y);
                max_z = max(max_z, show_points_[i][j].z);

                min_x = min(min_x, show_points_[i][j].x);
                min_y = min(min_y, show_points_[i][j].y);
                min_z = min(min_z, show_points_[i][j].z);
            }
        }

        for (int i = 0; i < show_points_.size(); i++) {
            for (int j = 0; j < show_points_[i].size(); j++) {
                if (!globalFunction::floatEqual(max_x, min_x)) {
                    show_points_[i][j].x = (show_points_[i][j].x - min_x) / (max_x - min_x) * 1.0f + (-0.5f);
                } else {
                    show_points_[i][j].x = 0.0f;
                }
                if (!globalFunction::floatEqual(max_y, min_y)) {
                    show_points_[i][j].y = (show_points_[i][j].y - min_y) / (max_y - min_y) * 1.0f + (-0.5f);
                } else {
                    show_points_[i][j].y = 0.0f;
                }
                if (!globalFunction::floatEqual(max_z, min_z)) {
                    show_points_[i][j].z = (show_points_[i][j].z - min_z) / (max_z - min_z) * 1.0f + (-0.5f);
                } else {
                    show_points_[i][j].z = 0.0f;
                }
            }
        }
    }

    // for interpolation
    int n_ = -1;
    int m_ = -1;

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

    unsigned int num_points_on_line_ = 10;

    vector<float> s_;
    vector<float> t_;
    vector<float> u_;
    vector<float> v_;
};

#endif //BSPLINE_BSPLINE_H