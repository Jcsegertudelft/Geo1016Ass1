/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"

using namespace easy3d;

bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,  /// output: focal length (i.e., K[0][0]).
        double& fy,  /// output: focal length (i.e., K[1][1]).
        double& cx,  /// output: x component of the principal point (i.e., K[0][2]).
        double& cy,  /// output: y component of the principal point (i.e., K[1][2]).
        double& s,   /// output: skew factor (i.e., K[0][1]), which is s = -alpha * cot(theta).
        Matrix33& R, /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t) /// output：a 3D vector encoding camera translation.
{

    // check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    unsigned int n_points_3d = points_3d.size();
    unsigned int n_points_2d = points_2d.size();
    if (n_points_3d < 6 || n_points_2d != n_points_3d) {
        return false;
    }

    // construct the P matrix (so P * m = 0).
    Matrix P (2 * n_points_3d, 12);
    for (unsigned int i = 0; i < n_points_3d; i++){
        unsigned int i_row = 2*i;
        Vector3D point_world = points_3d[i];
        Vector2D point_picture = points_2d[i];
        P.set_row(i_row, {
            point_world[0],
            point_world[1],
            point_world[2],
            1,
            0,
            0,
            0,
            0,
            -point_picture[0]*point_world[0],
            -point_picture[0]*point_world[1],
            -point_picture[0]*point_world[2],
            -point_picture[0]
        });
        P.set_row(i_row + 1, {
            0,
            0,
            0,
            0,
            point_world[0],
            point_world[1],
            point_world[2],
            1,
            -point_picture[1]*point_world[0],
            -point_picture[1]*point_world[1],
            -point_picture[1]*point_world[2],
            -point_picture[1]});
    }

    //solving for M where M = K * [ R, t] using SVD decomposition

    // From the obtained P matrix, Pm=0.
    // m -> Rows
    // n -> Columns

    int m_rows = 2 * n_points_3d;
    int n_cols = 12;


    // Computing the the SVD decomposition of Matrix P
    Matrix U(m_rows, m_rows, 0.0);   // initialized with 0s
    Matrix S(m_rows, n_cols, 0.0);   // initialized with 0s
    Matrix V(n_cols, n_cols, 0.0);   // initialized with 0s

    svd_decompose(P, U, S, V);

    int count = 0;
    Matrix34 m;
    for (int i = 0;i<3;i++) {
        for (int j = 0;j<4;j++) {
            m[i][j] = V(count,11);
            count++;
        }
        if (count > 11)
            break;

    }

    std::cout<<std::endl << "m Matrix : "<< m << std::endl;

    // Evaluation of the Calibration
    std::cout << std::endl << " Cross Checking m Matrix by reprojecting the points and calculating error: \n";
    std::cout <<"Point\t"<<"U(Ori.)\t"<<" V(Ori.)\t"<<" u(Repr.)\t"<<" v(Repr. )\t"<<"Error"<<std::endl;
    double total_error = 0.0;
    for (unsigned int i = 0; i < n_points_3d; i++)
    {
        // 3D point in homogeneous coordinates
        Matrix X(4,1,0.0);
        X.set_column(0,{
            points_3d[i][0],
            points_3d[i][1],
            points_3d[i][2],
            1
        });
        // Projection using Matrix m
        Matrix proj = m * X;

        // Perspective division
        double u = proj[0][0] / proj[2][0];
        double v = proj[1][0] / proj[2][0];

        // Original image point
        double u_original = points_2d[i][0];
        double v_original = points_2d[i][1];

        // Reprojection error
        double err = sqrt((u - u_original)*(u - u_original) + (v - v_original)*(v - v_original));
        total_error = err + total_error;
        std::cout<<"\t"<<i<<"\t\t"<<u_original<<"\t\t"<<v_original<<"\t\t"<<u<<"\t\t"<<v<<"\t\t"<<err<<std::endl;
    }
    double mean_error = total_error / n_points_3d;

    std::cout <<"\n Mean Reprojection Error: " << mean_error << std::endl;


    // Denoting the  M matrix into [ A b ] format i.e A = 3 x 3, B = 3 x 1 Matrices

    Matrix A(3,3,0.0);
    Matrix B(3,1,0.0);

    for (int i = 0;i<3;i++) {
        for (int j = 0;j<3;j++) {
            A(i,j) = m[i][j];
        }
    }

    for (int i=0;i<3;i++) {
        B(i,0) = m[i][3];
    }


    Vector3D a_1(A[0][0], A[0][1], A[0][2]);
    Vector3D a_2(A[1][0], A[1][1], A[1][2]);
    Vector3D a_3(A[2][0], A[2][1], A[2][2]);

    Vector3D b_vector(B[0][0], B[1][0], B[2][0]);

    double rho = 1 / (a_3.length());
    cx=rho*rho * (dot(a_1,a_3));
    cy=rho*rho * (dot(a_2,a_3));

    std::cout << "Cx :" << cx << std::endl << "Cy :" << cy << std::endl;

    double numerator = -1 * dot((cross(a_1,a_3)),(cross(a_2,a_3)));
    double denominator =  length(cross(a_1,a_3))*length(cross(a_2,a_3));
    double c = numerator/denominator;
    c = std::max(-1.0, std::min(1.0, c));
    double theta = acos(c);    
    double theta_deg = theta * 180 / M_PI;

    std::cout << "Theta value : "<< theta << " Rad" <<  std::endl;
    std::cout << "Theta in degrees : "<< theta_deg << std::endl;

    double alpha = rho * rho * length(cross(a_1,a_3)) * sin(theta);
    double beta = rho * rho * length(cross(a_2,a_3)) * sin(theta);

    std::cout << "Alpha : " << alpha << std::endl << "Beta : "<< beta << std::endl;

    // Making the intrinsic matrix
    Matrix33 K;

    K.set_row(0,{alpha,-1*alpha * (1/tan(theta)), cx});
    K.set_row(1,{0,beta/sin(theta),cy});
    K.set_row(2,{0,0,1});

    std::cout << "\n Matrix K : "<< K;

    // Calculating the extrinsic parameters from A and b

    Vector3D r_1 = (cross(a_2,a_3))/length((cross(a_2,a_3)));
    Vector3D r_3 = rho * a_3;
    Vector3D r_2 = cross(r_3,r_1);

    R.set_row(0,r_1);
    R.set_row(1,r_2);
    R.set_row(2,r_3);

     Matrix h = rho *  inverse(K) * B;
     Vector3D temp(h[0][0], h[1][0], h[2][0]);

    //Checking for sign of rho through Z-value of translation
    if (temp.z() < 0)
    {
        rho = -rho;
        temp = -temp;
        r_3 = -r_3;
        r_2 = -r_2;
        R.set_row(2,r_3);
        R.set_row(1,r_2);
    }
     t = temp;

    std::cout << "rho : " << rho << std::endl;

     std::cout<<"\n Matrix R : " << R << std::endl;

     std::cout << "t : " << t << std::endl;

     fx = alpha;
     fy = beta/sin(theta);
     s = -1 * alpha * (1/tan(theta));

    std::cout << "fx : " << fx << std::endl << "fy : " << fy << std::endl << "s : " << s << std::endl;




    return true;
}

