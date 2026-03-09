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


    // Compute the SVD decomposition of A
    Matrix U(m_rows, m_rows, 0.0);   // initialized with 0s
    Matrix S(m_rows, n_cols, 0.0);   // initialized with 0s
    Matrix V(n_cols, n_cols, 0.0);   // initialized with 0s

    svd_decompose(P, U, S, V);

    // Checks to see if SVD went correctly
    // Check 1: U is orthogonal, so U * U^T must be identity
    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
    std::cout << "P - U * S * V^T: \n" << P - U * S * transpose(V) << std::endl;

    // Forming the m  matrix by taking the last column of V matrix

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

    /*
    std::cout << " Cross Checking m Matrix by reprojecting \n";

    std::cout << " Points 3 D " << points_3d[0];

    Matrix test(4, 1, 0.0);
    test.set_column(0,{2,4,0,1});

    Matrix k = m *test;
    std::cout<<std::endl<< k << std::endl;

    // Divide by perspective
    float pers=k[2][0];

    std::cout<<pers<<std::endl;
    for (int i = 0;i<3;i++) {
        k[i][0]=k[i][0]/pers;
    }
    std::cout<< k << std::endl;

    */

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

    Vector3D b_vector(B[0][0], B[1][0], A[2][0]);

    double rho = 1 / (a_3.length());
    cx=rho*rho * (dot(a_1,a_3));
    cy=rho*rho * (dot(a_2,a_3));

    std::cout << "Cx :" << cx << std::endl << "Cy :" << cy << std::endl;

    double numerator = -1 * dot((cross(a_1,a_3)),(cross(a_2,a_3)));
    double denominator =  length(cross(a_1,a_3))*length(cross(a_2,a_3));
    double theta = acos(numerator/denominator);
    double theta_deg = theta * 180 / M_PI;

    std::cout << "Theta value : "<<theta << " Rad" <<  std::endl;
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

