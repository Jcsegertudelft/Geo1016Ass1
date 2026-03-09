## About the Repository

This repository contains the implementation of a **camera calibration pipeline** developed as part of the GEO1016 course at TU Delft. 
The project focuses on estimating the parameters of a pinhole camera model using known correspondences between **3D world coordinates** and their **2D image projections**.

The code implements the mathematical steps required to recover the **camera projection matrix**, which describes how points in the real world are mapped onto the image plane of a camera. 
Using these correspondences, the algorithm constructs a system of linear equations and solves it using **Singular Value Decomposition (SVD)**. 
From the resulting projection matrix, the intrinsic camera parameters (which describe the internal properties of the camera) and the extrinsic parameters (which describe the camera’s position and orientation in space) are extracted.

The repository includes the full workflow required to perform camera calibration:

- Constructing the linear system from 3D–2D point correspondences
- Solving for the camera projection matrix using SVD
- Decomposing the projection matrix into intrinsic and extrinsic components
- Computing camera parameters such as focal lengths, principal point, and skew
- Estimating the rotation and translation of the camera in the world coordinate system
- Reprojecting the 3D points back into the image plane to evaluate calibration accuracy

To validate the implementation, the computed camera model is used to **reproject the input 3D points onto the image plane**, and the difference between the predicted and observed pixel coordinates is measured using the **reprojection error**. 
This provides a quantitative measure of how accurately the estimated camera parameters model the real imaging process.

Overall, this repository demonstrates the practical implementation of fundamental concepts in **computer vision and photogrammetry**, particularly the estimation of camera parameters from geometric constraints.
