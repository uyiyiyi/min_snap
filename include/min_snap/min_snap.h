#include <ros/ros.h>
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include "OsqpEigen/OsqpEigen.h"
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/Twist.h>
#include <cmath>
#include <stdexcept>


using namespace std;
using namespace Eigen;


template<typename _Matrix_Type_> 
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = 
    numeric_limits<double>::epsilon()) 
{  
    JacobiSVD< _Matrix_Type_ > svd(a ,ComputeThinU | ComputeThinV);  
    double tolerance = epsilon * max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);  
    return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint(); 
}   

int factorial(int n)    
{    
    if(n < 0)    
        return(-1); /*Wrong value*/      
    if(n == 0)    
        return(1);  /*Terminating condition*/    
    else    
    {    
        return(n*factorial(n-1));        
    }
}


class min_snap
{
private:
    ros::Subscriber waypt_path_sub, end_of_path_sub;
    // ros::Publisher traj_pub;
    ros::Publisher traj_pub;
    static const int wp_num = 5;
    double time[wp_num];
    double max_velocity = 10;
    double total_time = 5;
    int derivative_order, seg_num, p_order, p_num, m, n;
public:
    MatrixXd path = MatrixXd::Zero(3, wp_num);
    VectorXd Solution_x, Solution_y, Solution_z;
    // path.push_back(wp1);
    double v0[2] = {0, 0};
    double a0[2] = {0, 0};
    double v1[2] = {0, 0};
    double a1[2] = {0, 0};
    double T = 5;
    
    min_snap(ros::NodeHandle nh);
    ~min_snap();
    void solve_bAp();
    void solve_Nseg_bAp();
    void uniform_time_arrange(MatrixXd path_);
    int solveQP(MatrixXd Hessian);
    void visTraj(VectorXd Sol_x, VectorXd Sol_y, VectorXd Sol_z);
    void waypoint_cb(const geometry_msgs::Twist::ConstPtr & msg);
    double PolyEval(const VectorXd& coeffs, double x);
    Vector3d CalcVector(const VectorXd& coeffsX, const VectorXd& coeffsY, const VectorXd& coeffsZ, double t);
};

min_snap::min_snap(ros::NodeHandle nh)
{
    waypt_path_sub = nh.subscribe("waypoint", 1000, &min_snap::waypoint_cb, this);
    traj_pub = nh.advertise<nav_msgs::Path>("trajectory", 1);

    path.col(0) = Vector3d(0, 0, 0);
    path.col(1) = Vector3d(1, 2, 0);
    path.col(2) = Vector3d(2, -1, 0);
    path.col(3) = Vector3d(4, 8, 0);
    path.col(4) = Vector3d(5, 2, 0);
    uniform_time_arrange(path);
    solve_Nseg_bAp();
    // visTraj(Solution_x, Solution_y, Solution_z);

    // solve_bAp();
    // cout << " path: " << endl << path << endl;
    // path.push_back(wp1);
    // path.push_back(wp2);
    // path.push_back(wp3);
    // path.push_back(wp4);
    // path.push_back(wp5);
    // cout << " path: " << path.size() << endl;
    // time_arrange(path);
}

min_snap::~min_snap()
{
}

void min_snap::solve_Nseg_bAp()
{
    // p(t) = p0 + p1 * t + p2 * t2 + ... + pn * tn = ∑ pi * ti
    derivative_order = 3; // pos = 0, vel = 1, acc = 2, jerk = 3, snap = 4;
    seg_num = path.row(0).size() - 1;
    p_order = 2 * derivative_order - 1; // Polynomial order, for jerk is 5, for snap is 7
    p_num = p_order + 1;
    m = 4 * wp_num - 2;
    n = p_num * seg_num;
    MatrixXd Q = MatrixXd::Zero(p_num * seg_num, p_num * seg_num);
    MatrixXd Qi =MatrixXd::Zero(p_num, p_num);
    VectorXd t = VectorXd::Ones(p_num);

    for(int seg = 0; seg < seg_num; seg++)
    {
        for (int i = derivative_order; i < p_num; i++)
        {
            for (int j = derivative_order; j < p_num; j++)
            {
                // Qi(i, j) = (factorial(i + 1) * factorial(j + 1) * pow(time[seg + 1] - time[seg], i + j + 3 - 2 * derivative_order)) / (factorial(i + 1 - derivative_order) * factorial(j + 1 - derivative_order) * (i + j + 3 - 2 * derivative_order));
                Qi(i, j) = (factorial(i + 1) * factorial(j + 1) * (pow(time[seg + 1], i + j + 3 - 2 * derivative_order) - pow(time[seg], i + j + 3 - 2 * derivative_order))) / (factorial(i + 1 - derivative_order) * factorial(j + 1 - derivative_order) * (i + j + 3 - 2 * derivative_order));
            }
        }
        Q.block(seg * p_num, seg * p_num, p_num, p_num) = Qi;
    }

    SparseMatrix<double> hessian(n, n);      //P: n*n正定矩阵,必须为稀疏矩阵SparseMatrix
    hessian = Q.sparseView();
    // cout << hessian << endl;
    
    // LLT<MatrixXd> lltOfA(Q); // compute the Cholesky decomposition of A 
    // cout << lltOfA.info() << endl;
    // if(lltOfA.info() == NumericalIssue) 
    // { 
    //     throw runtime_error("Possibly non semi-positive definitie matrix!"); 
    // }

    
    VectorXd gradient = VectorXd::Zero(n);                    //Q: n*1向量
    SparseMatrix<double> linearMatrix(m, n); //A: m*n矩阵,必须为稀疏矩阵SparseMatrix
    VectorXd lowerBound_x(m), lowerBound_y(m), lowerBound_z(m);                  //L: m*1下限向量
    VectorXd upperBound_x(m), upperBound_y(m), upperBound_z(m);                  //U: m*1上限向量
    
    MatrixXd Aeq = MatrixXd::Zero(m, n);
    MatrixXd Ai = MatrixXd::Zero(1, p_num);
    MatrixXd Ai_ = MatrixXd::Zero(1, p_num);
    for(int i = 0; i < seg_num; i++)
    {
        for(int j = 0; j < p_num; j++)
        {
            Ai(0, j) = pow(time[i], j);
            Ai_(0, j) = pow(time[i + 1], j);
            
        }
        Aeq.block(2 * i, p_num * i, 1, p_num) = Ai;
        Aeq.block(2 * i + 1, p_num * i, 1, p_num) = Ai_;
        lowerBound_x[2 * i] = path.col(i).x();
        upperBound_x[2 * i] = path.col(i).x();
        lowerBound_x[2 * i + 1] = path.col(i + 1).x();
        upperBound_x[2 * i + 1] = path.col(i + 1).x();
        
        lowerBound_y[2 * i] = path.col(i).y();
        upperBound_y[2 * i] = path.col(i).y();
        lowerBound_y[2 * i + 1] = path.col(i + 1).y();
        upperBound_y[2 * i + 1] = path.col(i + 1).y();
        
        lowerBound_z[2 * i] = path.col(i).z();
        upperBound_z[2 * i] = path.col(i).z();
        lowerBound_z[2 * i + 1] = path.col(i + 1).z();
        upperBound_z[2 * i + 1] = path.col(i + 1).z();
        
    }
    MatrixXd Veli = MatrixXd::Zero(1, 2 * p_num);
    MatrixXd Acci = MatrixXd::Zero(1, 2 * p_num);
    for(int i = 0; i < wp_num; i++)
    {
        for (int j = 0; j < p_num; j++)
        {
            if(j == 0)
            {
                Veli(0, j) = 0;
                Veli(0, j + p_num) = 0;
            }
            else
            {
                Veli(0, j) = j * pow(time[i], j - 1);
                Veli(0, j + p_num) = -j * pow(time[i], j - 1);
                if(i == 0)
                {
                    Veli(0, j) = j * pow(time[i], j - 1);
                    Veli(0, j + p_num) = 0;
                }
                if(i == wp_num - 1)
                {
                    Veli(0, j) = 0;
                    Veli(0, j + p_num) = j * pow(time[i], j - 1);
                }
            }
            if(j == 0 || j == 1)
            {
                Acci(0, j) = 0;
                Acci(0, j + p_num) = 0;
            }
            else
            {
                Acci(0, j) = j * (j - 1) * pow(time[i], j - 2);
                Acci(0, j + p_num) = -j * (j - 1) * pow(time[i], j - 2);
                if(i == 0)
                {
                    Acci(0, j) = j * (j - 1) * pow(time[i], j - 2);
                    Acci(0, j + p_num) = 0;
                }
                if(i == wp_num - 1)
                {
                    Acci(0, j) = 0;
                    Acci(0, j + p_num) = j * (j - 1) * pow(time[i], j - 2);
                }
            }
        }
        if(i == 0)
        {

            Aeq.block(2 * seg_num, 0, 1, 2 * p_num) = Veli;
            Aeq.block(2 * seg_num + 1, 0, 1, 2 * p_num) = Acci;
            lowerBound_x[2 * seg_num] = 0;
            upperBound_x[2 * seg_num] = 0;
            lowerBound_x[2 * seg_num + 1] = 0;
            upperBound_x[2 * seg_num + 1] = 0;
            lowerBound_y[2 * seg_num] = 0;
            upperBound_y[2 * seg_num] = 0;
            lowerBound_y[2 * seg_num + 1] = 0;
            upperBound_y[2 * seg_num + 1] = 0;
            lowerBound_z[2 * seg_num] = 0;
            upperBound_z[2 * seg_num] = 0;
            lowerBound_z[2 * seg_num + 1] = 0;
            upperBound_z[2 * seg_num + 1] = 0;
        }
        else if(i == wp_num - 1)
        {
            Aeq.block(m - 2, p_num * (seg_num - 2), 1, 2 * p_num) = Veli;
            Aeq.block(m - 1, p_num * (seg_num - 2), 1, 2 * p_num) = Acci;
            lowerBound_x[m - 2] = 0;
            upperBound_x[m - 2] = 0;
            lowerBound_x[m - 1] = 0;
            upperBound_x[m - 1] = 0;
            lowerBound_y[m - 2] = 0;
            upperBound_y[m - 2] = 0;
            lowerBound_y[m - 1] = 0;
            upperBound_y[m - 1] = 0;
            lowerBound_z[m - 2] = 0;
            upperBound_z[m - 2] = 0;
            lowerBound_z[m - 1] = 0;
            upperBound_z[m - 1] = 0;
        }
        else
        {
            Aeq.block(2 * i + 2 * seg_num, p_num * (i - 1), 1, 2 * p_num) = Veli;
            Aeq.block(2 * i + 2 * seg_num + 1, p_num * (i - 1), 1, 2 * p_num) = Acci;
            lowerBound_x[2 * i + 2 * seg_num] = 0;
            upperBound_x[2 * i + 2 * seg_num] = 0;
            lowerBound_x[2 * i + 2 * seg_num + 1] = 0;
            upperBound_x[2 * i + 2 * seg_num + 1] = 0;
            lowerBound_y[2 * i + 2 * seg_num] = 0;
            upperBound_y[2 * i + 2 * seg_num] = 0;
            lowerBound_y[2 * i + 2 * seg_num + 1] = 0;
            upperBound_y[2 * i + 2 * seg_num + 1] = 0;
            lowerBound_z[2 * i + 2 * seg_num] = 0;
            upperBound_z[2 * i + 2 * seg_num] = 0;
            lowerBound_z[2 * i + 2 * seg_num + 1] = 0;
            upperBound_z[2 * i + 2 * seg_num + 1] = 0;
        }
    }
    // cout << Aeq << endl;
    linearMatrix = Aeq.sparseView();
    // cout << linearMatrix << endl;
    // cout << "lowerBound_x = " << lowerBound_x << endl;
    // cout << "lowerBound_y = " << lowerBound_y << endl;
    // cout << "lowerBound_z = " << lowerBound_z << endl;


    // instantiate the solver
    OsqpEigen::Solver solver_x, solver_y, solver_z;
    // settings
    solver_x.settings()->setVerbosity(false);
    solver_x.settings()->setWarmStart(true);

    // set the initial data of the QP solver
    solver_x.data()->setNumberOfVariables(n);   //变量数n
    solver_x.data()->setNumberOfConstraints(m); //约束数m
    solver_x.data()->setHessianMatrix(hessian);
    solver_x.data()->setGradient(gradient);
    solver_x.data()->setLinearConstraintsMatrix(linearMatrix);
    solver_x.data()->setLowerBound(lowerBound_x);
    solver_x.data()->setUpperBound(upperBound_x);

    // instantiate the solver
    solver_x.initSolver();
    

    // solve the QP problem
    solver_x.solveProblem();
    Solution_x = solver_x.getSolution();
    // cout << "QPSolution_x" << endl
    //           << Solution_x << endl;


    // settings
    solver_y.settings()->setVerbosity(false);
    solver_y.settings()->setWarmStart(true);

    // set the initial data of the QP solver
    solver_y.data()->setNumberOfVariables(n);   //变量数n
    solver_y.data()->setNumberOfConstraints(m); //约束数m
    solver_y.data()->setHessianMatrix(hessian);
    solver_y.data()->setGradient(gradient);
    solver_y.data()->setLinearConstraintsMatrix(linearMatrix);
    solver_y.data()->setLowerBound(lowerBound_y);
    solver_y.data()->setUpperBound(upperBound_y);

    // instantiate the solver
    solver_y.initSolver();
    

    // solve the QP problem
    solver_y.solveProblem();
    Solution_y = solver_y.getSolution();
    // cout << "QPSolution_y" << endl
    //           << Solution_y << endl;

    // settings
    solver_z.settings()->setVerbosity(false);
    solver_z.settings()->setWarmStart(true);

    // set the initial data of the QP solver
    solver_z.data()->setNumberOfVariables(n);   //变量数n
    solver_z.data()->setNumberOfConstraints(m); //约束数m
    solver_z.data()->setHessianMatrix(hessian);
    solver_z.data()->setGradient(gradient);
    solver_z.data()->setLinearConstraintsMatrix(linearMatrix);
    solver_z.data()->setLowerBound(lowerBound_z);
    solver_z.data()->setUpperBound(upperBound_z);

    // instantiate the solver
    solver_z.initSolver();
    

    // solve the QP problem
    solver_z.solveProblem();
    Solution_z = solver_z.getSolution();
    // cout << "QPSolution_z" << endl
    //           << Solution_z << endl;
}

void min_snap::solve_bAp()
{
    /* x(t) = pi * t^i + ... + p5 * t^5 + p4 * t^4 + p3 * t^3 + p2 * t^2 + p1 * t + p0
    for 5th order, i = 5, size of vector to be determined p is i + 1 */
    int order = 6;
    double t0 = 0;
    double t1 = 2;
    VectorXd b = VectorXd::Zero(6); // Boundary condition: [x0, x1, v0, v1, a0, a1]
    b(1) = 2;
    b(3) = 1;
    VectorXd p(order + 1); // for order = 6, p(7)
    MatrixXd A(6, order + 1); // dimension is: (Boundary condition size: 6) * (order + 1))
    
    A.row(0) << pow(t0, 6), pow(t0, 5), pow(t0, 4), pow(t0, 3), pow(t0, 2), t0, 1;
    A.row(1) << pow(t1, 6), pow(t1, 5), pow(t1, 4), pow(t1, 3), pow(t1, 2), t1, 1;
    A.row(2) << 6 * pow(t0, 5), 5 * pow(t0, 4), 4 * pow(t0, 3), 3 * pow(t0, 2), 2 * t0, 1, 0;
    A.row(3) << 6 * pow(t1, 5), 5 * pow(t1, 4), 4 * pow(t1, 3), 3 * pow(t1, 2), 2 * t1, 1, 0;
    A.row(4) << 30 * pow(t0, 4), 20 * pow(t0, 3), 12 * pow(t0, 2), 6 * t0, 2, 0, 0;
    A.row(5) << 30 * pow(t1, 4), 20 * pow(t1, 3), 12 * pow(t1, 2), 6 * t1, 2, 0, 0;
    FullPivLU<MatrixXd> luA(A);
    int rank = luA.rank();
    cout << "Rank(A) = " << rank << endl;
    cout << "A = " << endl << A << endl;
    // cout << "A^-1 = " << endl << A.inverse() << endl;
    cout << "A^-1 = " << endl << pseudoInverse(A) << endl;
    // p = A.inverse() * b;
    p = pseudoInverse(A) * b;
    cout << "p = [" << p.transpose() << ']' << endl;
    /* The number of equations is less than the number of unknowns,
    so there are infinite solutions. 
    The solution here is the minimum norm solution. i.e. min||p|| */
    double x1 = A.row(1) * p;
    cout << "x1 = " << x1 << endl;
    double v1 = A.row(3) * p;
    cout << "v1 = " << v1 << endl;
}

void min_snap::uniform_time_arrange(MatrixXd path_)
{
    int point_num = path_.row(0).size();
    // cout << path_ << endl;
    // cout << point_num << endl;
    int seg_num = point_num - 1;
    double dist[seg_num], total_dist = 0;
    time[0] = 0;
    for(int i = 0; i < seg_num; i++)
    {
        dist[i] = (path_.col(i + 1) - path_.col(i)).norm();
        total_dist += dist[i];
        // cout << "dist: " << dist[i] << endl;
    }
    for(int i = 0; i < point_num; i++)
    {
        time[i + 1] = time[i] + dist[i] / total_dist * total_time;
        
    }
    // for(int i = 0; i < point_num; i++)
    // {
    //     cout << "time: " << time[i] << endl;
    // }
}

int min_snap::solveQP(MatrixXd Hessian)
{
    // allocate QP problem matrices and vectores
    
    return 0;
}

void min_snap::visTraj(VectorXd Sol_x, VectorXd Sol_y, VectorXd Sol_z)
{
    // 定义可视化消息
    visualization_msgs::Marker marker;
    visualization_msgs::MarkerArray marker_array;
    marker.header.frame_id = "map";
    marker.header.stamp = ros::Time::now();
    marker.ns = "trajectory";
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.orientation.w = 1.0;
    marker.id = 0;
    marker.type = visualization_msgs::Marker::SPHERE;
    marker.scale.x = 0.1;
    marker.color.r = 1.0;
    marker.color.a = 1.0;

    nav_msgs::Path traj;
    traj.header.frame_id = "/map";
    traj.header.stamp = ros::Time::now();


    for (size_t i = 0; i < seg_num; i++)
    {
        // 定义多项式系数
        VectorXd coeffsX(p_num), coeffsY(p_num), coeffsZ(p_num);
        for (size_t j = 0; j < p_num; j++)
        {
            coeffsX(j) = Sol_x(i * p_num + j);
            coeffsY(j) = Sol_y(i * p_num + j);
            coeffsZ(j) = Sol_z(i * p_num + j);
        }
        // cout << coeffsX.transpose() << ", " << coeffsY.transpose() << ", " << coeffsZ.transpose() << endl;
        // 定义起始时间和结束时间
        double startTime, endTime, timeStep;
        startTime = time[i];
        endTime = time[i + 1];
        timeStep = 0.05;
        // cout << "startTime = " << startTime << ", endTime = " << endTime << endl;
        for (double t = startTime; t <= endTime; t += timeStep)
        {
            Vector3d v = CalcVector(coeffsX, coeffsY, coeffsZ, t);
            geometry_msgs::PoseStamped p;
            p.header.stamp = ros::Time::now();
            p.pose.position.x = v.x();
            p.pose.position.y = v.y();
            p.pose.position.z = v.z();
            // p.x = v.x();
            // p.y = v.y();
            // p.z = v.z();
            // cout << v.transpose() << endl;
            traj.poses.push_back(p);
        }
    }
    traj_pub.publish(traj);
    ros::spin();
}

void min_snap::waypoint_cb(const geometry_msgs::Twist::ConstPtr & msg)
{
    // TODO
}

// 计算多项式函数值的函数
double min_snap::PolyEval(const VectorXd& coeffs, double x) {
    
    double result = 0.0;
    // for (int i = coeffs.size() - 1; i >= 0; --i) {
    //     result = result * x + coeffs[i];
    // }

    for (int i = 0; i < coeffs.size(); i++) {
        result += pow(x, i) * coeffs(i);
    }
    return result;
}

// 计算三维向量的函数
Vector3d min_snap::CalcVector(const VectorXd& coeffsX, const VectorXd& coeffsY, const VectorXd& coeffsZ, double t) {
    Vector3d v;
    v.x() = PolyEval(coeffsX, t);
    v.y() = PolyEval(coeffsY, t);
    v.z() = PolyEval(coeffsZ, t);
    return v;
}