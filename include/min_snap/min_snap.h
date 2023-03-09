#include <ros/ros.h>
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Core> 
#include <Eigen/SVD>
#include "OsqpEigen/OsqpEigen.h"
#include <stdexcept>


template<typename _Matrix_Type_> 
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = 
    std::numeric_limits<double>::epsilon()) 
{  
    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);  
    double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);  
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
    ros::Publisher traj_pub;
    static const int wp_num = 5;
    double time[wp_num];
    double max_velocity = 10;
    double total_time = 5;
public:
    
    Eigen::MatrixXd path = Eigen::MatrixXd::Zero(dim, wp_num);
    Eigen::VectorXd Solution_x, Solution_y, Solution_z;
    // path.push_back(wp1);
    double v0[2] = {0, 0};
    double a0[2] = {0, 0};
    double v1[2] = {0, 0};
    double a1[2] = {0, 0};
    double T = 5;
    
    min_snap(/* args */);
    ~min_snap();
    void solve_bAp();
    void solve_Nseg_bAp();
    void uniform_time_arrange(Eigen::MatrixXd path_);
    int solveQP(Eigen::MatrixXd Hessian);
    void getTraj(Eigen::VectorXd Sol_x, Eigen::VectorXd Sol_y, Eigen::VectorXd Sol_z);
    void waypoint_cb(const geometry_msgs::Twist::ConstPtr & msg);
};

min_snap::min_snap(ros::NodeHandle nh)
{
    waypt_path_sub = nh.subscribe("waypoint", 1000, &min_snap::waypoint_cb, this);
    // path.col(0) = Eigen::Vector2d(0, 0);
    // path.col(1) = Eigen::Vector2d(1, 2);
    // path.col(2) = Eigen::Vector2d(2, -1);
    // path.col(3) = Eigen::Vector2d(4, 8);
    // path.col(4) = Eigen::Vector2d(5, 2);
    // uniform_time_arrange(path);
    // solve_Nseg_bAp();

    // solve_bAp();
    // std::cout << " path: " << std::endl << path << std::endl;
    // path.push_back(wp1);
    // path.push_back(wp2);
    // path.push_back(wp3);
    // path.push_back(wp4);
    // path.push_back(wp5);
    // std::cout << " path: " << path.size() << std::endl;
    // time_arrange(path);
}

min_snap::~min_snap()
{
}

void min_snap::solve_Nseg_bAp()
{
    // p(t) = p0 + p1 * t + p2 * t2 + ... + pn * tn = ∑ pi * ti
    int derivative_order = 3; // pos = 0, vel = 1, acc = 2, jerk = 3, snap = 4;
    int seg_num = path.row(0).size() - 1;
    int p_order = 2 * derivative_order - 1; // Polynomial order, for jerk is 5, for snap is 7
    int p_num = p_order + 1;
    int m = 4 * wp_num - 2;
    int n = p_num * seg_num;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(p_num * seg_num, p_num * seg_num);
    Eigen::MatrixXd Qi = Eigen::MatrixXd::Zero(p_num, p_num);
    Eigen::VectorXd t = Eigen::VectorXd::Ones(p_num);

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

    Eigen::SparseMatrix<double> hessian(n, n);      //P: n*n正定矩阵,必须为稀疏矩阵SparseMatrix
    hessian = Q.sparseView();
    // std::cout << hessian << std::endl;
    
    // Eigen::LLT<Eigen::MatrixXd> lltOfA(Q); // compute the Cholesky decomposition of A 
    // std::cout << lltOfA.info() << std::endl;
    // if(lltOfA.info() == Eigen::NumericalIssue) 
    // { 
    //     throw std::runtime_error("Possibly non semi-positive definitie matrix!"); 
    // }

    
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(n);                    //Q: n*1向量
    Eigen::SparseMatrix<double> linearMatrix(m, n); //A: m*n矩阵,必须为稀疏矩阵SparseMatrix
    Eigen::VectorXd lowerBound_x(m), lowerBound_y(m), lowerBound_z(m);                  //L: m*1下限向量
    Eigen::VectorXd upperBound_x(m), upperBound_y(m), upperBound_z(m);                  //U: m*1上限向量
    
    Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(m, n);
    Eigen::MatrixXd Ai = Eigen::MatrixXd::Zero(1, p_num);
    Eigen::MatrixXd Ai_ = Eigen::MatrixXd::Zero(1, p_num);
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
        
        // lowerBound_z[2 * i] = path.col(i).z();
        // upperBound_z[2 * i] = path.col(i).z();
        // lowerBound_z[2 * i + 1] = path.col(i + 1).z();
        // upperBound_z[2 * i + 1] = path.col(i + 1).z();
        
    }
    Eigen::MatrixXd Veli = Eigen::MatrixXd::Zero(1, 2 * p_num);
    Eigen::MatrixXd Acci = Eigen::MatrixXd::Zero(1, 2 * p_num);
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
            // Veli.block(0, p_num, 1, p_num) = Eigen::MatrixXd::Zero(1, p_num);
            // Acci.block(0, p_num, 1, p_num) = Eigen::MatrixXd::Zero(1, p_num);
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
            // lowerBound_z[2 * seg_num] = 0;
            // upperBound_z[2 * seg_num] = 0;
            // lowerBound_z[2 * seg_num + 1] = 0;
            // upperBound_z[2 * seg_num + 1] = 0;
        }
        else if(i == wp_num - 1)
        {
            // Veli.block(0, 0, 1, p_num) = Eigen::MatrixXd::Zero(1, p_num);
            // Acci.block(0, 0, 1, p_num) = Eigen::MatrixXd::Zero(1, p_num);
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
            // lowerBound_z[m - 2] = 0;
            // upperBound_z[m - 2] = 0;
            // lowerBound_z[m - 1] = 0;
            // upperBound_z[m - 1] = 0;
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
            // lowerBound_z[2 * i + 2 * seg_num] = 0;
            // upperBound_z[2 * i + 2 * seg_num] = 0;
            // lowerBound_z[2 * i + 2 * seg_num + 1] = 0;
            // upperBound_z[2 * i + 2 * seg_num + 1] = 0;
        }
    }
    // std::cout << Aeq << std::endl;
    linearMatrix = Aeq.sparseView();
    // std::cout << linearMatrix << std::endl;
    // std::cout << lowerBound_x << std::endl;

    // instantiate the solver
    OsqpEigen::Solver solver;
    // settings
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // set the initial data of the QP solver
    solver.data()->setNumberOfVariables(n);   //变量数n
    solver.data()->setNumberOfConstraints(m); //约束数m
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(lowerBound_x);
    solver.data()->setUpperBound(upperBound_x);

    // instantiate the solver
    solver.initSolver();
    

    // solve the QP problem
    solver.solveProblem();
    Solution_x = solver.getSolution();
    std::cout << "QPSolution_x" << std::endl
              << Solution_x << std::endl;

    solver.data()->setLowerBound(lowerBound_y);
    solver.data()->setUpperBound(upperBound_y);
    solver.solveProblem();
    Solution_y = solver.getSolution();
    std::cout << "QPSolution_y" << std::endl
              << Solution_y << std::endl;

    // solver.data()->setLowerBound(lowerBound_z);
    // solver.data()->setUpperBound(upperBound_z);
    // solver.solveProblem();
    // Solution_y = solver.getSolution();
    // std::cout << "QPSolution_z" << std::endl
    //           << Solution_z << std::endl;
}

void min_snap::solve_bAp()
{
    /* x(t) = pi * t^i + ... + p5 * t^5 + p4 * t^4 + p3 * t^3 + p2 * t^2 + p1 * t + p0
    for 5th order, i = 5, size of vector to be determined p is i + 1 */
    int order = 6;
    double t0 = 0;
    double t1 = 2;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(6); // Boundary condition: [x0, x1, v0, v1, a0, a1]
    b(1) = 2;
    b(3) = 1;
    Eigen::VectorXd p(order + 1); // for order = 6, p(7)
    Eigen::MatrixXd A(6, order + 1); // dimension is: (Boundary condition size: 6) * (order + 1))
    
    A.row(0) << pow(t0, 6), pow(t0, 5), pow(t0, 4), pow(t0, 3), pow(t0, 2), t0, 1;
    A.row(1) << pow(t1, 6), pow(t1, 5), pow(t1, 4), pow(t1, 3), pow(t1, 2), t1, 1;
    A.row(2) << 6 * pow(t0, 5), 5 * pow(t0, 4), 4 * pow(t0, 3), 3 * pow(t0, 2), 2 * t0, 1, 0;
    A.row(3) << 6 * pow(t1, 5), 5 * pow(t1, 4), 4 * pow(t1, 3), 3 * pow(t1, 2), 2 * t1, 1, 0;
    A.row(4) << 30 * pow(t0, 4), 20 * pow(t0, 3), 12 * pow(t0, 2), 6 * t0, 2, 0, 0;
    A.row(5) << 30 * pow(t1, 4), 20 * pow(t1, 3), 12 * pow(t1, 2), 6 * t1, 2, 0, 0;
    Eigen::FullPivLU<Eigen::MatrixXd> luA(A);
    int rank = luA.rank();
    std::cout << "Rank(A) = " << rank << std::endl;
    std::cout << "A = " << std::endl << A << std::endl;
    // std::cout << "A^-1 = " << std::endl << A.inverse() << std::endl;
    std::cout << "A^-1 = " << std::endl << pseudoInverse(A) << std::endl;
    // p = A.inverse() * b;
    p = pseudoInverse(A) * b;
    std::cout << "p = [" << p.transpose() << ']' << std::endl;
    /* The number of equations is less than the number of unknowns,
    so there are infinite solutions. 
    The solution here is the minimum norm solution. i.e. min||p|| */
    double x1 = A.row(1) * p;
    std::cout << "x1 = " << x1 << std::endl;
    double v1 = A.row(3) * p;
    std::cout << "v1 = " << v1 << std::endl;
}

void min_snap::uniform_time_arrange(Eigen::MatrixXd path_)
{
    int point_num = path_.row(0).size();
    // std::cout << path_ << std::endl;
    // std::cout << point_num << std::endl;
    int seg_num = point_num - 1;
    double dist[seg_num], total_dist = 0;
    time[0] = 0;
    for(int i = 0; i < seg_num; i++)
    {
        dist[i] = (path_.col(i + 1) - path_.col(i)).norm();
        total_dist += dist[i];
        // std::cout << "dist: " << dist[i] << std::endl;
    }
    for(int i = 0; i < point_num; i++)
    {
        time[i + 1] = time[i] + dist[i] / total_dist * total_time;
        
    }
    // for(int i = 0; i < point_num; i++)
    // {
    //     std::cout << "time: " << time[i] << std::endl;
    // }
}

int min_snap::solveQP(Eigen::MatrixXd Hessian)
{
    // allocate QP problem matrices and vectores
    
    return 0;
}

void min_snap::getTraj(Eigen::VectorXd Sol_x, Eigen::VectorXd Sol_y, Eigen::VectorXd Sol_z)
{
    Eigen::Vector3d 
}