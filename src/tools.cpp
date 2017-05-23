#include <iostream>
#include "tools.h"

#define DIV0_LIMIT_CHECK (0.0000001)

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // Check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if((estimations.size() != ground_truth.size()) || (estimations.size() == 0))
  {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // Accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i)
  {
    VectorXd residual = estimations[i] - ground_truth[i];

    // Coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
}

  // Calculate the mean
  rmse = rmse/estimations.size();

  // Calculate the squared root
  rmse = rmse.array().sqrt();

  // Return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);

  // Recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // Calculate some terms
  double c1 = (px * px) + (py * py);
  // Check division by zero
  if(fabs(c1) < DIV0_LIMIT_CHECK)
  {
    c1 = DIV0_LIMIT_CHECK;
  }
  double c2 = sqrt(c1);
  double c3 = (c1 * c2);

  // Compute the Jacobian matrix
  Hj <<  (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
         py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
