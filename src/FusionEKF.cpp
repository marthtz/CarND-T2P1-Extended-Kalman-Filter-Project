#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

#define PI_X2 (2*M_PI)

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


// Set acceleration noise
static double noise_ax = 9;
static double noise_ay = 9;


/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  // LIDAR Measurement covariance matrix
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  // RADAR measurement covariance matrix
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    // Initialize in case there's no RADAR or LIDAR measurement 
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Get all 3 RADAR measurements
      double rho_measured = measurement_pack.raw_measurements_[0];
      double phi_measured = measurement_pack.raw_measurements_[1];
      double rhodot_measured = measurement_pack.raw_measurements_[2];

      // Convert from polar to cartesian coordinates
      double x = rho_measured * cos(phi_measured);
      double y = rho_measured * sin(phi_measured);

      // According to tips and tricks, there's not enough information to set
      // vx and vy with a RADAR initialization,. Therefore, set to 0.
      double vx = 0;// rhodot_measured * cos(phi_measured);
      double vy = 0;// rhodot_measured * sin(phi_measured);

      // Set state vector with first measurement
      ekf_.x_ << x,  y, vx, vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity

      // Get both LIDAR measurements - set vx and vy to 0
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // Set first timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // Print init values
    //cout << ekf_.x_ << end
    //cout << previous_timestamp_ << endl;

    // Init state covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    // Init state transition matrix F
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // Compute the time elapsed between the current and previous measurements
  // dt is expressed in seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  // Integrate dt in transition matrix F
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Set process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  // Start a new prediction
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    Tools tools;
    Eigen::VectorXd current_measurement;
    current_measurement = VectorXd(3);

    // Get current measurment data for RADAR
    //double rho_measured = measurement_pack.raw_measurements_[0];
    double phi_measured = measurement_pack.raw_measurements_[1];
    //double rhodot_measured = measurement_pack.raw_measurements_[2];

    // Normalize phi as measurement data show values larger and smaller than PI
    if (phi_measured > M_PI)
    {
      phi_measured -= PI_X2;
    }
    else if (phi_measured < -M_PI)
    {
      phi_measured += PI_X2;
    }

    // Set measurement vector with current values
    current_measurement << measurement_pack.raw_measurements_[0],
                           phi_measured,
                           measurement_pack.raw_measurements_[2];

    // Calculate Jacobian matrix of current state vector
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    // Set to use RADAR covariance matrix
    ekf_.R_ = R_radar_;

    // Start update with current measurement
    ekf_.UpdateEKF(current_measurement);
  }
  else
  {
    // Laser updates
    // Set to use LIDAR H and R matrices
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    // Start update with current measurement
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}
