#include "FusionEKF.h"
#include "Eigen/Dense"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  // measurement matrix - laser
  H_laser_ << 1, 0, 0, 0, 0, 1, 0, 0;
  // measurement matrix - radar
  // Hj_;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    // ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates
      //         and initialize state.
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_rate = measurement_pack.raw_measurements_[2];
      float position_x = rho * cos(phi);
      float position_y = rho * sin(phi);
      float velocity_x = rho_rate * cos(phi);
      float velocity_y = rho_rate * sin(phi);
      ekf_.x_ << position_x, position_y, velocity_x, velocity_y;
      H_in = tools.CalculateJacobian(const VectorXd &ekf_.x_);
      ekf_.R_ = R_radar_;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1], 0, 0;
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
    }

    // check if x, y are zeros
    if (ekf_.x_[0] < 0.0001) {
      ekf_.x_[0] = 0.0001;
    }
    if (ekf_.x_[1] < 0.0001) {
      ekf_.x_[1] = 0.0001;
    }
    // initialize object covariance matrix P
    MatrixXd P_in = MatrixXd(4, 4);
    P_in << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1000.0, 0.0, 0.0,
        0.0, 0.0, 1000.0;
    ekf_.P_ = P_in;
    // initialize process covariance matrix Q
    MatrixXd Q_in = MatrixXd(4, 4);
    Q_in << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0;
    ekf_.Q_ = Q_in;
    // initialize time
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed
   * time. Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // calculate the time elapse
  double t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  // the initial transition matrix F_
  F_in = MatrixXd(4, 4);
  F_in << 1.0, 0.0, t, 0.0, 0.0, 1.0, 0.0, t, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0;
  ekf_.F_ = F_in;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    kf_.R_ = R_radar_;
    kf_.H_ = tools.CalculateJacobian(const VectorXd &ekf_.x_);
    ;
    kf_.UpdateEKF(const VectorXd &z);
  } else {
    // TODO: Laser updates
    kf_.R_ = R_laser_;
    kf_.H_ = H_laser_;
    kf_.Update(const VectorXd &z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
