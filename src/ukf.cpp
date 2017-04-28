#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  //P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5, 5);      // initial covariance is identity matrix

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8;     // was ridiculous value of 30...

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;   // was ridiculous value of 30...

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  /*
    Missing initialisations are:
      n_x_
      n_aug_
      lambda_
      weights_
  */
  n_x_ = 5;       //set state dimension
  n_aug_ = 7;     //set augmented dimension
  lambda_ = 3 - n_aug_;  //define spreading parameter

  //set vector for weights - this depends only on n_aug_ and lambda_, so can calculate once
  VectorXd weights_ = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*
  The steps are as follows:
    1. Check initialised
    2. Predict
    3. Update
    Note that the process needs to be aware of differences between lidar/radar
  */

  std::cout << "Processing..." << std::endl;  // test code
  std::cout << "Measurement: " << std::endl << meas_package.raw_measurements_ << std::endl;
  std::cout << "Sensor: " << meas_package.sensor_type_ << std::endl;

  // Check the initialisation
  if (is_initialized_ != true) {
    // need to initialise the x_ state vector
    // depending on whether the measurement is radar or lidar (radar needs polar co-ord transform)

    std::cout << "Initialising..." << std::endl;  // test code
    //x_.fill(0.0);    // test code

    // TO-DO: what else do we need in here???
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // This is a lidar measurement - initialise as such, using px,py from the measurement
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(0), 0.0, 0.0, 0.0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // This is a radar measurement - initialise as such, converting to state space
      const double rho_init = meas_package.raw_measurements_(0);
      const double phi_init = meas_package.raw_measurements_(1);
      const double px_init = rho_init * std::cos(phi_init);
      const double py_init = rho_init * std::sin(phi_init);
      x_ << px_init, py_init, 0.0, 0.0, 0.0;
    }

    // concerned about an empty x_ matrix - check for this and don't initialise if so...
    if (std::fabs(x_.norm()) < dblZeroThreshold) {
      // norm is below zero threshold
      return;       // early return if this is the case - and not setting the initialisation flag
    }

    previous_timestamp_ = meas_package.timestamp_;
    
    std::cout << "Timestamp: " << previous_timestamp_ << std::endl;   // test code

    is_initialized_ = true;

    return;         // early return if only initialising the object

  }

  // Predict (motion model)
  const long long new_timestamp = meas_package.timestamp_;      // get the measurement timestamp
  const double delta_t = (new_timestamp - previous_timestamp_) * dblMicroSectoSec;

  std::cout << "Time delta: " << delta_t << std::endl;          // test code

  Prediction(delta_t);                                          // Make the prediction

  previous_timestamp_ = new_timestamp;

  // Update (sensor models)
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    // this is a laser measurement and we are processing those
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // this is a radar measurement and we are processing those
    UpdateRadar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  std::cout << "Prediction called, delta_t = " << delta_t << std::endl;     // test code

  /*
  NM notes: the prediction process uses a single motion model (inc. velocity estimate)
  Had to stare at this for a while to convince myself we can do this even when not using radar measurements
  (the Kalman Filter estimates velocity for us...)
  */

  // First create Sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  GenerateSigmaPoints(&Xsig_pred);
  
  std::cout << "Created Sigma Points = " << Xsig_pred << std::endl;     // test code

  // What next??

  return;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  std::cout << "LIDAR update called, raw data = " << meas_package.raw_measurements_ << std::endl;     // test code

  return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  std::cout << "RADAR update called, raw data = " << meas_package.raw_measurements_ << std::endl;     // test code
  return;

  //set measurement dimension, radar can measure r, phi, and r_dot - could do this once only for efficiency
  int n_z = 3;

  // get predicted sigma points in state space  Xsig_pred
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  // TO-DO

  // get predicted state mean                   x_  - this is known

  // get predicted state covariance             P_  - this is known

  // get sigma points in measurement space      Zsig
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

                                             // extract values for better readibility
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);
    double v = Xsig_pred(2, i);
    double yaw = Xsig_pred(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1, i) = atan2(p_y, p_x);                                 //phi
    Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  // get vectore of mean predicted measurement  z_pred
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // get predicted measurement covariance       S
  MatrixXd S = MatrixXd(n_z, n_z);
  //MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
                                             //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix       - ** can calc R just once?? **
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;
  S = S + R;

  // get incoming radar measurement             z
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;

  // create matrix for cross correlation        Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

                                             //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}

/*
  Sigma point generation function
*/
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  // state dimension, lamba spreading param, state vector x_ and state covariance P_ are all set already

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0) = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsig.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }

  //write result
  *Xsig_out = Xsig;
}