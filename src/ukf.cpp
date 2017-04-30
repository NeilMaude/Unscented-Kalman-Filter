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
      Xsig_pred_
  */
  n_x_ = 5;       //set state dimension
  n_aug_ = 7;     //set augmented dimension
  lambda_ = 3 - n_aug_;  //define spreading parameter

  //set vector for weights - this depends only on n_aug_ and lambda_, so can calculate just the once
  weights_ = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  n_z_laser_ = 2;   // set lidar measurement space dimension
  n_z_radar_ = 3;   // set radar measurement space dimension

  // Can calc measurement noise covariance matrices just once in constructor
  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  DONE:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*
  The steps are as follows:
    1. Check initialised
    2. Predict (generate Sigma points, predict Sigma points, predict mean and covariance)
    3. Update (predict measurement, update state)
    Note that the process needs to be aware of differences between lidar/radar
  */

  /*
  // test code
  std::cout << "Processing..." << std::endl;  
  std::cout << "Measurement: " << std::endl << meas_package.raw_measurements_ << std::endl;
  std::cout << "Sensor: " << meas_package.sensor_type_ << std::endl;
  */

  // Check the initialisation
  if (is_initialized_ != true) {
    // need to initialise the x_ state vector
    // depending on whether the measurement is radar or lidar (radar needs polar co-ord transform)

    // test code
    //std::cout << "Initialising..." << std::endl;  
    
    // Initialise initial state from the first measurement, depending on the sensor type
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // This is a lidar measurement - initialise as such, using px,py from the measurement
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(0), 0.0, 0.0, 0.0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // This is a radar measurement - initialise as such, converting to state space
      // Note: there is velocity data from the radar measurement, but not using this here
      const double rho_init = meas_package.raw_measurements_(0);
      const double phi_init = meas_package.raw_measurements_(1);
      const double px_init = rho_init * std::cos(phi_init);
      const double py_init = rho_init * std::sin(phi_init);
      x_ << px_init, py_init, 0.0, 0.0, 0.0;
    }

    // Concerned about an empty x_ matrix - check for this
    // Don't initialise if so ... will get it next time if that is the case
    if (std::fabs(x_.norm()) < dblZeroThreshold) {
      // norm is below zero threshold
      return;       // early return if this is the case - and not setting the initialisation flag
    }

    // Save the timestamp for the next iteration
    previous_timestamp_ = meas_package.timestamp_;
    
    // test code
    //std::cout << "Timestamp: " << previous_timestamp_ << std::endl;   

    is_initialized_ = true;

    return;         // early return from ProcessMeasurement, if only initialising the object

  }

  // Predict (motion model)
  const long long new_timestamp = meas_package.timestamp_;      // get the measurement timestamp
  const double delta_t = (new_timestamp - previous_timestamp_) * dblMicroSectoSec;

  // test code
  //std::cout << "Time delta: " << delta_t << std::endl;          

  Prediction(delta_t);                                          // Make the prediction

  previous_timestamp_ = new_timestamp;                          // Save the timestamp for next iteration

  // Update, using correct sensor model, if enabled
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
  DONE:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix (P_).
  */
  
  // test code
  //std::cout << "Prediction called, delta_t = " << delta_t << std::endl;     

  /*
  NM notes: the prediction process uses a single motion model (inc. velocity estimate)
  Had to stare at this for a while to convince myself we can do this the same way on each iteration,
  regardless of the sensor type (Kalman Filters don't use sensor data until the update step).
  */

  // Create Augmented Sigma Points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  GenerateAugmentedSigmaPoints(&Xsig_aug);

  // test code
  //std::cout << "Created Augmented Sigma Points = " << std::endl << Xsig_aug << std::endl;     

  // Predict these Sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  SigmaPointPrediction(Xsig_aug, delta_t, &Xsig_pred);

  // test code
  //std::cout << "Predicted Sigma Points = " << std::endl << Xsig_pred << std::endl;     

  // Predict mean and covariance and update x_, P_
  PredictMeanAndCovariance(Xsig_pred);

  // Save the predicted Sigma points (Xsig_pred) to use in the update process
  // This will save generating more Sigma points for the measurement prediction
  Xsig_pred_ = Xsig_pred;

  return;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  DONE:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // test code
  //std::cout << "LIDAR update called, raw data = " << meas_package.raw_measurements_ << std::endl;     

  // First predict the measurement

  // get predicted sigma points in state space  Xsig_pred_ - this is retained from the prediction step
  // get predicted state mean                   x_  - this is known
  // get predicted state covariance             P_  - this is known

  // get sigma points in measurement space, Zsig, by transforming predicted Sigma points saved in predict step
  // Measurement space dimension is 2 (px, py), taken from n_z_laser_
  MatrixXd Zsig = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //sigma point predictions in process space
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    //sigma point predictions in measurement space
    Zsig(0, i) = px;
    Zsig(1, i) = py;
  }

  // calc mean predicted measurement vector z_pred
  VectorXd z_pred = VectorXd(n_z_laser_);
  z_pred.fill(0.0);
  //mean predicted measurement
  z_pred = Zsig * weights_;

  // calc measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_laser_, n_z_laser_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 Sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  S = S + R_laser_;

  // Have now found the mean prediction and covariance

  // Now do the update process

  // get the incoming measurements
  VectorXd z = VectorXd(n_z_laser_);
  z = meas_package.raw_measurements_;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_laser_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  
  //print result - test code
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

  //NIS Lidar Update
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

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

  // test code
  //std::cout << "RADAR update called, raw data = " << meas_package.raw_measurements_ << std::endl;     
  
  //set measurement dimension, radar can measure r, phi, and r_dot - could do this once only for efficiency
  //int n_z = 3;

  // First predict the measurement

  // get predicted sigma points in state space  Xsig_pred_ - this is retained from the prediction step
  // get predicted state mean                   x_  - this is known
  // get predicted state covariance             P_  - this is known

  // Get sigma points in measurement space, Zsig, by transforming previous predicted Sigma points into measurement space
  // Measurement space dimension is 3, as held in n_z_radar_
  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    // extract values from predictions - will be in the state space
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;
    // save into measurement space
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                             //r
    Zsig(1, i) = atan2(p_y, p_x);                                         //phi
    Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);     //r_dot
  }

  // calc mean predicted measurement vector z_pred
  VectorXd z_pred = VectorXd(n_z_radar_);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // get predicted measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    // calculate S
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  S = S + R_radar_;

  // Have now got the mean prediction and covariance
  
  // Now do the update process

  // get incoming radar measurement, z
  VectorXd z = VectorXd(n_z_radar_);
  z << meas_package.raw_measurements_;

  // create matrix for cross correlation, Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {          //2n+1 Sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
    // Calculate Tc
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

  //print result - test code
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

  //NIS Update
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/*
Augmented Sigma point generation function
*/
void UKF::GenerateAugmentedSigmaPoints(MatrixXd* Xsig_out) {
  // state dimension, lamba spreading param, state vector x_ and state covariance P_ are all set already

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;             // why are these zero always?
  x_aug(6) = 0;
  //std::cout << "x_aug:" << std::endl << x_aug << std::endl;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_*std_a_;
  P_aug(6, 6) = std_yawdd_*std_yawdd_;
  //std::cout << "P_aug:" << std::endl << P_aug << std::endl;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;
}

/*
Sigma point prediction, using augmented Sigma points input, outputting non-augmented prediction
*/
void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t, MatrixXd* Xsig_out) {

  //create matrix with predicted sigma points as columns - will be 5 * 15
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // loop over the columns (each column represents a sigma point)
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //std::cout << "Column: " << i << std::endl << Xsig_aug.col(i) << std::endl;

    // retrieve the values
    double px = Xsig_aug.col(i)(0);
    //std::cout << "px = " << px << std::endl;
    double py = Xsig_aug.col(i)(1);
    double vel = Xsig_aug.col(i)(2);
    double rho = Xsig_aug.col(i)(3);
    double rhodot = Xsig_aug.col(i)(4);
    double vel_err = Xsig_aug.col(i)(5);
    double rhovel_err = Xsig_aug.col(i)(6);

    double new_px = 0;
    double new_py = 0;
    double new_vel = 0;
    double new_rho = 0;
    double new_rhodot = 0;

    // some repeated values used in the calcs - do these just once per call to this function
    double rho_dt = rho + (delta_t * rhodot);
    double rho_sin = sin(rho);
    double rho_cos = cos(rho);
    double vk_over_rho = vel / rhodot;

    if (rhodot >= dblZeroThreshold) {
      // this is the case for non-zero rhovel (i.e. there is a rate of change of bearing/not straight line motion)
      new_px = px + (vk_over_rho * (sin(rho_dt) - rho_sin)) + (0.5 * delta_t * delta_t * rho_cos * vel_err);
      new_py = py + (vk_over_rho * (rho_cos - cos(rho_dt))) + (0.5 * delta_t * delta_t * rho_sin * vel_err);
      new_vel = vel + (delta_t * vel_err);
      new_rho = rho + (delta_t * rhodot) + (0.5 * delta_t * delta_t *  rhovel_err);
      new_rhodot = rhodot + (delta_t * rhovel_err);

      //  std::cout << "Non-zero rhodot, new values: " << new_px << " : " <<
      //    new_py << " : " << new_vel << " : " << new_rho << " : " << new_rhodot << std::endl;
    }
    // check for zeros
    if (rhodot < dblZeroThreshold) {
      // this is the case for zero rhovel (zero rate of bearing change)
      new_px = px + (vel * rho_cos * delta_t) + (0.5 * delta_t * delta_t * rho_cos * vel_err);
      new_py = py + (vel * rho_sin * delta_t) + (0.5 * delta_t * delta_t * rho_sin * vel_err);
      new_vel = vel + (delta_t * vel_err);
      new_rho = rho + (delta_t * rhodot) + (0.5 * delta_t * delta_t * rhovel_err);
      new_rhodot = rhodot + (delta_t * rhovel_err);

      //  std::cout << "Zero rhodot, new values: " << new_px << " : " <<
      //    new_py << " : " << new_vel << " : " << new_rho << " : " << new_rhodot << std::endl;
    }
    VectorXd Xpred = VectorXd(n_x_);
    Xpred << new_px, new_py, new_vel, new_rho, new_rhodot;
    Xsig_pred.col(i) = Xpred;
  }

  //print result
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;     // test code

  //write result
  *Xsig_out = Xsig_pred;

}

/*
Predict mean and covariance from a set of predicted Sigma points
*/
void UKF::PredictMeanAndCovariance(MatrixXd Xsig_pred) {

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  //std::cout << "weights_ = " << weights_ << std::endl;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    
    /*std::cout << "Debug: " << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "i = " << i << std::endl;
    std::cout << "weights_(i) = " << weights_(i) << std::endl;
    std::cout << "Xsig_pred.col(i) = " << Xsig_pred.col(i) << std::endl;
    */

    x = x + weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose();

  }

  /*
  //print result - test code
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;
  */ 

  //store result into state vector x_ and state covariance matrix P_
  x_ = x;
  P_ = P;
}