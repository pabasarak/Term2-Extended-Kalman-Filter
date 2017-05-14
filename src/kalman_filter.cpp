#include "kalman_filter.h"
#include "tools.h"
#include <stdlib.h>
#include <iostream>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

# define PI           3.14159265358979323846
KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  
  	x_ = F_ * x_;
	//cout << "EKF: x_" << endl;
	MatrixXd Ft = F_.transpose();
	//cout << "EKF: Ft" << endl;
	P_ = F_ * P_ * Ft + Q_;
	//cout << "EKF: P_" << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
   MatrixXd hx =VectorXd(3);
  if(((px * px) < 0.00001) || ((py * py) < 0.00001))
	{
	cout << "Error while converting to polar coordinates: Division by Zero" << endl;
	hx <<0,
		0,
		0;
	//return;
	}else{
 
	//cout << "px " << px << "  py "<< py << endl;
	double rho = sqrt(px * px + py * py);
	const double phi = atan2(py, px);
	const double rho_dot = (px * vx + py * vy) / (rho);
  hx <<sqrt(px*px + py*py),
			atan2(py,px),
			(px*vx + py*vy)/(sqrt(px*px + py*py));
			
	hx << rho, phi, rho_dot;
	//cout << hx << endl;
	}
	 Tools tools = Tools();
	MatrixXd Hj_ = tools.CalculateJacobian(x_);
	VectorXd y = z - hx;
	// normalizing y(1)
	while (y(1) < -PI)
		y(1) += 2 * PI;
		while (y(1) > PI)
		y(1) -= 2 * PI;
	MatrixXd Ht = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj_) * P_;
	
}
