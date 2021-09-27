#pragma once

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <chrono>
#include <string.h>
#include <iostream>

#include "imu_ekf_ros/initRequest.h"
#include "sensor_msgs/Imu.h"
#include "std_msgs/Int32MultiArray.h"
#include "geometry_msgs/Vector3.h"
#include "geometry_msgs/Quaternion.h"
#include <tf/transform_broadcaster.h>
#include <cstdlib>

#include "ros/ros.h"

namespace EKF
{
	void debug(auto str)
	{
		std::cout << str << std::endl;
	}
}

// timer class to time parts of
class Timer 
{
public:
	Timer()
	{
		last = std::chrono::steady_clock::now();
	}
	float PrintDt(std::string label)
	{
		std::cout << label << ": " << CalculateTimeDiff() << std::endl;
	}

private:
	std::chrono::steady_clock::time_point last;
	float CalculateTimeDiff()
	{
		const auto old = last;
		last = std::chrono::steady_clock::now();
		const  std::chrono::duration<float> frame_time = last-old;
		return frame_time.count();
	}
};

// struct which contains constant parameters used in algorithm
struct EKF_struct {
		// global variables in filter
		Eigen::Matrix<double,12,12> Q = Eigen::Matrix<double,12,12>::Zero(); // noise matrix for IMU
		Eigen::Matrix<double,3,3> Ra = Eigen::Matrix<double,3,3>::Zero(); // noise matrix for accelerometer
		double g; // gravity
		double dt; // time step
		double num_data; // number of data collected for initialization
		Timer timer;
};

class IMU_EKF
{
public:

	IMU_EKF(ros::NodeHandle & n);
	void initialize_ekf(ros::NodeHandle &n);

	Eigen::Matrix<double,10, 1> get_state() const
	{
		return m_state;
	}

private:

	// State vector [orientation, gyro bias, accel bias]
	Eigen::Matrix<double,10, 1> m_state; // state
	Eigen::Matrix<double, 9, 9> m_cov; // covariance
	Eigen::Matrix<double,3,1> m_g_pred; // predicted gravity
	Eigen::Matrix<double,3,1> m_g_pred_sum; // sum of all predicted gravities
	Eigen::Matrix<double,3,1> m_g_true; // true gravity vector
	Eigen::Matrix<double,3,3> Re; // measurement noise matrix
	EKF_struct m_filter; // filter object

	// ROS related subscribers
	ros::Subscriber m_imu_subscriber;
	ros::Subscriber m_sun_sensor_subscriber;
	ros::Publisher m_orientation_pub;

	// PI
	static constexpr double PI = 2*acos(0.0);

	// number of consecutive accelerometer measurments to declare robot is stationary
	const int NUM_STATIONARY = 125;

	// acceleration threshold to detect if robot is stationary
	static constexpr double ACCEL_THRESH = 0.1; // m/s^2

	// is the rover stationary
	bool m_rover_stationary = false; 

	// stationary counts
	int m_accel_counter = 0;

	// noise parameters for random walk
	const double m_lambda_g = -1.0/100;
	const double m_lambda_a = -1.0/100;

	enum SENSOR_TYPE
	{
		SUN_SENSOR = 1,
		ACCELEROMETER = 2
	};

	void sun_sensor_callback(const geometry_msgs::Vector3::ConstPtr& msg);

	void imu_callback(const sensor_msgs::Imu::ConstPtr& msg);

	// F matrix. dx_dot = F*dx + G*w 
	void computeF(const Eigen::Matrix<double,3,1> &f_i, 
		const Eigen::Matrix<double,3,3> &R_body_to_nav_next, 
		Eigen::Matrix<double,9,9> &F) const;

	// G matrix. dx_dot = F*dx + G*w 
	void computeG(const Eigen::Matrix<double,3,3> &R_body_to_nav_next, 
		Eigen::Matrix<double,9,12> &G) const;

	// discrete state transition matrix and noise matrix
	void computePhiAndQdk(const Eigen::Matrix<double,3,1> &f_i, 
		const Eigen::Matrix<double,3,3> &R_body_to_nav_next, 
		Eigen::Matrix<double,9,9> &Phi, 
		Eigen::Matrix<double,9,9> &Qdk) const;

	// measurement update using gravity prediction
	void stationaryMeasurementUpdate(const Eigen::Matrix<double,3,3> & R_body_to_nav);

	// general Kalman filter
	void EKF(const Eigen::MatrixXd & H, 
		const Eigen::MatrixXd & R, 
		const Eigen::MatrixXd & z, 
		const SENSOR_TYPE sensor_type);

	// convert vector to 3x3 skew symmetric matrix
	void to_skew(const Eigen::Matrix<double,3,1> &v, Eigen::Matrix<double,3,3> &m) const
	{
		m << 0, -v(2), v(1),
		v(2), 0, -v(0),
		-v(1), v(0), 0;
	}
};
