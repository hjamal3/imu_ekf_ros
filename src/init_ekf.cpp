#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "imu_ekf_ros/initRequest.h"
#include "math.h"
#include "ros/ros.h"
#include "sensor_msgs/Imu.h"
#include "std_msgs/Float64.h"
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>

// debugging
#include <iostream>

// data collection parameters
int num_data = 0;
int imu_counter = 0;

// data storage elements
Eigen::Vector3d sum_accel;
Eigen::Vector3d sum_gyro;

// public node handle pointer 
ros::NodeHandle* n_ptr;
ros::Subscriber* sub_imu_ptr;
ros::ServiceServer service;

bool handle_init_ekf(imu_ekf_ros::initRequest::Request  &req, imu_ekf_ros::initRequest::Response &res)
{
	// compute initial orientation
	Eigen::Vector3d g_b = sum_accel / num_data;

	// initial roll (phi) and pitch (theta)
	double phi = atan2(-g_b[1],-g_b[2]);
	double theta = atan2(g_b[0], sqrt(g_b[1]*g_b[1] + g_b[2]*g_b[2]));

	// set initial yaw to zero
	double psi = 0;

	// q is navigation to body transformation: R_bi
	// YPR: R_ib = R(yaw)R(pitch)R(Roll)
	// RPY: R_bi = R(-Roll)R(-Pitch)R(-yaw)	
	Eigen::Quaternion<double> q = Eigen::AngleAxisd(psi, Eigen::Vector3d::UnitZ())
    * Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitY())
    * Eigen::AngleAxisd(phi, Eigen::Vector3d::UnitX());
    q = q.inverse(); 

	// compute gyroscope biases
	Eigen::Vector3d gyro_biases = sum_gyro / num_data;

	// store in response
	res.gyro_bias[0].data = gyro_biases[0];
	res.gyro_bias[1].data = gyro_biases[1];
	res.gyro_bias[2].data = gyro_biases[2];
	res.init_orientation.x = q.x();
	res.init_orientation.y = q.y();
	res.init_orientation.z = q.z();
	res.init_orientation.w = q.w();

	ROS_INFO("init_ekf: processed response");

	if (true)
	{
		std::cout << gyro_biases << std::endl;
		std::cout << q << std::endl;
	}
}

void imu_callback(const sensor_msgs::Imu::ConstPtr& msg)
{
	if (imu_counter < num_data)
	{
		// get accelerometer data
		geometry_msgs::Vector3 a = msg->linear_acceleration;

		// get gyroscope data
		geometry_msgs::Vector3 w = msg->angular_velocity;

		// add to matrix
		sum_accel -= Eigen::Vector3d(a.x,a.y,a.z);
		sum_gyro += Eigen::Vector3d(w.x,w.y,w.z);

		// increment counter
		imu_counter++;

	} else 
	{
		// stop receiving new imu data
		(*sub_imu_ptr).shutdown();

		// create the service
		service = (*n_ptr).advertiseService("/initialize_ekf", handle_init_ekf);
	}
}


// service handle

int main(int argc, char **argv)
{
	ROS_INFO("init_ekf node started.");

	ros::init(argc, argv, "init_imu_ekf_node");

	// create node handle and pointer
	ros::NodeHandle n;
	n_ptr = &n;

	// get number of data items to average from parameter server
	n.param("num_data", num_data, 500);

	// imu callback
	ros::Subscriber sub_imu = n.subscribe("/imu/data", 10, imu_callback);
	sub_imu_ptr = &sub_imu;
	ros::spin();
}
