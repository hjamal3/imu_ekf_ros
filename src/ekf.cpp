// HJ: August 3, 2020
// EKF adopted from earlier python code, fuses IMU and encoder data 
// ***I didn't use typedefs for clarity for readers
#include "ekf.h"
#include "imu_ekf_ros/initRequest.h"
#include "ros/ros.h"
#include "sensor_msgs/Imu.h"
#include "std_msgs/Int32MultiArray.h"
#include "std_msgs/Float64.h"
#include "geometry_msgs/Quaternion.h"
#include <tf/transform_broadcaster.h>
#include <cstdlib>

ros::Publisher quat_pub;
static int counter = 0;

void debug(auto str)
{
	std::cout << str << std::endl;
}

// convert vector to 3x3 skew symmetric matrix
void to_skew(const Eigen::Matrix<double,3,1> &v, Eigen::Matrix<double,3,3> &m)
{
	m << 0, -v(2), v(1),
	v(2), 0, -v(0),
	-v(1), v(0), 0;
}

void computeF(const Eigen::Matrix<double,3,1> &f_i, const Eigen::Matrix<double,3,3> &R_body_to_nav_next, Eigen::Matrix<double,9,9> &F)
{
	// store relevant values in F
	static Eigen::Matrix<double,3,3> f_i_skew(3,3);
	to_skew(f_i,f_i_skew);
	const double lambda_g = -1.0/100;
	const double lambda_a = -1.0/100;
	F.block<3,3>(0,3) = -1*R_body_to_nav_next; // where does the -1 come from???
	F.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Identity()*lambda_g; // Fg
	F.block<3,3>(6,6) = Eigen::Matrix<double,3,3>::Identity()*lambda_a; // Fa
}

void computeG(const Eigen::Matrix<double,3,3> &R_body_to_nav_next, Eigen::Matrix<double,9,12> &G)
{
	G.block<3,3>(0,0) = -1*R_body_to_nav_next;
	G.block<3,3>(0,6) = -1*R_body_to_nav_next;
	G.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Identity();
	G.block<3,3>(6,9) = Eigen::Matrix<double,3,3>::Identity();
}

void computePhiAndQdk(const Eigen::Matrix<double,3,1> &f_i, const Eigen::Matrix<double,3,3> &R_body_to_nav_next, Eigen::Matrix<double,9,9> &Phi, Eigen::Matrix<double,9,9> &Qdk)
{
	// van loan algorithm: https://www.cs.cornell.edu/cv/ResearchPDF/computing.integrals.involving.Matrix.Exp.pdf
	// formulate larger matrix
	// compute F matrix: F is 9 x 9
	static Eigen::Matrix<double,9,9> F = Eigen::Matrix<double,9,9>::Zero();;
	computeF(f_i, R_body_to_nav_next, F);

	// compute G. G is 9 x 12.
	static Eigen::Matrix<double,9,12> G = Eigen::Matrix<double,9,12>::Zero();
	computeG(R_body_to_nav_next, G);

	// empty matrix
	static Eigen::Matrix<double,18,18> A = Eigen::Matrix<double,18,18>::Zero();

	// fill in A matrix according to algorithm
	A.block<9,9>(0,0) = -F;
	A.block<9,9>(9,9) = F.transpose();
	A.block<9,9>(0,9) = G*(filter.Q)*G.transpose();
	A = A*filter.dt;

	// matrix exponential
	Eigen::Matrix<double,18,18> B = A.exp();
	Phi = B.block<9,9>(9,9).transpose();
	Qdk = Phi*B.block<9,9>(0,9);
}

void stationaryMeasurementUpdate(const Eigen::Matrix<double,3,3> & R_body_to_nav)
{
	//debug("stationary update!");
	rover_stationary = false;

	// stationary update matrix: Ha is 3 x 9
	static Eigen::Matrix<double,3,9> Ha = Eigen::Matrix<double,3,9>::Zero();

	// known gravity
	Eigen::Matrix<double,3,1> g_meas(0,0,filter.g);
	static Eigen::Matrix<double,3,3> g_skew(3,3);
	to_skew(g_meas, g_skew);
	Ha.block<3,3>(0,0) = -1*g_skew;
	Ha.block<3,3>(0,6) = R_body_to_nav;

	// measurement and noise matrices
	static Eigen::MatrixXd y_pred = Eigen::Matrix<double,3,1>::Zero(); // predicted measurement
	static Eigen::MatrixXd y_meas = Eigen::Matrix<double,3,1>::Zero(); // actual measurement 

	// predicted measurement is gravity 
	y_pred << g_pred(0), g_pred(1), g_pred(2);

	// actual measurement is actual gravity
	y_meas << g_meas(0), g_meas(1), g_meas(2);

	// predicted gravity: g_pred
	Eigen::Matrix<double,3,1> z = y_meas - y_pred;

	// call filter... TO DO. 3 or 2???
	EKF(Ha.block<2,9>(0,0),filter.Ra.block<2,2>(0,0),z(Eigen::seq(0,1)));
}

// input: measured heading of rover with respect to inertial z axis
void sun_sensor_callback(const std_msgs::Float64::ConstPtr& msg)
{
	// process ros data and compute unit vector. in reality this converts sun sensor readings to a 3D sun vecor
	// for now this will just give ground truth sun for testing
	double heading = msg->data;
	double x = cos(heading);
	double y = sin(heading);

	// process current orientation and measure predicted heading (x-axis)
	Eigen::Quaternion<double> b = Eigen::Quaternion<double>(state(0),state(1),state(2),state(3));
	Eigen::Vector3d s_b;
	s_b << 1, 0, 0; // you predict the sun straight ahead
	Eigen::Vector3d s_i = b.inverse()._transformVector(s_b);

	// stationary update matrix: Hs is 3 x 9
	static Eigen::Matrix<double,3,9> Hs = Eigen::Matrix<double,3,9>::Zero();

	// // noise matrix
	static Eigen::Matrix<double,3,3> Rs(3,3);
	Rs << 0.2, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.2;

	// known sensor reading
	Eigen::Matrix<double,3,1> s_meas(x,y,0);
	Eigen::Matrix<double,3,3> s_skew(3,3);
	to_skew(s_meas, s_skew);
	Hs.block<3,3>(0,0) = -1*s_skew;

	// // variables to fill -> I kept it like this for reviewer clarity. TODO refactor to shorten.
	Eigen::MatrixXd H; // measurement matrix
	Eigen::MatrixXd R; // noise matrix
	Eigen::MatrixXd y_pred; // predicted measurement
	Eigen::MatrixXd y_meas; // actual measurement 

	// measurement and noise matrices
	H = Hs;
	R = Rs;
	y_pred = Eigen::Matrix<double,3,1>::Zero();
	y_meas = Eigen::Matrix<double,3,1>::Zero();

	// predicted heading  
	y_pred << s_i(0), s_i(1), s_i(2);

	// measurement (from ephemeris)
	y_meas << s_meas(0), s_meas(1), s_meas(2);

	// in reality i would use the sensor reading to predict the unit vector and compare 
	// it to what is in the solar model. this is just for testing.

	// predicted gravity: g_pred
	Eigen::Matrix<double,3,1> z = y_meas - y_pred;

	// only correct heading herror
	EKF(H.block<3,9>(0,0),R.block<3,3>(0,0),z.block<3,1>(0,0));
	// EKF(H.block<1,9>(2,0),R.block<1,1>(0,0),z.block<1,1>(2, 0));
	// std::cout << H.block<1,9>(2,0) << std::endl;
	std::cout << cov(2,2) << std::endl;
}

// updates state. general to measurements given appropriately sized 
void EKF(const Eigen::MatrixXd & H, const Eigen::MatrixXd & R, const Eigen::MatrixXd & z)
{

	// compute Kalman gain
	Eigen::MatrixXd K = cov*H.transpose()*((R+H*cov*H.transpose()).inverse()); // fix auto later

	// compute state correction
	Eigen::Matrix<double,9,1> dx = K*z;

	// 	rotation update 
	Eigen::Matrix<double,3,1> ro = dx(Eigen::seq(0,2));
	static Eigen::Matrix<double,3,3> P(3,3);
	to_skew(ro,P);
	// predicted orientation

	Eigen::Quaternion<double> b_q = Eigen::Quaternion<double>(state[0],state[1],state[2],state[3]);
	Eigen::Matrix<double,3,3> R_nav_to_body = b_q.toRotationMatrix();
	Eigen::Matrix<double,3,3> R_nav_to_body_next = R_nav_to_body*(Eigen::Matrix<double,3,3>::Identity() - P); // (10.67)

	// reorthogonalize matrix (make it into a rotation matrix)... use SVD
	Eigen::JacobiSVD<Eigen::Matrix<double,3,3>> svd(R_nav_to_body_next,Eigen::DecompositionOptions::ComputeFullU | Eigen::DecompositionOptions::ComputeFullV);
	R_nav_to_body_next = svd.matrixU()*((svd.matrixV()).transpose());

	// 	compute next quaternion
	Eigen::Quaternion<double> b_next(R_nav_to_body_next);
	if (abs(ro[0]) > 0.1 || abs(ro[1]) > 0.1 || abs(ro[2]) > 0.1)
	{
		//debug("angle update high");
	}

	// update state
	state.block<4,1>(0,0) << b_next.w(), b_next.x(), b_next.y(), b_next.z(); // update quaternion: [w, x, y, z]
	state.block<3,1>(4,0) += dx(Eigen::seq(3,5)); // update accel bias
	state.block<3,1>(7,0) += dx(Eigen::seq(6,8)); // update gyro bias

	// update covariance
	cov = (Eigen::Matrix<double,9,9>::Identity()-K*H)*(cov);

	// symmetrify 
	cov = (cov + cov.transpose())/2.0;

	if (false)
	{
		Eigen::Quaternion<double> b_next_body_to_nav = b_next.inverse();
		static tf::TransformBroadcaster br;
		tf::Transform transform;
		transform.setOrigin( tf::Vector3(0, 0, 0));
		tf::Quaternion q(b_next_body_to_nav.x(),b_next_body_to_nav.y(),b_next_body_to_nav.z(),b_next_body_to_nav.w());
		transform.setRotation(q);
		br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "map", "IMU"));
	}
}

// imu callback
void imu_callback(const sensor_msgs::Imu::ConstPtr& msg)
{
	// Uncomment to time things
	//filter.timer.PrintDt("Ignore"); 
	//filter.timer.PrintDt("IMU"); 

	// IMU data
	geometry_msgs::Vector3 w = msg->angular_velocity;
	geometry_msgs::Vector3 f = msg->linear_acceleration;
	Eigen::Matrix<double,3,1> w_nb(w.x,w.y,w.z);
	Eigen::Matrix<double,3,1> f_b(f.x,f.y,f.z);

	/* Update State */

	// current orientation
	Eigen::Matrix<double,4,1> b = state(Eigen::seq(0,3));
	Eigen::Quaternion<double> b_prev = Eigen::Quaternion<double>(b(0),b(1),b(2),b(3)); // w, x, y, z

	// current accel bias
	Eigen::Matrix<double,3,1> x_g = state(Eigen::seq(4,6));

	// current gyro bias
	Eigen::Matrix<double,3,1> x_a = state(Eigen::seq(6,8));

	// subtract out gyroscope bias. also w_bn = (-w_nb)
	Eigen::Matrix<double,3,1> w_bn = -1*(w_nb-x_g);
	double w_norm = w_bn.norm();

	// subtract out accelerometer bias
	f_b = f_b - x_a;

	// differential rotation: [w, x, y, z]
	Eigen::Matrix<double,3,1> xyz = sin(w_norm*filter.dt/2.0)*w_bn/w_norm;
	Eigen::Quaternion<double> db = Eigen::Quaternion<double>(cos(w_norm*filter.dt/2.0), xyz[0], xyz[1], xyz[2]);

	// update orientation
	Eigen::Quaternion<double> b_next = db*b_prev;

	// get average quaternion by interpolation
	Eigen::Quaternion<double> b_avg = b_prev.slerp(0.5,b_next);

	// b is the nav to body transformation. we need body to nav transformation -> invert 
	Eigen::Quaternion<double> b_body_to_nav_avg = b_avg.inverse(); // for specific force (5.9 Principles of GNSS book)

	// rotate specific force into inertial frame
	Eigen::Matrix<double,3,1> f_i = b_body_to_nav_avg._transformVector(f_b);

	// gravity vector
	Eigen::Matrix<double,3,1> g_vec(0,0,filter.g);

	// get acceleration in inertial frame. (acceleration of body wrt inertial frame in inertial frame)
	Eigen::Matrix<double,3,1> a_i = f_i - g_vec;

	// store in state -> this is time propagation step. 
	state << b_next.w(),b_next.x(),b_next.y(),b_next.z(), x_g(0),x_g(1),x_g(2), x_a(0),x_a(1),x_a(2);

	/* Update covariance */
	Eigen::Matrix<double,3,3> R_body_to_nav_next = b_next.inverse().toRotationMatrix();

	// compute state transition matrix Phi
	// compute Qdk (discrete noise). Qdk is 9 x 9
	static Eigen::Matrix<double,9,9> Phi(9,9);
	static Eigen::Matrix<double,9,9> Qdk(9,9);
	computePhiAndQdk(f_i, R_body_to_nav_next, Phi, Qdk);

	// update covariance (15x15)
	cov = Phi*cov*Phi.transpose()+Qdk;

	/* Measurement update using accelerometer to correct roll and pitch */
	// if 50 accelerometer readings were close enough to the gravity vector, robot is stationary
	//  check if current measurement is stationary
	if (abs(f_b.norm() - filter.g) < stat_acces_thresh) // tuned. test: you should never be able to hold the imu in your hand and have an update.
	{
		accel_counter++;
		g_pred_sum += R_body_to_nav_next*(f_b); // R*(x_a-y_a) TODO: CHECK THIS

	} else {
		accel_counter = 0;
		g_pred_sum = Eigen::Matrix<double,3,1>::Zero();
		rover_stationary = false;
	}
	// if n consecutive stationary, use accel_data
	if (accel_counter == num_stat_measurements)
	{
		// predict gravity in navigation frame and store prediction in global variable.
		g_pred = g_pred_sum/num_stat_measurements; // averaging
		accel_counter = 0;
		g_pred_sum << 0,0,0;
		static int stat_ctr = 0;
		std::cout << "stationary update: " << stat_ctr << std::endl;
		stat_ctr++;
		stationaryMeasurementUpdate(R_body_to_nav_next);
	}

	if (true)
	{
		Eigen::Quaternion<double> b_next_body_to_nav = b_next.inverse();
		static tf::TransformBroadcaster br;
		tf::Transform transform;
		transform.setOrigin( tf::Vector3(0, 0, 0));
		tf::Quaternion q(b_next_body_to_nav.x(),b_next_body_to_nav.y(),b_next_body_to_nav.z(),b_next_body_to_nav.w());
		transform.setRotation(q);
		br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "map", "IMU"));
		geometry_msgs::Quaternion msg;
		msg.x = b_next_body_to_nav.x();
		msg.y = b_next_body_to_nav.y();
		msg.z = b_next_body_to_nav.z();
		msg.w = b_next_body_to_nav.w();
		quat_pub.publish(msg);
	}

	counter++;
}

void initialize_ekf(ros::NodeHandle &n)
{
	ROS_INFO("ekf: waiting for initialization service");

	// create client for service
	ros::ServiceClient client = n.serviceClient<imu_ekf_ros::initRequest>("/initialize_ekf");

	// instantiate service class
	imu_ekf_ros::initRequest srv;

	// publish angles
	quat_pub = n.advertise<geometry_msgs::Quaternion>("/quat", 0);

	// call the service
	if (!client.waitForExistence(ros::Duration(-1)))
	{
		ROS_ERROR("ekf: initialize_ekf didn't send data");
	} else 
	{
		ROS_INFO("ekf: connected to initialize service");
	}
	if (client.call(srv))
	{
		ROS_INFO("ekf: initialize_ekf responded with data.");
		// store received data
		geometry_msgs::Quaternion b = srv.response.init_orientation;
		Eigen::Vector3d x_g = {srv.response.gyro_bias[0].data, srv.response.gyro_bias[1].data, srv.response.gyro_bias[2].data};

		// filter rate parameters
		int num_data, hz; // number of data points used to initialize, imu hz
		n.param("num_data",num_data,125);
		n.param("imu_hz",hz,125);
		filter.dt = 1.0/(double)hz;
		filter.num_data = (double)num_data;
		double T = num_data/hz;  //number of measurements over rate of IMU

		// initialize noise terms
		double sigma_xg, sigma_nug, sigma_xa, sigma_nua;
		n.param<double>("sigma_xg",sigma_xg,0.00000290); // Gyro (rate) random walk
		n.param<double>("sigma_nug",sigma_nug,0.00068585); // rad/s/rt_Hz, Gyro white noise
		n.param<double>("sigma_xa",sigma_xa,0.00001483);  // Accel (rate) random walk m/s3 1/sqrt(Hz)
		n.param<double>("sigma_nua",sigma_nua,0.00220313); // accel white noise

		// noise matrix for IMU (Q)
		for (int i = 0; i < 3; i++)
		{
			filter.Q(i,i) = sigma_nua*sigma_nua;
			filter.Q(3+i,3+i) = sigma_xg*sigma_xg;
			filter.Q(6+i,6+i) = sigma_nug*sigma_nug;
			filter.Q(9+i,9+i) = sigma_xa*sigma_xa;
		}

		// noise matrix for accelerometer (Ra)
		filter.Ra = Eigen::Matrix<double,3,3>::Identity(3,3)*(sigma_nua*sigma_nua);

		// gravity vector
		n.param<double>("g",filter.g,9.80665);

		// initialize state: [b, x_a, x_g] = [quaternion, gyro bias, accel bias],  size 10
		state << b.w,b.x,b.y,b.z,x_g[0],x_g[1],x_g[2], 0,0,0;

		// initialize covariance
		cov.block<2,2>(0,0) = (sigma_nua/filter.g)*(sigma_nua/filter.g)/T*Eigen::Matrix<double,2,2>::Identity();
		cov(2,2) = 10.0;
		cov.block<3,3>(3,3) = (sigma_nug)*(sigma_nug)/T*Eigen::Matrix<double,3,3>::Identity();

		// TODO: finish this part
		//cov(0,7) =  cov[0,0]*cov[7,7]
		//cov.block
		//cov[0:2,7:9] = np.diag([math.sqrt(cov[0,0]*cov[7,7]),math.sqrt(cov[1,1]*cov[8,8])]) # absolutely no idea
		//cov[6:,0:3] = np.transpose(cov[0:3,6:])
	}
	else
	{
		ROS_ERROR("ekf: Failed to call service initialize_ekf.");
		ros::shutdown();
	}
}


int main(int argc, char **argv)
{
	ROS_INFO("EKF node started.");

	ros::init(argc, argv, "imu_ekf_node");

	ros::NodeHandle n;

	// initialize ekf
	initialize_ekf(n);

	// imu callback
	ros::Subscriber sub_imu = n.subscribe("/imu/data", 1000, imu_callback);

	// sun sensor callback
	ros::Subscriber sub_sunsensor = n.subscribe("/sun_sensor", 1, sun_sensor_callback);

	ros::spin();

	std::cout << counter << std::endl;

	return 0;
}
