// HJ: August 3, 2020
// EKF adopted from earlier python code, fuses IMU and encoder data 
// ***I didn't use typedefs for clarity for readers
#include "ekf.h"

void IMU_EKF::computeF(const Eigen::Matrix<double,3,1> &f_i, 
	const Eigen::Matrix<double,3,3> &R_body_to_nav_next, 
	Eigen::Matrix<double,9,9> &F) const
{
	// store relevant values in F
	static Eigen::Matrix<double,3,3> f_i_skew(3,3);
	to_skew(f_i,f_i_skew);
	F.block<3,3>(0,3) = -1*R_body_to_nav_next; // where does the -1 come from???
	F.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Identity()*m_lambda_g; // Fg
	F.block<3,3>(6,6) = Eigen::Matrix<double,3,3>::Identity()*m_lambda_a; // Fa
}

void IMU_EKF::computeG(const Eigen::Matrix<double,3,3> &R_body_to_nav_next, 
	Eigen::Matrix<double,9,12> &G) const
{
	G.block<3,3>(0,0) = -1*R_body_to_nav_next;
	G.block<3,3>(0,6) = -1*R_body_to_nav_next;
	G.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Identity();
	G.block<3,3>(6,9) = Eigen::Matrix<double,3,3>::Identity();
}

void IMU_EKF::computePhiAndQdk(const Eigen::Matrix<double,3,1> &f_i, 
	const Eigen::Matrix<double,3,3> &R_body_to_nav_next, 
	Eigen::Matrix<double,9,9> &Phi, 
	Eigen::Matrix<double,9,9> &Qdk) const
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
	A.block<9,9>(0,9) = G*(m_filter.Q)*G.transpose();
	A = A*m_filter.dt;

	// matrix exponential
	Eigen::Matrix<double,18,18> B = A.exp();
	Phi = B.block<9,9>(9,9).transpose();
	Qdk = Phi*B.block<9,9>(0,9);
}

void IMU_EKF::stationaryMeasurementUpdate(const Eigen::Matrix<double,3,3> & R_body_to_nav)
{
	//EKF::debug("stationary update!");
	m_rover_stationary = false;

	// stationary update matrix: Ha is 3 x 9
	static Eigen::Matrix<double,3,9> Ha = Eigen::Matrix<double,3,9>::Zero();

	// known gravity
	Eigen::Matrix<double,3,1> g_meas(0,0,m_filter.g);
	static Eigen::Matrix<double,3,3> g_skew(3,3);
	to_skew(g_meas, g_skew);
	Ha.block<3,3>(0,0) = -1*g_skew;
	Ha.block<3,3>(0,6) = R_body_to_nav;

	// measurement and noise matrices
	static Eigen::MatrixXd y_pred = Eigen::Matrix<double,3,1>::Zero(); // predicted measurement
	static Eigen::MatrixXd y_meas = Eigen::Matrix<double,3,1>::Zero(); // actual measurement 

	// predicted measurement is gravity 
	y_pred << m_g_pred(0), m_g_pred(1), m_g_pred(2);

	// actual measurement is actual gravity
	y_meas << g_meas(0), g_meas(1), g_meas(2);

	// predicted gravity: m_g_pred
	Eigen::Matrix<double,3,1> z = y_meas - y_pred;

	// call filter... TO DO. 3 or 2???
	EKF(Ha.block<2,9>(0,0),m_filter.Ra.block<2,2>(0,0),z(Eigen::seq(0,1)), ACCELEROMETER);
}

// input: a 3D vector in the inertial frame of the rover
// for testing, using the surveyed x-axis of the robot
// in reality, the measured sun ray
void IMU_EKF::sun_sensor_callback(const geometry_msgs::Vector3::ConstPtr& msg)
{
	std::cout << "Sun sensor update" << std::endl;
	// input is a 3D sun ray in the inertial frame
	Eigen::Quaternion<double> b = Eigen::Quaternion<double>(m_state(0),m_state(1),m_state(2),m_state(3));
	Eigen::Quaternion<double> b_body_to_nav = b.inverse();

	// stationary update matrix: Hs is 3 x 9
	static Eigen::Matrix<double,3,9> Hs = Eigen::Matrix<double,3,9>::Zero();

	// noise matrix
	static Eigen::Matrix<double,3,3> Rs(3,3);
	const double angle_uncertainty = 0.1; // rad^2
	Rs << angle_uncertainty/180.0*PI, 0.0, 0.0, 0.0, angle_uncertainty/180.0*PI, 0.0, 0.0, 0.0, angle_uncertainty/180.0*PI;

	// known sensor reading
	Eigen::Matrix<double,3,1> s_pred = b_body_to_nav._transformVector({1,0,0}); // predict the x-axis. in place of true sensor reading
	Eigen::Matrix<double,3,1> s_meas = {msg->x, msg->y, msg->z};
	s_meas = s_meas/s_meas.norm(); // normalize the model output

	static Eigen::Matrix<double,3,3> s_skew(3,3);
	to_skew(s_meas, s_skew);
	Hs.block<3,3>(0,0) = -1*s_skew;

	// no corrections to tilt and roll
	Hs(0,1) = 0.0;
	Hs(1,0) = 0.0;

	static Eigen::MatrixXd y_pred = Eigen::Matrix<double,3,1>::Zero();
	static Eigen::MatrixXd y_meas = Eigen::Matrix<double,3,1>::Zero();

	// predicted heading  
	y_pred << s_pred(0), s_pred(1), s_pred(2);

	// measurement (from ephemeris)
	y_meas << s_meas(0), s_meas(1), s_meas(2);

	// predicted gravity: m_g_pred
	Eigen::Matrix<double,3,1> z = y_meas - y_pred;

	// only correct heading error
	EKF(Hs.block<2,9>(0,0),Rs.block<2,2>(0,0),z.block<2,1>(0,0), SUN_SENSOR);
}

// updates state. general to measurements given appropriately sized 
void IMU_EKF::EKF(const Eigen::MatrixXd & H, const Eigen::MatrixXd & R, const Eigen::MatrixXd & z, const SENSOR_TYPE sensor_type)
{

	// compute Kalman gain
	Eigen::MatrixXd K = m_cov*H.transpose()*((R+H*m_cov*H.transpose()).inverse()); // fix auto later

	// compute state correction
	Eigen::Matrix<double,9,1> dx = K*z;

	// rotation update 
	Eigen::Matrix<double,3,1> ro = dx(Eigen::seq(0,2));
	static Eigen::Matrix<double,3,3> P(3,3);

	// only correct the z-axis (yaw)
	if (sensor_type == SUN_SENSOR)
	{
		ro(0) = 0;
		ro(1) = 0;
	}
	to_skew(ro,P);

	// predicted orientation
	Eigen::Quaternion<double> b_q = Eigen::Quaternion<double>(m_state[0],m_state[1],m_state[2],m_state[3]);
	Eigen::Matrix<double,3,3> R_nav_to_body = b_q.toRotationMatrix();
	Eigen::Matrix<double,3,3> R_nav_to_body_next = R_nav_to_body*(Eigen::Matrix<double,3,3>::Identity() - P); // (10.67)

	// reorthogonalize matrix (make it into a rotation matrix)... use SVD
	Eigen::JacobiSVD<Eigen::Matrix<double,3,3>> svd(R_nav_to_body_next,Eigen::DecompositionOptions::ComputeFullU | Eigen::DecompositionOptions::ComputeFullV);
	R_nav_to_body_next = svd.matrixU()*((svd.matrixV()).transpose());

	// 	compute next quaternion
	Eigen::Quaternion<double> b_next(R_nav_to_body_next);

	// update state
	m_state.block<4,1>(0,0) << b_next.w(), b_next.x(), b_next.y(), b_next.z(); // update quaternion: [w, x, y, z]
	m_state.block<3,1>(4,0) += dx(Eigen::seq(3,5)); // update accel bias
	m_state.block<3,1>(7,0) += dx(Eigen::seq(6,8)); // update gyro bias

	// update covariance
	m_cov = (Eigen::Matrix<double,9,9>::Identity()-K*H)*(m_cov);

	// symmetrify 
	m_cov = (m_cov + m_cov.transpose())/2.0;
}

// imu callback
void IMU_EKF::imu_callback(const sensor_msgs::Imu::ConstPtr& msg)
{
	// Uncomment to time things
	//m_filter.timer.PrintDt("Ignore"); 
	//m_filter.timer.PrintDt("IMU"); 

	// IMU data
	geometry_msgs::Vector3 w = msg->angular_velocity;
	geometry_msgs::Vector3 f = msg->linear_acceleration;
	Eigen::Matrix<double,3,1> w_nb(w.x,w.y,w.z);
	Eigen::Matrix<double,3,1> f_b(f.x,f.y,f.z);

	/* Update State */

	// current orientation
	Eigen::Matrix<double,4,1> b = m_state(Eigen::seq(0,3));
	Eigen::Quaternion<double> b_prev = Eigen::Quaternion<double>(b(0),b(1),b(2),b(3)); // w, x, y, z

	// current accel bias
	Eigen::Matrix<double,3,1> x_g = m_state(Eigen::seq(4,6));

	// current gyro bias
	Eigen::Matrix<double,3,1> x_a = m_state(Eigen::seq(6,8));

	// subtract out gyroscope bias. also w_bn = (-w_nb)
	Eigen::Matrix<double,3,1> w_bn = -1*(w_nb-x_g);
	double w_norm = w_bn.norm();

	// subtract out accelerometer bias
	f_b = f_b - x_a;

	// differential rotation: [w, x, y, z]
	Eigen::Matrix<double,3,1> xyz = sin(w_norm*m_filter.dt/2.0)*w_bn/w_norm;
	Eigen::Quaternion<double> db = Eigen::Quaternion<double>(cos(w_norm*m_filter.dt/2.0), xyz[0], xyz[1], xyz[2]);

	// update orientation
	Eigen::Quaternion<double> b_next = db*b_prev;

	// get average quaternion by interpolation
	Eigen::Quaternion<double> b_avg = b_prev.slerp(0.5,b_next);

	// b is the nav to body transformation. we need body to nav transformation -> invert 
	Eigen::Quaternion<double> b_body_to_nav_avg = b_avg.inverse(); // for specific force (5.9 Principles of GNSS book)

	// rotate specific force into inertial frame
	Eigen::Matrix<double,3,1> f_i = b_body_to_nav_avg._transformVector(f_b);

	// get acceleration in inertial frame. (acceleration of body wrt inertial frame in inertial frame)
	Eigen::Matrix<double,3,1> a_i = f_i - m_g_true;

	// store in state -> this is time propagation step. 
	m_state << b_next.w(),b_next.x(),b_next.y(),b_next.z(), x_g(0),x_g(1),x_g(2), x_a(0),x_a(1),x_a(2);

	/* Update covariance */
	Eigen::Matrix<double,3,3> R_body_to_nav_next = b_next.inverse().toRotationMatrix();

	// compute state transition matrix Phi
	// compute Qdk (discrete noise). Qdk is 9 x 9
	static Eigen::Matrix<double,9,9> Phi(9,9);
	static Eigen::Matrix<double,9,9> Qdk(9,9);
	computePhiAndQdk(f_i, R_body_to_nav_next, Phi, Qdk);

	// update covariance (15x15)
	m_cov = Phi*m_cov*Phi.transpose()+Qdk;

	/* Measurement update using accelerometer to correct roll and pitch */
	// if 50 accelerometer readings were close enough to the gravity vector, robot is stationary
	//  check if current measurement is stationary
	if (abs(f_b.norm() - m_filter.g) < ACCEL_THRESH) // tuned. test: you should never be able to hold the imu in your hand and have an update.
	{
		m_accel_counter++;
		m_g_pred_sum += R_body_to_nav_next*(f_b); // R*(x_a-y_a) TODO: CHECK THIS

	} else {
		m_accel_counter = 0;
		m_g_pred_sum = Eigen::Matrix<double,3,1>::Zero();
		m_rover_stationary = false;
	}
	// if n consecutive stationary, use accel_data
	if (m_accel_counter == NUM_STATIONARY)
	{
		// predict gravity in navigation frame and store prediction in global variable.
		m_g_pred = m_g_pred_sum/NUM_STATIONARY; // averaging
		m_accel_counter = 0;
		m_g_pred_sum << 0,0,0;
		stationaryMeasurementUpdate(R_body_to_nav_next);
	}

	// Publish TF and orientation
	Eigen::Quaternion<double> b_next_body_to_nav = b_next.inverse();
	static tf::TransformBroadcaster br;
	tf::Transform transform;
	transform.setOrigin( tf::Vector3(0, 0, 0));
	tf::Quaternion q(b_next_body_to_nav.x(),b_next_body_to_nav.y(),b_next_body_to_nav.z(),b_next_body_to_nav.w());
	transform.setRotation(q);
	br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "map", "IMU"));
	geometry_msgs::Quaternion quat_msg;
	quat_msg.x = b_next_body_to_nav.x();
	quat_msg.y = b_next_body_to_nav.y();
	quat_msg.z = b_next_body_to_nav.z();
	quat_msg.w = b_next_body_to_nav.w();
	m_orientation_pub.publish(quat_msg);
}

void IMU_EKF::initialize_ekf(ros::NodeHandle &n)
{
	ROS_INFO("ekf: waiting for initialization service");

	// create client for service
	ros::ServiceClient client = n.serviceClient<imu_ekf_ros::initRequest>("/initialize_ekf");

	// instantiate service class
	imu_ekf_ros::initRequest srv;

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

		// initial orientation
		geometry_msgs::Quaternion b = srv.response.init_orientation;

		// initial gyro biases
		Eigen::Vector3d x_g = {srv.response.gyro_bias[0].data, 
			srv.response.gyro_bias[1].data, 
			srv.response.gyro_bias[2].data};

		// filter rate parameters
		int num_data = 0; // number of data points used to initialize orientation
		int hz = 0; // imu freqency
		n.param("num_data",num_data,125);
		n.param("imu_hz",hz,125);
		m_filter.dt = 1.0/(double)hz;
		m_filter.num_data = (double)num_data;
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
			m_filter.Q(i,i) = sigma_nua*sigma_nua;
			m_filter.Q(3+i,3+i) = sigma_xg*sigma_xg;
			m_filter.Q(6+i,6+i) = sigma_nug*sigma_nug;
			m_filter.Q(9+i,9+i) = sigma_xa*sigma_xa;
		}

		// noise matrix for accelerometer (Ra)
		m_filter.Ra = Eigen::Matrix<double,3,3>::Identity(3,3)*(sigma_nua*sigma_nua);

		// gravity vector
		n.param<double>("g",m_filter.g,9.80665);

		// initialize gravity vector
		m_g_true << 0, 0, m_filter.g;

		// initialize state: [b, x_a, x_g] = [quaternion, gyro bias, accel bias],  size 10
		m_state << b.w,b.x,b.y,b.z,x_g[0],x_g[1],x_g[2], 0,0,0;

		// initialize covariance
		m_cov.block<2,2>(0,0) = (sigma_nua/m_filter.g)*(sigma_nua/m_filter.g)/T*Eigen::Matrix<double,2,2>::Identity();
		m_cov(2,2) = 1000*PI/180.0; // yaw has large uncertainty initially
		m_cov.block<3,3>(3,3) = (sigma_nug)*(sigma_nug)/T*Eigen::Matrix<double,3,3>::Identity();

		// TODO: finish this part
		//m_cov(0,7) =  m_cov[0,0]*m_cov[7,7]
		//m_cov.block
		//m_cov[0:2,7:9] = np.diag([math.sqrt(cov[0,0]*cov[7,7]),math.sqrt(cov[1,1]*cov[8,8])]) # absolutely no idea
		//m_cov[6:,0:3] = np.transpose(cov[0:3,6:])
	}
	else
	{
		ROS_ERROR("ekf: Failed to call service initialize_ekf.");
		ros::shutdown();
	}
}


IMU_EKF::IMU_EKF(ros::NodeHandle & n)
{
	// imu callback
	m_imu_subscriber =  n.subscribe("/imu/data", 1000, &IMU_EKF::imu_callback, this);

	// sun sensor callback
	m_sun_sensor_subscriber = n.subscribe("/sun_sensor", 1, &IMU_EKF::sun_sensor_callback, this);

	// orientation publisher
	m_orientation_pub = n.advertise<geometry_msgs::Quaternion>("/quat", 0, this);

	// initialize ekf
	initialize_ekf(n);
}

int main(int argc, char **argv)
{
	ROS_INFO("EKF node started.");

	ros::init(argc, argv, "imu_ekf_node");

	ros::NodeHandle n;

	IMU_EKF imu_ekf(n);

	ros::spin();

	return 0;
}
