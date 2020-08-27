#!/usr/bin/env python
# HJ: 5/19/20
# ROS 
import rospy
from sensor_msgs.msg import Imu
from geometry_msgs.msg import Vector3Stamped
from imu_ekf_ros.srv import *

# other
import numpy as np
import math
from pyquaternion import Quaternion
from scipy.linalg import expm
from helper import to_skew

#debugging
from geometry_msgs.msg import PoseStamped
from geometry_msgs.msg import Quaternion as quatmsg
import tf

counter = 0
# accelerometer
accel_counter = 0


# TODO: Finish EKF with only accel
def measurement_update():

	# access current state
	global state
	b = state[0]
	x_g = state[1]
	x_a = state[2]
	b_q = Quaternion(array=b)
	R_body_to_nav = np.transpose(b_q.rotation_matrix)

	# measurement matrix gyroscope
	g_const = 9.81
	Ha = np.array([[0,g_const,0,0,0,0],[-g_const, 0,0,0,0,0],[0,0,0,0,0,0]]) 	
	Ha = np.hstack((Ha,R_body_to_nav))

	# measurement matrices
	H = Ha

	# previous covariance: 9x9
	global cov

	# compute individual sensor variances
	global sigma_xg, sigma_nug, sigma_xa, sigma_nua

	# noise matrix for measurement
	R = np.array([[sigma_nua**2,0,0],[0,sigma_nua**2,0],[0,0,sigma_nua**2]])

	# compute Kalman gain

	K = cov.dot(H.T).dot(np.linalg.inv(R+H.dot(cov).dot(H.T)))

	# access current sensor readings and stack them
	global g_pred

	y_pred = g_pred

	# known gravity vector
	g = np.array([0,0,9.80])
	y = g

	# residual z
	z = y-y_pred

	# compute state correction
	dx = K.dot(z)

	# rotation update 
	p = dx[:3]
	P = to_skew(p)
	R_nav_to_body = R_body_to_nav.T
	R_nav_to_body_next = R_nav_to_body.dot(np.identity(3) - P) # (10.67)

	# reorthogonalize matrix... required by PyQuaternion... use SVD
	U, S, Vh = np.linalg.svd(R_nav_to_body_next)
	R_nav_to_body_next = U.dot(Vh)

	b_next = Quaternion(matrix=R_nav_to_body_next).elements

	# bias update
	dx_g = dx[3:6]
	dx_a = dx[6:]

	# update state 
	state[0] = b_next # update quaternion
	state[1] = state[1] + dx_g # update gyro bias
	state[2] = state[2] + dx_a  # update accel bias

	# update covariance
	cov = (np.identity(9)-K.dot(H)).dot(cov)


# update orientation here. use RK3
def imu_callback(data):

	# sample time. watch out for truncation
	global dt

	# get angular velocity and accelerometer data
	w = np.array([data.angular_velocity.x, data.angular_velocity.y,data.angular_velocity.z])
	a = np.array([data.linear_acceleration.x,data.linear_acceleration.y,data.linear_acceleration.z])

	# current state: [b, x_g, x_a], [quaternion, gyro bias, accel bias]
	global state
	b = state[0]
	x_g = state[1]
	x_a = state[2]

	# subtract out gyroscope bias. measured angular velocity is w_nb. Rot body wrt nav in body frame.
	# we need w_bn = (-w_nb)
	w = -1*(w - x_g)
	w_norm = np.linalg.norm(w)

	# get db, the differential rotation (quaternion)
	# quaternion is stored as [w, x, y, z]
	db = np.concatenate(([math.cos(w_norm*dt/2)],math.sin(w_norm*dt/2)*w/w_norm))

	# convert to PyQuaternion format
	b_prev = Quaternion(array=b)
	db = Quaternion(array=db)

	# multiply previous quaternion by differential to get latest
	b_next = db*b_prev

	# update state (convert to numpy array) -> this is time propagation step. 
	state[0] = b_next.elements

	# publish Quaternion and also tf
	if True:
		b_vis = b_next.inverse
		p = quatmsg()
		p.w = b_vis[0]
		p.x = b_vis[1]
		p.y = b_vis[2]
		p.z = b_vis[3]
		pub = rospy.Publisher('AHRS_EKF_quat', quatmsg, queue_size=1)
		pub.publish(p)

		# useful for visualization
		br = tf.TransformBroadcaster()
		# world to odom transformation
		br.sendTransform((0, 0, 0),
		(b_vis[1],b_vis[2],b_vis[3],b_vis[0]),
		rospy.Time.now(),
		"IMU",
		"ENU") # NED to quat

	# get rotation matrix from b (body to navigation frame) (takes points in body frame and puts them in nav frame) Rnb
	R_body_to_nav = np.transpose(b_next.rotation_matrix)

	# compute F. F is 9x9. xdot = Fx + Gw
	F = np.zeros((9,9))
	lambda_g = 1.0/1000 # correlation time of gyro (from stochastic bias model)
	lambda_a = 1.0/1000 # correlation time of accel (from stochastic bias model)
	F_g = -1.0/1000*np.identity(3)
	F_a = -1.0/1000*np.identity(3)
	F[:3,3:6] = -R_body_to_nav
	F[3:6,3:6] = F_g
	F[6:,6:] = F_a

	# compute state transition matrix (Phi) (9x9)
	Phi = expm(F*dt)

	# Compute G. G is 9 x 12. might get rid of the noise due to earth rotation later
	G = np.zeros((9,12))
	G[:3,:3] = -R_body_to_nav
	G[:3,6:9] = -R_body_to_nav
	G[3:6,3:6] = np.identity(3)
	G[6:,9:] = np.identity(3)

	# compute noise matrix Q. Q is 12 x 12.
	Q = np.zeros((12,12))
	sigma_in = 0.0000

	# noise
	global sigma_xg, sigma_nug, sigma_xa, sigma_nua

	# initialize Q
	Q[:3,:3] = sigma_in**2*np.identity(3)
	Q[3:6,3:6] = sigma_xg**2*np.identity(3)
	Q[6:9,6:9] = sigma_nug**2*np.identity(3)
	Q[9:12,9:12] = sigma_xa**2*np.identity(3)

	# compute Qdk. Qdk is 9 x 9
	Qdk = G.dot(Q).dot(G.T)*dt

	# get previous covariance
	global cov

	# update covariance (9x9)
	cov = Phi.dot(cov).dot(Phi.T)+Qdk

	# predict gravity in navigation frame and store prediction in global variable. used in filter.
	# subtract out 
	global g_pred
	g_pred = R_body_to_nav.dot(x_a - a)

	# periodically normalize the quaternion
	global counter
	if counter == 1000:
		state[0] = state[0]/np.linalg.norm(state[0])
		counter = 0

	# if 50 accelerometer readings were close enough to the gravity vector, robot is stationary
	# -> do measurement update
	global accel_counter

	# consider doing some sort of averaging rather than the latest g_pred, like in the initialization
	if np.abs(np.linalg.norm(a) - 9.8021) < 0.05: # if your imu is still
		accel_counter += 1
	else:
		accel_counter = 0
	if accel_counter == 200:
		print("Measurement Update")
		measurement_update()
		accel_counter = 0

	# increment counter
	counter += 1


# call initalization node
def initalize_ahrs_client():

	# wait for service to become active
	rospy.wait_for_service('initalize_ahrs')

	try:
		# create service method
		initalize_ahrs = rospy.ServiceProxy('initalize_ahrs', initRequest)

		# request the service
		print("ahrs_ekf: Starting initialization.")
		resp = initalize_ahrs()
		print("ahrs_ekf: Received initalization data.")

		# access service response
		b = resp.init_orientation
		gyro_bias = resp.gyro_bias

		# state: [b, x_g, x_a], [quaternion, gyro bias, accel bias]
		# quaternion is stored as [w, x, y, z]
		b0 = np.array([b.w,b.x,b.y,b.z])

		# numpy gyro and accel bias
		x_g = np.array([gyro_bias[0],gyro_bias[1],gyro_bias[2]])
		x_a = np.zeros(3)

		# store in global variable numpy state:
		global state
		state = np.array([b0,x_g, x_a])

		# initialize covariance
		global cov
		cov = np.identity(9)

		# imu hz
		imu_hz = rospy.get_param('imu_hz', 200)
		global dt
		dt = 1.0/imu_hz

		# TODO: make these noise terms in some launch file or something
		global sigma_xg, sigma_nug, sigma_xa, sigma_nua

		sigma_xg = 0.00000290 # Gyro (rate) random walk
		sigma_nug = 0.00068585   # rad/s/rt_Hz, Gyro white noise
		sigma_xa = 0.00001483 # Accel (rate) random walk m/s3 1/sqrt(Hz)
		sigma_nua = 0.00220313 # accel white noise
		T = 1000.0/200 # number of measurements over rate of IMU
		g = 9.8 # gravity m/s/s

		cov[:3,:3] = np.diag([sigma_nua/g,sigma_nua/g,0])**2/T
		cov[3:6,3:6] = np.identity(3)*sigma_nug**2/T
		cov[0:2,7:9] = np.diag([math.sqrt(cov[0,0]*cov[7,7]),math.sqrt(cov[1,1]*cov[8,8])]) # absolutely no idea
		cov[6:,0:3] = np.transpose(cov[0:3,6:])

	except rospy.ServiceException, e:
		print "Service call failed: %s"%e


if __name__ == "__main__":
	rospy.init_node('ahrs_ekf')

	# initialize orientation
	print('Initalizing IMU...')
	initalize_ahrs_client()

	# start callback
	rospy.Subscriber("/imu/data", Imu, imu_callback)

	rospy.spin()
