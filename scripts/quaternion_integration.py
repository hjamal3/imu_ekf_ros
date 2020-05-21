#!/usr/bin/env python
# HJ: 5/19/20
# ROS 
import rospy
from sensor_msgs.msg import Imu
from geometry_msgs.msg import Vector3Stamped
from imu_ahrs.srv import *

# other
import numpy as np
import math
from pyquaternion import Quaternion
from scipy.linalg import expm
from helper import to_skew

#debugging
from geometry_msgs.msg import PoseStamped
from geometry_msgs.msg import Quaternion as quatmsg
counter = 0
import tf


def measurement_update_with_mag():

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

	# measurement matrix magnetometer
	# TODO. What is the true north??
	Hm = np.array([0,0,1,0,0,0,0,0,0])

	# measurement matrices
	H = np.vstack((Ha,Hm))

	# previous covariance: 9x9
	global cov

	# compute individual sensor variances
	sigma_m = 0.00001
	sigma_a = 0.00001

	# noise matrix for measurement
	R = np.array([[sigma_a**2,0,0,0],[0,sigma_a**2,0,0],[0,0,sigma_a**2,0],[0,0,0,sigma_m**2]])
	# operate on R!!!! transform to proper frame!


	# compute Kalman gain
	K = cov.dot(H.T).dot(np.linalg.inv(R+H.dot(cov).dot(H.T)))

	# access current sensor readings and stack them
	global g_pred, m_pred

	try: # sometimes in first round the mag data comes before the imu data. if that's the case ignore.
		y_pred = np.concatenate((g_pred,m_pred[1]))
	except: 
		return

	# stack known gravity and magnetometer vector
	g = np.array([0,0,9.81])
	m = np.array([0]) 
	y = np.concatenate((g,m))

	# residual z
	z = y-y_pred

	# compute state correction
	dx = K.dot(z)

	# rotation update 
	p = dx[:3]
	P = to_skew(p)
	R_nav_to_body = R_body_to_nav.T
	R_nav_to_body_next = R_nav_to_body.dot(np.identity(3) - P) # (10.67)
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

def mag_callback(data):

	# get mag data
	mag = np.array([data.vector.x,data.vector.y,data.vector.z])

	# normalize magnetometer reading
	mag = mag/np.linalg.norm(mag)

	# get rotation matrix from b (body to navigation frame)
	global state
	b = state[0]
	b = Quaternion(array=b)
	R_body_to_nav = np.transpose(b.rotation_matrix)

	# predict magnetometer in navigation frame and store prediction in global variable. used in filter
	global m_pred
	m_pred = R_body_to_nav.dot(mag)

	# run measurement update once you get magnetometer data since it's slower
	#measurement_update_with_mag()


# update orientation here. use RK3
def imu_callback(data):

	# sample time. watch out for truncation
	dt = 1.0/200

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

	# visualize propagated quaternion
	if True:
		b_vis = b_next.inverse
		p = PoseStamped()
		p.header.frame_id = 'quat'
		pub = rospy.Publisher('pose', PoseStamped, queue_size=1)
		pub.publish(p)
		br = tf.TransformBroadcaster()
		# world to odom transformation. Ros uses the transformation between frames not points. The transformation
		# that takes the world frame to odom frame. 
		br.sendTransform((0, 0, 0),
		(b_vis[1],b_vis[2],b_vis[3],b_vis[0]),
		rospy.Time.now(),
		"quat",
		"world")

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
		# TODO: This might take a lot of time, so find a faster approximation
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
	Pxg       = 2e-6;        # rad^2/s^2, ss bias cov
	sigma_xg  = math.sqrt(2*lambda_g*Pxg)      # rad/s/s/rt_Hz, bias drift rate
	sigma_nug = 2.2e-3   # rad/s/rt_Hz, angle drift rate (gyro meas. noise)
	Pxa       = 2e-4;        #m^2/s^4, ss bias cov
	sigma_xa  = math.sqrt(2*lambda_a*Pxa) # m/s/s/s/rt_Hz, bias drift rate

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
	global g_pred
	g_pred = R_body_to_nav.dot(x_a - a)

	# periodically normalize the quaternion
	global counter
	if counter == 200:
		state[0] = state[0]/np.linalg.norm(state[0])

	# ... print the difference from unit norm and see


	# a check to see if filter going at 200 Hz. 
	# global counter
	# # if counter == 200:
	# # 	print('Check')
	# # 	counter = 0


	counter += 1


# call initalization node
def initalize_ahrs_client():

	# wait for service to become active
	rospy.wait_for_service('initalize_ahrs')

	try:
		# create service method
		initalize_ahrs = rospy.ServiceProxy('initalize_ahrs', initRequest)

		# request the service
		print("quat_integration: Starting initialization.")
		resp = initalize_ahrs()
		print("quat_integration: Recieved initalization data.")

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


	except rospy.ServiceException, e:
		print "Service call failed: %s"%e


if __name__ == "__main__":
	rospy.init_node('quat_integration')

	# initialize orientation
	print('Initalizing IMU...')
	initalize_ahrs_client()

	# start callback
	rospy.Subscriber("/imu/data", Imu, imu_callback)
	rospy.Subscriber("/imu/mag", Vector3Stamped, mag_callback)

	# TODO: add any parameters here like IMU hz and what not.
	# create a list of configurable parameters and maybe add a config file

	rospy.spin()
