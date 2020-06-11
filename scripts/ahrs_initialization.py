#!/usr/bin/env python
# HJ: 5/18/20
import rospy
from sensor_msgs.msg import Imu
from geometry_msgs.msg import Vector3Stamped, Quaternion
from imu_ahrs.srv import *
import math
import numpy as np
from tf.transformations import quaternion_from_euler
import tf

n_imu = 0
sum_accel = np.zeros(3)
sum_gyro = np.zeros(3)

def imu_callback(data):

	# get accel data
	a = data.linear_acceleration

	# get gyro data: angular velocity
	w = data.angular_velocity

	global n_imu, sum_accel, sum_gyro, num_data

	if n_imu < num_data:

		# add latest accel and gyro data
		sum_accel -= [a.x, a.y, a.z]

		sum_gyro += [w.x, w.y, w.z]

		n_imu += 1
	else:

		rospy.Service('initalize_ahrs', initRequest, handle_init_ahrs)

		# shut down subscription to imu data
		global sub_imu

		sub_imu.unregister()


def handle_init_ahrs(req):
# def computeVals():

	global sum_accel, sum_gyro, num_data

	# estimate gravity vector 
	g_b = sum_accel/num_data

	# compute initial roll (phi) and pitch (theta)
	phi = math.atan2(-g_b[1],-g_b[2])
	theta = math.atan2(g_b[0], math.sqrt(g_b[1]**2 + g_b[2]**2))

	# set yaw (psi) as zero since no absolute heading available
	psi = 0

	# q is navigation to body transformation: R_bi
	# YPR: R_ib = R(yaw)R(pitch)R(Roll)
	# RPY: R_bi = R(-Roll)R(-Pitch)R(-yaw)
	q = quaternion_from_euler(-phi, -theta, -psi,'sxyz')

	quat = Quaternion()
	quat.x = q[0]
	quat.y = q[1]
	quat.z = q[2]
	quat.w = q[3] # tf convention

	# compute gyroscope biases
	gyro_avg = sum_gyro/num_data
	gyro_biases = [gyro_avg[0], gyro_avg[1], gyro_avg[2]]

	print("ahrs_initialization_server: Initial YPR: ", psi, theta, phi)
	print("ahrs_initialization_server: Gyro Biases ", gyro_biases)
	return initRequestResponse(quat, gyro_biases)


if __name__ == "__main__":

	# TODO: close node upon calibration
	rospy.init_node('ahrs_initialization_server')

	# get number of points to average
	global num_data
	num_data = rospy.get_param('num_data', 1000)

	# Imu callback
	global sub_imu
	sub_imu = rospy.Subscriber("/imu/data", Imu, imu_callback)

	rospy.spin()
