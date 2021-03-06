# imu_ekf_ros
A C++ and python ROS package that fuses the accelerometer and gyroscope of an IMU to estimate attitude.  
<p align="center">
	  <img src="/results/screencap.png" />
</p>
After catkin_make and compiling the scripts, cd into the launch folder and type:
roslaunch cpp_ekf.launch for the C++ version (better and more up to date).
roslaunch ekf.launch for the Python version (probably broken).

-This node subscribes to sensor_msgs/IMU messages on the topic '/imu/data', containing accelerometer and gyroscope data.  

-This node publishes a quaternion to 'AHRS_EKF_quat'.  

-To visualize, run rviz, create an Axes and change the reference frame to 'ENU'. Create another Axes and change the reference frame to 'IMU'. 

Wait for 5 seconds for the initialization procedure. Update the following noise terms inside code for your IMU:

		sigma_xg # Gyro (rate) random walk  
		sigma_nug # Gyro white noise  
		sigma_xa # Accel (rate) random walk   
		sigma_nua # Accel white noise  

Tested with Xsens MTI-20 and Sensonor STIM300.
Primary reference is 'Aided Navigation: GPS with High Rate Sensors' by Jay A. Farrell, chapter 10.
