# imu_ahrs_ros
A general ROS package that fuses the accelerometer and gyroscope of an IMU to estimate orientation. 

After catkin_make and compiling the scripts, cd into the launch folder and type:
roslaunch ahrs.launch  

Your IMU publishes sensor_msgs/IMU messages on the topic '/imu/data', containing accelerometer and gyroscope data.

To visualize, run rviz, create an Axes and change the reference frame to 'NED'. Create another Axes and change the reference frame to 'quat'.

![GitHub Logo](/results/screencap.png)
Format: ![Alt Text](url)
