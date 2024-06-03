import statistics
from bagpy import bagreader
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import norm


#Stationary Data
path = "k:\Stationary_data.bag"

bag = bagreader(path)
data = bag.message_by_topic('/imu')
readings = pd.read_csv(data)

w = np.array(readings['imu.orientation.w']) * (np.pi/180)
x = np.array(readings['imu.orientation.x']) * (np.pi/180)
y = np.array(readings['imu.orientation.y']) * (np.pi/180)
z = np.array(readings['imu.orientation.z']) * (np.pi/180)

#def euler_from_quaternion(x, y, z, w):
t0 = +2.0 * (w * x + y * z)
t1 = +1.0 - 2.0 * (x * x + y *y)
roll_x = np.degrees(np.arctan2(t0, t1))

t2 = +2.0 * (w * y - z * x)
t2 = np.where(t2>+1.0, +1.0,t2)
t2 = np.where(t2<-1.0, -1.0,t2)
pitch_y = np.degrees(np.arcsin(t2))

t3 = +2.0 * (w * z + x * y)
t4 = +1.0 - 2.0 * (y * y+ z * z)
yaw_z = np.degrees(np.arctan2(t3, t4))

time = np.array(readings['Time'])

#angular velocity
f, ax = plt.subplots(2, 3)
f.subplots_adjust(hspace=0.6)

ax[0,0].plot(time, np.array(readings['imu.angular_velocity.x']))
ax[0,1].plot(time, np.array(readings['imu.angular_velocity.y']))
ax[0,2].plot(time, np.array(readings['imu.angular_velocity.z']))

ax[1,0].hist(readings['imu.angular_velocity.x'], bins= 40, density=True)
ax[1,1].hist(readings['imu.angular_velocity.y'], bins= 40, density=True)
ax[1,2].hist(readings['imu.angular_velocity.z'], bins= 40, density=True)

mu_x, std_x = norm.fit(np.array(readings['imu.angular_velocity.x']))
mu_y, std_y = norm.fit(np.array(readings['imu.angular_velocity.y']))
mu_z, std_z = norm.fit(np.array(readings['imu.angular_velocity.z']))

p_x = norm.pdf(40, mu_x, std_x)
p_y = norm.pdf(40, mu_y, std_y)
p_z = norm.pdf(40, mu_z, std_z)

ax[0,0].set_xlabel('Time (Seconds)')
ax[0,0].set_ylabel('Angular Velocity_X (rad/sec)')
ax[0,0].set_title('Time vs Angular Velocity_X')
ax[0,1].set_xlabel('Time (Seconds)')
ax[0,1].set_ylabel('Angular Velocity_Y (rad/sec)')
ax[0,1].set_title('Time vs Angular Velocity_Y')
ax[0,2].set_xlabel('Time (Seconds)')
ax[0,2].set_ylabel('Angular Velocity_Z (rad/sec)')
ax[0,2].set_title('Time vs Angular Velocity_Z')

ax[1,0].set_xlabel('Angular Velocity_X (rad/sec)')
ax[1,0].set_ylabel('Frequency')
ax[1,0].set_title('Angular Velocity_X (rad/sec) vs Frequency')
ax[1,1].set_xlabel('Angular Velocity_Y (rad/sec)')
ax[1,1].set_ylabel('Frequency')
ax[1,1].set_title('Angular Velocity_Y (rad/sec) vs Frequency')
ax[1,2].set_xlabel('Angular Velocity_Z (rad/sec)')
ax[1,2].set_ylabel('Frequency')
ax[1,2].set_title('Angular Velocity_Z (rad/sec) vs Frequency')

#linear acc
f, ax = plt.subplots(2, 3)
f.subplots_adjust(hspace=0.6)

ax[0,0].plot(time, np.array(readings['imu.linear_acceleration.x']))
ax[0,1].plot(time, np.array(readings['imu.linear_acceleration.y']))
ax[0,2].plot(time, np.array(readings['imu.linear_acceleration.z']))

ax[1,0].hist(readings['imu.linear_acceleration.x'], bins= 40, density=True)
ax[1,1].hist(readings['imu.linear_acceleration.y'], bins= 40, density=True)
ax[1,2].hist(readings['imu.linear_acceleration.z'], bins= 40, density=True)

mu_x, std_x = norm.fit(np.array(readings['imu.linear_acceleration.x']))
mu_y, std_y = norm.fit(np.array(readings['imu.linear_acceleration.y']))
mu_z, std_z = norm.fit(np.array(readings['imu.linear_acceleration.z']))

p_x = norm.pdf(40, mu_x, std_x)
p_y = norm.pdf(40, mu_y, std_y)
p_z = norm.pdf(40, mu_z, std_z)

ax[0,0].set_xlabel('Time (Seconds)')
ax[0,0].set_ylabel('Linear Acceleration_X (m/s\u00b2)')  
ax[0,0].set_title('Time vs Linear Acceleration_X')
ax[0,1].set_xlabel('Time (Seconds)')
ax[0,1].set_ylabel('Linear Acceleration_Y (m/s\u00b2)')
ax[0,1].set_title('Time vs Linear Acceleration_Y')
ax[0,2].set_xlabel('Time (Seconds)')
ax[0,2].set_ylabel('Linear Acceleration_Z (m/s\u00b2)')
ax[0,2].set_title('Time vs Linear Acceleration_Z')

ax[1,0].set_xlabel('Linear Acceleration_X (m/s\u00b2))')
ax[1,0].set_ylabel('Frequency')
ax[1,0].set_title('Linear Acceleration_X (m/s\u00b2) vs Frequency')
ax[1,1].set_xlabel('Linear Acceleration_Y (m/s\u00b2))')
ax[1,1].set_ylabel('Frequency')
ax[1,1].set_title('Linear Acceleration_Y (m/s\u00b2) vs Frequency')
ax[1,2].set_xlabel('Linear Acceleration_Z (m/s\u00b2)')
ax[1,2].set_ylabel('Frequency')
ax[1,2].set_title('Linear Acceleration_Z (m/s\u00b2) vs Frequency')


#magnetic field
f, ax = plt.subplots(2, 3)
f.subplots_adjust(hspace=0.6)

ax[0,0].plot(time, np.array(readings['mag_field.magnetic_field.x']))
ax[0,1].plot(time, np.array(readings['mag_field.magnetic_field.y']))
ax[0,2].plot(time, np.array(readings['mag_field.magnetic_field.z']))

ax[1,0].hist(readings['mag_field.magnetic_field.x'], bins= 40, density=True)
ax[1,1].hist(readings['mag_field.magnetic_field.y'], bins= 40, density=True)
ax[1,2].hist(readings['mag_field.magnetic_field.z'], bins= 40, density=True)

mu_x, std_x = norm.fit(np.array(readings['mag_field.magnetic_field.x']))
mu_y, std_y = norm.fit(np.array(readings['mag_field.magnetic_field.y']))
mu_z, std_z = norm.fit(np.array(readings['mag_field.magnetic_field.z']))

p_x = norm.pdf(40, mu_x, std_x)
p_y = norm.pdf(40, mu_y, std_y)
p_z = norm.pdf(40, mu_z, std_z)

ax[0,0].set_xlabel('Time (Seconds)')
ax[0,0].set_ylabel('Magnetic Field ')  
ax[0,0].set_title('Time vs Magnetic Field in X')
ax[0,1].set_xlabel('Time (Seconds)')
ax[0,1].set_ylabel('Magnetic Field')
ax[0,1].set_title('Time vs Magnetic Field in Y')
ax[0,2].set_xlabel('Time (Seconds)')
ax[0,2].set_ylabel('Magnetic Field)')
ax[0,2].set_title('Time vs Magnetic Field in Z')

ax[1,0].set_xlabel('Magnetic Field in X')
ax[1,0].set_ylabel('Frequency')
ax[1,0].set_title('Magnetic Field in X vs Frequency')
ax[1,1].set_xlabel('Magnetic Field in Y')
ax[1,1].set_ylabel('Frequency')
ax[1,1].set_title('Magnetic Field in Y vs Frequency')
ax[1,2].set_xlabel('Magnetic Field in Z')
ax[1,2].set_ylabel('Frequency')
ax[1,2].set_title('Magnetic Field in Z vs Frequency')


#roll, pitch, yaw
f, ax = plt.subplots(2, 3)
f.subplots_adjust(hspace=0.6)

ax[0,0].plot(time, roll_x, label = 'Time VS roll_x')
ax[0,1].plot(time, pitch_y, label = 'Time VS pitch_y')
ax[0,2].plot(time, yaw_z, label = 'Time VS yaw_z')

ax[1, 0].hist(roll_x, bins=40)
ax[1, 1].hist(pitch_y, bins=40)
ax[1, 2].hist(yaw_z, bins=40)

ax[0,0].set_xlabel('Time (Seconds)')
ax[0,0].set_ylabel('roll_x (degrees)')
ax[0,0].set_title('Time vs roll_x')
ax[0,1].set_xlabel('Time (Seconds)')
ax[0,1].set_ylabel('pitch_y (degrees)')
ax[0,1].set_title('Time vs pitch_y')
ax[0,2].set_xlabel('Time (Seconds)')
ax[0,2].set_ylabel('yaw_z (degrees)')
ax[0,2].set_title('Time vs yaw_z')

ax[1, 0].set_xlabel('roll_x (degrees)')
ax[1, 0].set_ylabel('Frequency')
ax[1, 0].set_title('roll_x vs Frequency')

ax[1, 1].set_xlabel('pitch_y (degrees)')
ax[1, 1].set_ylabel('Frequency')
ax[1, 1].set_title('pitch_y vs Frequency')

ax[1, 2].set_xlabel('yaw_z')
ax[1, 2].set_ylabel('Frequency')
ax[1, 2].set_title('yaw_z vs Frequency')

plt.show()


#MEAN CALCULATION OF RPY
print('Mean & Standard Deviation of RPY:')
print('mean = ',statistics.mean(roll_x))
print('mean = ',statistics.mean(pitch_y))
print('mean = ',statistics.mean(yaw_z))
print('standard deviation = ',statistics.stdev(roll_x))
print('standard deviation = ',statistics.stdev(pitch_y))
print('standard deviation = ',statistics.stdev(yaw_z))

#MEAN CALCULATION OF ANGULAR VELOCITY
print('\nMean & Standard Deviation of Angular Velocity:')
for i in ['imu.angular_velocity.x', 'imu.angular_velocity.y', 'imu.angular_velocity.z']:
    print('mean = ',np.array(readings[i]).mean())
    print('standard deviation = ',np.array(readings[i]).std())

#MEAN CALCULATION OF LINEAR ACCELERATION
print('\nMean & Standard Deviation of Linear Acceleration:')
for i in ['imu.linear_acceleration.x','imu.linear_acceleration.x','imu.linear_acceleration.x']:
    print('mean = ',np.array(readings[i]).mean())
    print('standard deviation = ',np.array(readings[i]).std())
    print('Median = ',np.median(np.array(readings[i])))

#MEAN CALCULATION OF MAGNETIC FIELD
print('\nMean & Standard Deviation of Magnetic Field:')
for i in ['mag_field.magnetic_field.x', 'mag_field.magnetic_field.y', 'mag_field.magnetic_field.z']:
    print('mean = ',np.array(readings[i]).mean())
    print('standard deviation = ',np.array(readings[i]).std())