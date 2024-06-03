# IMU-Noise-Analysis-and-Performance-Evaluation

Overview
This repository contains the code and data analysis for the project "IMU Noise Analysis and Performance Evaluation." The project aims to characterize the noise and performance of an Inertial Measurement Unit (IMU) using Allan Variance analysis and other statistical measures.

Table of Contents
Introduction
Allan Variance Analysis
Location A Analysis
Short Stationary Data Analysis
Moving Data Analysis
Results and Discussion
Noise Characteristics
Performance Comparison
Environmental Noise Sources
Conclusion
References
Introduction
This project is part of the EECE5554 – Robot Sensing and Navigation course. The primary objective is to perform a detailed analysis of the noise characteristics and performance of an IMU.

Allan Variance Analysis
Location A Analysis
The Allan Variance analysis for the 5-hour data collected at Location A shows that the gyro is dominated by white noise at low averaging times and flicker noise at higher averaging times. The X-axis shows the most stability, followed by the Z and Y axes.

Short Stationary Data Analysis
Data collected in a quiet classroom shows minimal noise in the angular velocity, linear acceleration, and magnetic field measurements, indicating a stable environment.

Moving Data Analysis
Data collected while the IMU was moved in multiple directions shows spikes correlating with the movements, providing insights into the sensor's response to dynamic conditions.

Results and Discussion
Noise Characteristics
The stationary data exhibits minimal noise, likely due to thermal, mechanical, or electrical sources.
Allan Variance plots help in identifying different types of noise, such as white noise and flicker noise.
Performance Comparison
The IMU performs accurately with moderate levels of error when compared to the VN100 datasheet.
Environmental Noise Sources
Mechanical, electrical, and environmental noises are potential sources that affect the measurements.
Systematic errors are consistent while unsystematic errors are random and unpredictable.
Conclusion
The IMU noise and performance evaluation project provides a comprehensive analysis of the sensor's behavior under different conditions. The Allan Variance analysis is crucial in understanding the noise characteristics and improving sensor performance.

References
Complete Video of Data Collection
VN100 Datasheet
