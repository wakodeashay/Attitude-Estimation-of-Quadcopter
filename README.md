# Attitude-Estimation-of-Quadcopter

This repository contains implementation of State Estimation using Kalman Filter, Extended Kalman Filter and Unscented Kalman Filter
xdot=A(x)*x+Gamma*u+d
y=G(x)+v
1. main.m : Contains code to three filters kf, ekf and ukf
2. kG.m : function returns output for given state
3. kCmat.m : Jacobian of G(x)
4. Areal.m : function which return A(x)
5. cov.m : Calculates the covariance of estimation error
6. Attitude Estimation of Quadcopter.pdf : Details of Implementation, and graphs
