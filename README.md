# A Robust Data-Driven Identification Algorithm for Nonlinear Dynamical Systems with Time Delays

This paper proposes a data-driven technique to identify nonlinear dynamical systems with
delays and estimate the time delays in the controllers' feedback. The technique relies on
only time serious record of input and output data and the system structure is unknown
while the controller is pre-designed with unknown time delays. Recently developed sparse
optimization algorithms are attracting a great attention in system identication community;
however, it's not tried for delayed systems. We build a computationally ecient and
robust sparse regression algorithm and focus on practical challenges for real world data.
We employ cross-validation techniques from machine learning for automatic model selection
and an algebraic operation for preprocessing signals. The algebraic operation minimizes
the requirement of derivatives estimation in identication procedure without any introducing
new parameters to be identied. It acts as a ltering noise and is independent from
initial conditions. To increase the accuracy of estimation, we join the sparse regression
with a bootstrapping resampling technique for real data with high level of noise. We use
Taylor expansion to parameterize the time delays in the system equations. A nonlinear
Dung oscillator is simulated to motivate the development of the proposed estimation
technique. Experimental responses of time histories of a nonlinear rotary 
exible joint are
presented to validate the eciency and accuracy of the proposed method. The studies
reported herein indicate that the proposed method is quite promising for high-dimensional
dynamical nonlinear systems with multiple inputs and multiple outputs.
Keywords: Time delay estimation, Sparse regression, Algebraic data processing, System
identication, Bootstrapping resampling, Nonlinear dynamics
