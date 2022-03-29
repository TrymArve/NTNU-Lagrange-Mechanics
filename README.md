# Lagrange Mechanics (NTNU)
MATLAB support for effectively modelling your robotic system, testing your controllers, animating the simulations and saving it as videos.

Includes:
- "MakeLagrange(...)" - function that takes the positions of every mass described in terms of the generalized coordinates, then returns the kinetic and potential energies associated with every mass, and the system in total. Though, most importantly, it returns the Euler-Lagrange equations of motion. (see: "Inertia Wheel Pendulum" for an example of how to include inertia. Automatic inertia functionality will be added soon)

- A general dynamics integrator:

    "controller" - simply configure the controller struct to design your controller. (f.ex: controller.type = "PD"; controller.PD.Kp = 10; controller.PD.Kd = 1;)

    "parameters" - define some parameters such as length of arms, gravitational field, or any other necessary parameters of your system. (parameters = [g; L1; L2];)

    "mass"       - define the sizes of the masses in your system. (mass = [m1; m2];)

- "Simulate_EL;" - Script that automatically simulates your system, and plots the result. It then saves the data in "tsim", "xsim", "usim", "ddqsim", "ref".

- "Animate(...)" - Function that takes in the some object that you define, and animates them. Use this to see how your system moves yourself! This function can also save the animation as a video.

- Some other functions that are relevant for robotics.
