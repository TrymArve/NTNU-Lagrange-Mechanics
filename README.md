# Lagrange Mechanics (NTNU)
MATLAB functions for quickly deriving the Euler-Lagrange equations for your robotic system, simulating the system while testing your controllers and animating the simulations. One can also save the animations as videos.

|❗  OBS: The files in this repository rely on my General MATLAB Library: [General Library](https://github.com/TrymArve/General-Matlab-Library) |
|---|

Includes:
- "MakeLagrange(...)" - function that takes the positions of every mass described in terms of the generalized coordinates, then returns the kinetic and potential energies associated with every mass, and the system in total. Though, most importantly, it returns the Euler-Lagrange equations of motion. (see: "Inertia Wheel Pendulum" for an example of how to include inertia. Automatic inertia functionality will be added soon)

- A general dynamics integrator:

    -- "controller" - simply configure the controller struct to design your controller. F.ex:

        controller.type  = "PD"; 
        controller.PD.Kp = 10; 
        controller.PD.Kd = 1;

    -- "parameters" - define some parameters such as length of arms, gravitational field, or any other necessary parameters of your system. F.ex: 
    
        parameters = [g; L1; L2];

    -- "mass"       - define the sizes of the masses in your system. 
    
        mass = [m1; m2];

- "Simulate_EL;" - Script that automatically simulates your system, and plots the result. It then saves the data in "tsim", "xsim", "usim", "ddqsim", "ref".

- "Animate(...)" - Function that takes in the some object that you define, and animates them. Use this to see how your system moves yourself! This function can also save the animation as a video.

- Some other functions that are relevant for robotics.
