
close all; clc;

values.nq = length(S.q);

try
    mass;
catch
    error('SIMULATION FAILED: mass-vector is not defined. Define a vector names "mass", that contains the values of the masses and inertias in your system.(in the same order as they are defined in the system struct S)')
end
try
    parameters;
catch
    error('SIMULATION FAILED: parameter-vector is not defined. Define a vector names "parameters", that contains the values of the parameters in your system.(in the same order as they are defined in the system struct S)')
end
try
    controller;
catch
    disp('WARNING: "controller"-struct is not defined. Proceeding with controller.type = "off"')
    controller.type = "off";
end
try
    tf;
catch
    error('SIMULATION FAILED: "tf"(time finished/ duration) is not defined. (must have the name: "tf")')
end
try
    init_state;
catch
    error('SIMULATION FAILED: "init_state"(initial coordinates/velocities) is not defined. (must have the name: "init_state")')
end
try
    controller.ref(1);
catch
    disp('WARNING: "controller.ref". Proceeding with reference equal to zero')
    controller.ref  = @(t) zeros(values.nq,1)*(t-t);
    controller.dref = @(t) zeros(values.nq,1)*(t-t);
    controller.ddref = @(t) zeros(values.nq,1)*(t-t);
end


fprintf(' \n \r Simulating... ')
[tsim,xsim] = ode45(@(t,x) SystemDynamics_Integrator(t, x, mass, parameters, controller, values),[0,tf],init_state);
fprintf('done. \n \r')

% Make usim and ddqsim the same length as tsim and xsim:
ddqsim  = CleanSeries(ddqsim,tsim);
usim    = CleanSeries(usim,tsim);






disp('Resulting variables:')

text = ['tsim   - vector of timestamps.         size = nt x 1    (', num2str(size(tsim)), ')'];
disp(text)
text = ['xsim   - matrix of states(q and dq).   size = nt x 2*nq (', num2str(size(xsim)), ')'];
disp(text)
text = ['ddqsim - matrix of accelerations       size = nt x nq   (', num2str(size(ddqsim)), ')'];
disp(text)
text = ['usim   - matrix of control inputs      size = nt x nu   (', num2str(size(usim)), ')'];
disp(text)


% Plot Simulation
if controller.type == "off"
    Plot_Simulation(S,tsim,xsim);
elseif controller.type == "MPC"
    Plot_Simulation(S,tsim,xsim);
    Plot_Inputs(usim_mpc(:,1),usim_mpc(:,2:end),1);
else
    Ref = Plot_Simulation(S,tsim,xsim,controller);
    usim = Plot_Inputs(tsim,usim);
    text = ['Ref    - matrix of reference states    size = nt x nq   (', num2str(size(Ref)), ')'];
    disp(text)

end



% if controller.type == "MPC"
%     usim = [];
%     ind_old = 0;
%     for i = 1:length(usim_mpc(:,1))
%         [~,ind] = min(abs(tsim-usim_mpc(i,1)));
%         usim = [usim; ones(ind-ind_old,1)*usim_mpc(i,2:end)];
%         ind_old = ind;
%     end
% end














