function[dstate] = SystemDynamics_Integrator_2(t, state, mass, parameters, controller, values)

% t             - Time variable
% state         - Contains all state vaiables: [pos1;pos2;po3,...,vel1,vel2,vel3,...]
% mass          - Contains values of the various masses and inertia: ex: [m1,I1,m2,I2,...]
% parameters    - Contains the vaious parameters in the system: ex:
%                 [g,L1,L2,](the first value should be the gravitational field strength)
% controller    - A struct containing the fields:

%                  -"type" what controller are you using (required)
%                  -"ref","dref","ddref" reference state(and derivatives) (required for some controllers)
%                  -"FL" - set to a string: "on" to turn on Feeback Linearization
%                  -"transformation" to transform u --> u_transormed before acting on system
%                  -"B"  input signals - to - force on system (Q = B*u) (NOTE: usim contains the forces on the system, after B)

%                  - type-specific fields: 
%                     - types: "off","PD","KF","EXT" - selects type of controller

%                           - off: sets input to u = 0 (of approriate size)

%                           - PD : Uses a PD-controller, has sub-fields:
%                                    - "Kp": proportional gain                  (required)
%                                    - "Kd": derivative gain                    (required)
%                                    - "C" : measurement matrix - y = C*q      (required)
%                                            (e = C*(ref - q))

%                           - KF : Uses standard feedback structure: u = -K*q + F*r, uses sub-fields:
%                                    - "K" : state feedback gain-matrix (function handle or constant)        (required)
%                                    - "F" : reference feedforward gain-matrix (function handle or constant) (not required)

%                           - EXT: Uses inputs given by either: "controller.EXT.function" or "controller.EXT.timeseries"
%                                    - EXT.source: can be either "function" or "timeseries"  (required)
%                                           -Timeseries uses linear interpolation to estimate input value
%                                    - Example 1: 
%                                           controller.type = "EXT"
%                                           controller.EXT.source = "function"
%                                           controller.EXT.function = @(t) [   sin(t)   ; 
%                                                                           sin(t)/(t+1)]
%                                    - Example 2:
%                                           controller.type = "EXT"
%                                           controller.EXT.source = "timeseries"
%                                           controller.EXT.function = [ (1 x nt)array of timestamps;
%                                                                      (nu x nt)array of corresponding inputvalues]

%                     - ref   : a function of time describing the reference point as a function handle: ref(t) 
%                     - dref  : the derivative of ref  as a function handle
%                     - ddref : the derivative of dref as a function handle

% values        - A struct containing needed values and additional info/configuration:
%                     - "nq"     : number of states(length of q) (required)
%                     - "umax"   : the maximum values(absolute) the inputs signal can have. EX: values.umax = [200;150;50;0];
%                     - "ExtMCG" : uses function handles given via "values.M/C/G" as dynamics rather than the default MCG matlabfunction made by the "MakeLagrange" function
%                     - "M"      : function handle for Mass-matrix
%                     - "C"      : function handle for Coreolis-matrix
%                     - "G"      : function handle for Gravity-vector


q  = state(1:values.nq);
dq = state(values.nq+1:end);

%external MCG:
if isfield(values,'ExtMCG') && (values.ExtMCG == 1)
    M = values.M(q,dq);
    C = values.C(q,dq);
    G = values.G(q,dq);
else
    [M,C,G] = MCG(q,dq,mass,parameters);
end


%%%%% CONTROL:
if controller.type == "off"
        u = zeros(values.nq,1);




% PD:
elseif (controller.type == "PD")
    if ~isfield(controller,'Kdd')
        u = controller.PD.Kp*controller.PD.C*(controller.ref(t) - q) + controller.PD.Kd*controller.PD.C*(controller.dref(t) - dq);
    else
        u = controller.PD.C*controller.ddref(t) + controller.PD.Kp*controller.PD.C*(controller.ref(t) - q) + controller.PD.Kd*controller.PD.C*(controller.dref(t) - dq);
    end




% KF:
elseif (controller.type ==  "KF")
        if isa(controller.KF.K,'function_handle')
            u = -controller.KF.K(t,q,dq)*[q;dq];
        else
            u = -controller.KF.K*[q;dq];
        end
        if isfield(controller.KF,'F')
            if isa(controller.KF.F,'function_handle')
               u = u + controller.KF.F(t,q,dq)*controller.ref(t); 
            else
               u = u + controller.KF.F*controller.ref(t);
            end
        end





% External controller:
elseif controller.type ==  "EXT"
        if controller.EXT.source == "function"
            u = controller.EXT.function(t,q,dq);
        elseif controller.EXT.source == "timeseries"
            u = interp1(controller.EXT.timeseries(1,:),controller.EXT.timeseries(2:end,:),t);
        else
            error('ERROR: select valid source for External controller: "function" or "timeseries".')
        end






% MPC controller:  
elseif controller.type ==  "MPC"
        try
            t_prev = evalin('base', 'usim_mpc(end,1);');
        catch
        end
            
        if  (t==0) || (t-t_prev) >= controller.MPC.delta_u
            if t ~= 0
                % Set the initial guess equal to the previous solution:
                controller.MPC.q_guess      = evalin("base",'Reserved_q_opt;');
                controller.MPC.u_guess      = evalin("base",'Reserved_u_opt;');
                controller.MPC.q_coll_guess = evalin("base",'Reserved_q_coll_opt;');
                controller.MPC              = CLC_guess(controller.MPC);
                
                % Force the next optimal input to be close to the previous
                % second input:
                u = evalin('base', "Reserved_u_next;");
                controller.MPC = CLC_limit_initial_u(controller.MPC,u,controller.MPC.u_deviation);
            end
                 
            %%%%%%% Evaluate optimal input:
            controller.MPC.initial_state = [q;dq];
            controller.MPC = CLC_solve(controller.MPC);
            u = controller.MPC.u;
            u_next = controller.MPC.u_opt(2,:);
            assignin("base","Reserved_u_next",u_next);
            %%%%%%% Save the optimal solution:
            assignin("base","Reserved_q_opt",controller.MPC.q_opt);
            assignin("base","Reserved_u_opt",controller.MPC.u_opt);
            assignin("base","Reserved_q_coll_opt",controller.MPC.q_coll_opt);
            %%%%%%% Store first input to base work space:
            t_str = num2str(t);
            u_str = num2str(u(1));
            for ii = 2:length(u)
                u_str = [u_str,',', num2str(u(ii))];
            end
            if t == 0
                % usim_mpc:
                text = ['usim_mpc = [', t_str, ',' , u_str,'];'];
                evalin('base',text)
            else
                % usim_mpc:
                text = ['usim_mpc = cat(1,usim_mpc,[', t_str, ',' , u_str,']);'];
                evalin('base',text)
            end       
            if isfield(controller.MPC,'printMPCsteps') && (controller.MPC.printMPCsteps == 1)
                PlotNLP(controller.MPC)
                text = ['length(usim_mpc(:,1))'];
                block_num = evalin('base',text);
                title(['OL trajectory ', num2str(block_num)])
            end

        else
            u = evalin('base', "usim_mpc(end,2:end)';");
        end
        

        
        
else
        error('ERROR: Choose valid controller type')
end




%%%%%% Transformations:


if isfield(controller,'transformation') && (controller.transformation.enable == 1)
    u = controller.transformation.function(u,q,dq,t);
end


if isfield(controller,'B') && ~(controller.type == "off")
    Q = controller.B*u;
else
    Q = u;
end



if (isfield(controller,'FL')) && (controller.FL == "on")
    Q = (M*Q + C*dq + G);
    Q(not(any(controller.B,2))) = 0;
end


if isfield(values,'umax') && ~isstring(values.umax) && ~(controller.type == "off")
    if length(values.umax) > 1
        Qmax = controller.B*values.umax;
        for j = 1:length(Q)
            if abs(Q(j)) > Qmax(j)
                Q(j) = (1 - 2*(Q(j) < 0))*Qmax(j);
            end
        end
    else
        for j = 1:length(Q)
            if abs(Q(j)) > values.umax
                Q(j) = (1 - 2*(Q(j) < 0))*values.umax;
            end
        end
    end
end


%%%%%% Integrate:
ddq = M\(Q - C*dq - G);
dstate = [dq;ddq];




%%%%%% save stuff to workspace:
if isfield(controller,'B') && ~(controller.type == "off")
    u = Q(any(controller.B,2));  %don'save the forced that identically zero
else
    u = Q;
end

t_str = num2str(t);

% u:
u_str = num2str(u(1));
for ii = 2:length(u)
    u_str = [u_str,',', num2str(u(ii))];
end

% ddq:
ddq_str = num2str(ddq(1));
for ii = 2:values.nq
    ddq_str = [ddq_str, ',' num2str(ddq(ii))];
end

if t == 0
    % usim:
    text = ['usim = [', t_str, ',' , u_str,'];'];
    evalin('base',text)
    % ddqsim:
    text = ['ddqsim = [', t_str, ',' , ddq_str,'];'];
    evalin('base',text)
else
    % usim:
    text = ['usim = cat(1,usim,[', t_str, ',' , u_str,']);'];
    evalin('base',text)
    % ddqsim:
    text = ['ddqsim = cat(1,ddqsim,[', t_str, ',' , ddq_str,']);'];
    evalin('base',text)
end



end






















