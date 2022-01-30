function[dstate] = SystemDynamics_Integrator(t, state, mass, parameters, controller, values)

% t             - Time variable
% state         - Contains all state vaiable: [pos1;pos2;po3,...,vel1,vel2,vel3,...]
% mass          - Contains values of the various masses and inertia: ex: [m1,I1,m2,I2,...]
% parameters    - Contains the vaious parameters in the system: ex:
%                 [g,L1,L2,](the first value should be the gravitational field strength)
% controller    - A struct containing the field "type","ref","dref","ddref" and type-specific fields: 
%                     - types: "off","KF","PD","FL","EXT" - selects type of controller
%                           - off: sets input to u = 0 (of approriate size)
%                           - KF : Uses standard feedback structure: u = -K*q + F*r, uses sub-fields: (MIMO)
%                                    - "K" : state feedback gain-matrix
%                                    - "F" : reference feedforward gain-matrix
%                           - PD : Uses a PD-controller, has sub-fields: (SISO)
%                                    - "Kp": proportional gain
%                                    - "Kd": derivative gain
%                                    - Example: 
%                                           controller.type = "PD"
%                                           controller.PD.Kp = 5
%                                           controller.PD.Kd = 2
%                           - FLKF: Uses a Feedback Linearization controller with KF controller, has sub-fields: (MIMO)
%                                    - "K": state feedback gain-matrix
%                                    - "F": reference feedforward gain-matrix
%                                    - Example: 
%                                           controller.type = "FLKF"
%                                           controller.FLKF.K = [3 2 1; 6 5 2; 1 3 2]
%                                           controller.FLKF.F = [3 2; 6 5; 1 4]
%                           - FLPD: Uses a Feedback Linearization controller with PD, has sub-fields: (SISO)
%                                    - "Kp": proportional gain
%                                    - "Kd": derivative gain
%                                    - "Mask": Mask what signals are
%                                              corrected by PD control
%                                    - Example: 
%                                           controller.type = "FLPD"
%                                           controller.FLPD.Kp = 6
%                                           controller.FLPD.Kd = 1
%                                           controller.FLPD.Mask = [1;0;0] (if used on MIMO)
%                           - EXT: Uses inputs given by either:
%                                       "controller.EXT.function" or "controller.EXT.timeseries"
%                                    - EXT.source: can be either "function" or "timeseries"
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
%                     - ref   : a function of time describing the reference
%                            point as a function handle: ref(t)
%                     - dref  : the derivative of 
%                     - save: if controller.save == 1, then the controll input is stored in a global timeseries U.
% values        - A struct containing needed values and additional info/configuration: 
%                     - "nq"     : number of states(length of q)
%                     - "umax"   : the maximum values(absolute) the inputs signal can have. EX: values.Umax = [200;150;50;0];
%                     - "ExtMCG" : uses function handles given via "values.M/C/G" as dynamics rather than the default
%                     - "M"      : function handle for Mass-matrix
%                     - "C"      : function handle for Coreolis-matrix
%                     - "G"      : function handle for Gravity-vector


q  = state(1:values.nq);
dq = state(values.nq+1:end);

if isfield(values,'ExtMCG') && (values.ExtMCG == 1)
    M = values.M(q,dq);
    C = values.C(q,dq);
    G = values.G(q,dq);
else
    [M,C,G] = MCG(q,dq,mass,parameters);
end

if isfield(controller,'saveddq') && (controller.saveddq == 1)
    global DDQ
end

%CONTROL:
switch controller.type
    case "off"
        u = zeros(values.nq,1);
    case "KF"
        u = -controller.KF.K*[q;dq];
        if isfield(controller.KF,'F')
        u = u + controller.KF.F*controller.ref(t);
        end
    case "PD"
        u = controller.PD.Kp*(q - controller.ref(t)) + controller.PD.Kd*(dq - controller.dref(t));
        u = u.*controller.PD.Mask;
    case "FLKF"
        v = controller.ddref(t) - controller.FLKF.K*([q - controller.ref(t); dq - controller.dref(t)]) - controller.FLKF.F*controller.dref(t);
        u = (M*v + C*dq + G).*controller.FLKF.Mask;
    case "FLPD"
        v = controller.ddref(t) - controller.FLPD.Kp*(q - controller.ref(t)) - controller.FLPD.Kd*(dq - controller.dref(t));
        u = (M*v + C*dq + G).*controller.FLPD.Mask;
    case "EXT"
        if controller.EXT.source == "function"
            u = controller.EXT.function(q,dq,t);
        elseif controller.EXT.source == "timeseries"
            u = interp1(controller.EXT.timeseries(1,:),controller.EXT.timeseries(2:end,:),t);
        else
            error('ERROR: select valid source for External controller')
        end
    case "SF" %(super feedback - Antons Method for stabilizing periodic trajecotries)
        NS = interp1(controller.SF.NMG(1,:)',controller.SF.NMG(2:4,:)',t); %nominal s-value
        if t > 0; AS = controller.SF.AS(q,dq,DDQ(2:end,end)); %actual [s,ds,dds] value.
        else; AS = [0;0;0]; end
        if isfield(controller.SF,'integrate_W') && controller.SF.integrate_W ~= 0
            if controller.SF.init_s(1) ~= NS(1)
                [~,intW] = ode45(@(t,x) controller.SF.W(t,x),[controller.SF.init_s(1);AS(1)],0);
                W = intW(end);
            else
                W = 0;
            end
        else
            W = controller.SF.W(AS(1),controller.SF.init_s(1));
        end
        
        if isfield(controller.SF,'integrate_R') && controller.SF.integrate_R ~= 0 
            if controller.SF.init_s(1) ~= NS(1)
                [~,intR] = ode45(@(t,x) controller.SF.R(t,x),[controller.SF.init_s(1);AS(1)],0);
                R = intR(end);
            else
                R = 0;
            end
        else
            R = controller.SF.R(AS(1),controller.SF.init_s(1));
        end
        I = AS(2)^2 - (exp(-2*W))*(controller.SF.NMG(3,1)^2 - (exp(-2*W))*2*R);
        v = -controller.SF.K(NS(1),NS(2),NS(3))*[I;controller.SF.Y(q,AS(1));controller.SF.dY(dq,AS(1),AS(2))];        
        u = controller.SF.f(v,q,dq).*controller.SF.Mask;
    otherwise
        error('ERROR: Choose valid controller type')
end

if isfield(controller,'transformation') && (controller.transformation.enable == 1)
    u = controller.transformation.function(u,q,dq,t);
end

if isfield(values,'umax') && ~isstring(values.umax) && ~(controller.type == "off")
    ind = find(abs(u) >= values.umax);
    for i = 1:length(ind)
        u(ind(i)) = (1 - 2*(u(ind(i)) < 0))*values.umax(ind(i));
    end
end

if isfield(controller,'save') && controller.save == 1
    global U
    if isempty(U)
        U = [t;u];
    else
        U = cat(2,U,[t;u]);
    end
end


ddq = M\(u - C*dq - G);

if isfield(controller,'saveddq') && (controller.saveddq == 1)
    if isempty(DDQ)
        DDQ = [t;ddq];
    else
        DDQ = cat(2,DDQ,[t;ddq]);
    end
end


dstate = [dq;ddq];

end