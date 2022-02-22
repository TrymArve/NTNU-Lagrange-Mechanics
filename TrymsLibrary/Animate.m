function[Fig] = Animate(tsim,qsim,obj,config)

%Input 'help' to return a help-struct, instead of Animating;
if (nargin == 3)
    config.dummy = [];
elseif (nargin <= 1)
    if tsim == "example"
        %Fig = Example();  %out dated
    else
        Fig = Help();        
    end
    return;
end



% rotation matrices: (for arrows)
theta = pi/4+pi-pi/12;
R1 = [cos(theta) -sin(theta);
     sin(theta) cos(theta)];
theta = -theta;
R2 = [cos(theta) -sin(theta);
     sin(theta) cos(theta)];
arrow_arm = 0.2;


%Confugitation:
if ~isfield(config,'scale')
config.scale = 1;
end

if ~isfield(config,'figurelocation')
    config.figurelocation = [1 1];
end    
if ~isfield(config,'figureheight')
    config.figureheight = 1080;
end  
if ~isfield(config,'aspect')
    config.aspect = 5/4;
end
if ~isfield(config,'figurewidth')
    config.figurewidth = config.aspect*config.figureheight;
end

if ~isfield(config,'position')
    config.position = [config.figurelocation config.figurewidth config.figureheight];
end

%%%%% Make figure:
Fig = figure('Position', config.position,'visible','off');
%%%%%

if ~isfield(config,'frameheight') && ~isfield(config,'axis')
    error('WARNING: You must define either: "config.frameheight" or "config.axis"');
end

if ~isfield(config,'framecenter')
    config.framecenter = [0 0];
end  
if ~isfield(config,'framewidth')
    config.framewidth = config.frameheight*(config.figurewidth/config.figureheight);
end
centerx = config.framecenter(1);
centery = config.framecenter(2);
width   = config.framewidth;
height  = config.frameheight;

if ~isfield(config,'axis')
        config.axis = [-width/2 width/2 -height/2 height/2] + [centerx centerx centery centery];
end


if ~isfield(config,'simspeed')
    config.simspeed = 1;
end
if ~isfield(config,'tf')
    config.tf = tsim(end);
end
if ~isfield(config,'grid')
    config.grid = 'on';
end
if ~isfield(config,'enterToStart')
    config.enterToStart = 0;
end
if ~isfield(config,'ts')
    config.ts = 0;
end
%%% Video Properties:
if isfield(config,'video')
    if ~isfield(config.video,'resolution')
        config.video.resolution = 0.01;
    end
    if ~isfield(config.video,'profile')
        config.video.profile = 'Archival';
    end
else
    config.video.enable = "off";
end

% Set Up:
axis(config.axis)
SimSpeed = config.simspeed;
tf = config.tf;
ts = config.ts;
qsim = qsim*config.scale;


if config.video.enable == "on"
    disp('-----------')
    LengthOfAnimation = (tf-ts)/(SimSpeed*config.video.resolution) + 1; % add one for buffer(does not need to be accurate)
    movieVector(ceil(LengthOfAnimation)) = getframe;
    SecondsPerFrame = 72/556; %(uncompressed AVI)
    SecondsPerFrame = 86.5/501;%(uncompressed AVI)
    SecondsPerFrame = 405.7/2000;%(uncompressed AVI)
    EstimatedProcessingTime = ceil(LengthOfAnimation*SecondsPerFrame);
    disp(['Estimated processing time: ', num2str(EstimatedProcessingTime), 'seconds'])
    LengthOfAnimation = 0;
    fprintf('Processing Video....')
    ProcessingTime = tic;
elseif config.enterToStart == 0
    fprintf('Animating... ')
end

% Animation:
t_disp = ts/SimSpeed;
tic
start = 0;
 while t_disp < tf/SimSpeed
    clf

    %Interpolate to find state at t_disp:
    q   = interp1(tsim,qsim,SimSpeed*t_disp)';
    
    %Display time on x-axis:
    time_text = num2str(SimSpeed*t_disp);
    time_dot = find(time_text == '.');
    time_text = [time_text, '000000'];
    time_text = time_text(1:time_dot+3);
    xlabel(['t = ', time_text])
    
    %Animate
    F = fieldnames(obj);
    for i = 1:length(F)
        f = F{i};
        
        type = getfield(obj,f,'type');
        if iscell(type)
            L = length(type);
        else
            L = 1;
            type = {type};
        end
        
        if (isfield(getfield(obj,f),'disable')) && (getfield(obj,f,'disable') == 1)
        else
        
        for ii = 1:L
            
            if isfield(getfield(obj,f),'color')
                col = getfield(obj,f,'color');
                if isa(col,"function_handle")
                    col = col(q,t_disp*SimSpeed);
                end
            else
                col = 'k';
            end
            
            if isfield(getfield(obj,f),'linestyle')
                linestyle = getfield(obj,f,'linestyle');
            else
                linestyle = '-';
            end
            
            switch type{ii}
                
                case 'ball'
                    
                    h = getfield(obj,f,'c');
                    c = h(q);
                    r = getfield(obj,f,'r');
                    ball(c,r,col)

                case 'line'
                    
                    h = getfield(obj,f,'a');
                    a = h(q);
                    h = getfield(obj,f,'b');
                    b = h(q);
                    line(a,b,col)

                case 'point'
                    
                    h = getfield(obj,f,'p');
                    p = h(q);
                    point(p,col)
                    
                case 'arrow'
                    
                    h = getfield(obj,f,'head');
                    head = h(q);
                    h = getfield(obj,f,'tail');
                    tail = h(q);
                    arrow(head,tail,col,R1,R2,arrow_arm);

                case 'box'
                    
                    if isfield(getfield(obj,f),'def')
                        def = getfield(obj,f,'def');
                    else
                        def = 'CC';
                        disp('NB: box definition set to CC')
                    end
                    h = getfield(obj,f,'B1'); 
                    B1 = h(q);
                    h = getfield(obj,f,'B2'); 
                    B2 = h(q);                    
                    % (B1/B2 can represent: Center/Any corner(CC), Center/Diagonal(CD), Two Opposite Corners(2C))
                    box(B1,B2,def,col,linestyle);
                    
                case 'curve'
                    
                    h = getfield(obj,f,'x');
                    x = h(q);
                    h = getfield(obj,f,'y');
                    y = h(q);
                    if isfield(getfield(obj,f),'linestyle')
                        linestyle = getfield(obj,f,'linestyle');
                    else
                        linestyle = '-';
                    end
                    curve(x,y,col,linestyle);
                    
                case 'moment'

                    h = getfield(obj,f,'target');
                    target = h(q);
                    h = getfield(obj,f,'magnitude');
                    magnitude = h(q);
                    if isfield(getfield(obj,f),'thickness')
                        h = getfield(obj,f,'thickness');
                        thickness = h(q);    
                    else
                        thickness = 1;
                    end
                    if isfield(getfield(obj,f),'color2')
                        col2 = getfield(obj,f,'color2');
                    else
                        col2 = col1;
                    end
                    if isfield(getfield(obj,f),'minr')
                        minr = getfield(obj,f,'minr');   %minimum radius
                        moment(target,magnitude,thickness,col,col2,minr)
                    else
                    moment(target,magnitude,thickness,col,col2)
                    end

                case 'timed_ball'
                    
                    h = getfield(obj,f,'c');
                    c = h(q,t_disp*SimSpeed);
                    h = getfield(obj,f,'r');
                    r = h(q,t_disp*SimSpeed);
                    ball(c,r,col)
                    
                case 'timed_line'
                    
                    h = getfield(obj,f,'a');
                    a = h(q,t_disp*SimSpeed);
                    h = getfield(obj,f,'b');
                    b = h(q,t_disp*SimSpeed);
                    line(a,b,col)
                    
                case 'timed_point'
                    
                    h = getfield(obj,f,'p');
                    p = h(q,t_disp*SimSpeed);
                    point(p,col)
                    
                case 'timed_arrow'
                    
                    h = getfield(obj,f,'head');
                    head = h(q,t_disp*SimSpeed);
                    h = getfield(obj,f,'tail');
                    tail = h(q,t_disp*SimSpeed);
                    arrow(head,tail,col,R1,R2,arrow_arm);
                    
                case 'timed_curve'
                    
                    h = getfield(obj,f,'x');
                    x = h(q,t_disp*SimSpeed);
                    h = getfield(obj,f,'y');
                    y = h(q,t_disp*SimSpeed);
                    curve(x,y,col,linestyle);
                    
                case 'timed_box'
                    
                    if isfield(getfield(obj,f),'def')
                        def = getfield(obj,f,'def');
                    else
                        def = 'CC';
                        disp('NB: timed_box definition set to CC')
                    end
                    h = getfield(obj,f,'B1'); 
                    B1 = h(q,t_disp*SimSpeed);
                    h = getfield(obj,f,'B2'); 
                    B2 = h(q,t_disp*SimSpeed);                    
                    % (B1/B2 can represent: Center/Any corner(CC), Center/Diagonal(CD), Two Opposite Corners(2C))
                    box(B1,B2,def,col,linestyle);
                    
                case 'timed_moment'

                    h = getfield(obj,f,'target');
                    target = h(qt_disp*SimSpeed);
                    h = getfield(obj,f,'magnitude');
                    magnitude = h(qt_disp*SimSpeed);
                    h = getfield(obj,f,'thickness');
                    thickness = h(qt_disp*SimSpeed);    
                    if isfield(getfield(obj,f),'color2')
                        col2 = getfield(obj,f,'color2');
                    else
                        col2 = col1;
                    end
                    if isfield(getfield(obj,f),'minr')
                        minr = getfield(obj,f,'minr');   %minimum radius
                        moment(target,magnitude,thickness,col,col2,minr(t_disp*SimSpeed))
                    else
                    moment(target,magnitude,thickness,col,col2)
                    end
                otherwise
                    error(['WARNING: An invalid object type chosen: "', type{ii}, '"'])
            end
        end
        end
        grid(config.grid)
        axis(config.axis)
    end

        

    
    if (isfield(config,'video')) && (config.video.enable == "on")
        %%%%%%%%%% Save as gif:
        start = 1;
        LengthOfAnimation = LengthOfAnimation + 1;
        if mod(LengthOfAnimation,80) == 1
            fprintf(repmat('\b',1,3))
        elseif mod(LengthOfAnimation,20) == 1
            fprintf('.')
        end
        movieVector(LengthOfAnimation) = getframe;
        t_disp = t_disp + config.video.resolution + ts/SimSpeed;
        %%%%%%%%%%
    elseif ~start && config.enterToStart
        figure(Fig)
        dummy = input('Press Enter to start the animation');
        fprintf('Animating... ')
        start = 1;
        tic
    else
        %Don't make a gif, but display on screen
        figure(Fig)
        t_disp = toc + ts/SimSpeed;
    end



    
 end
 

if config.video.enable == "on"
    movieVector = movieVector(1:LengthOfAnimation);
    myWriter = VideoWriter('AnimationVideo',config.video.profile);
    myWriter.FrameRate = 1/config.video.resolution;
    disp('done.')
    F = fieldnames(config.video);
    for j = 1:length(F)
        f = F{j};
        try
            switch f    
                case 'Quality'
                    myWriter.Quality = getfield(config.video,f);
                case 'CompressionRatio'
                    myWriter.CompressionRatio = getfield(config.video,f);
                case 'FrameRate'
                    myWriter.FrameRate = getfield(config.video,f);
                case 'LosslessCompression'
                        myWriter.LosslessCompression = getfield(config.video,f);
            end
        catch
            disp(['WARNING: Could not set the property: "', f, '", this property may be incompatible with the video format you have selected.'])
        end
    end
                
    fprintf('writing... ')
    open(myWriter)
    writeVideo(myWriter, movieVector)
    close(myWriter)
    disp('done.')
    disp(' ')
    disp('Your Video:')
    disp(myWriter)
    disp('NB! If you are satisfied with the video, change the filename, such that it will not be overridden by the next video you make.')
    ProcessingTime = toc(ProcessingTime);
    disp(['The processing time for this video: ', num2str(ProcessingTime), 'seconds'])
    disp('-----------')
else
    disp('done.')
end






end


function[] = ball(c,r,color)
%figure(Fig)
a = 0:0.01:pi*2;
L = length(a);
x = NaN(1,L);
y = x;
    for i = 1:L
        x(i) = c(1) + r*cos(a(i));
        y(i) = c(2) + r*sin(a(i));
    end
    hold on
    plot(x,y,'Color',color);
end

function[] = line(a,b,color)
%figure(Fig)

hold on
plot([a(1) b(1)],[a(2) b(2)],'-','Color',color,'Marker','.')
end

function[] = point(p,color)
%figure(Fig)
hold on
plot(p(1),p(2),'.','Color',color)
end

function[] = box(B1,B2,def,color,linestyle)
%figure(Fig)
if def == "CC" %(B1 - Center of box, B2 - any corner)
    B1 = 2*B1 - B2; %(opposite corner to B2)
elseif def == "CD" %(B1 - Center, B2 - Diagonal)
    B1 = B1 + 0.5*B2*arm(pi/6); %(place B1 at Upper Right corner)
    B2 = B1 - B2*arm(pi/6);     %(place B2 at Bottom Left corner)
elseif def == "2C"
    % DEFAULT
    % (B1 - some corner)
    % (B2 - opposite corner)
else
    error('WARNING: object of type "box" has invalid value in field "def".')
end

hold on
plot([B1(1) B2(1)],[B1(2) B1(2)],'Linestyle',linestyle,'Color',color)
plot([B2(1) B2(1)],[B1(2) B2(2)],'Linestyle',linestyle,'Color',color)
plot([B2(1) B1(1)],[B2(2) B2(2)],'Linestyle',linestyle,'Color',color)
plot([B1(1) B1(1)],[B2(2) B1(2)],'Linestyle',linestyle,'Color',color)
end

function[] = arrow(head,tail,color,R1,R2,arrow_arm)

hold on
plot([head(1) tail(1)],[head(2) tail(2)],'-','Color',color,'Marker','none')

a = head - tail;
Ra = R1*a*arrow_arm;
La = R2*a*arrow_arm;
Ra = Ra + head;
La = La + head;

plot([head(1) Ra(1)],[head(2) Ra(2)],'-','Color',color,'Marker','none')
plot([head(1) La(1)],[head(2) La(2)],'-','Color',color,'Marker','none')
end

function[] = curve(x,y,color,linestyle)
hold on
plot(x,y,'Color',color,'LineStyle',linestyle)
end


function[] = moment(target,magnitude,thickness,color,color2,minr)
if nargin < 6
    minr = 0;
end

a = -(0:0.01:(pi*2-1));
NEG = (magnitude < 0);
if NEG
    color = color2;
    a = -a;
end

L = length(a);
x = NaN(1,L);
y = x;
    for i = 1:L
        x(i) = target(1) + (minr + abs(magnitude))*cos(a(i));
        y(i) = target(2) + (minr + abs(magnitude))*sin(a(i));
    end
    hold on
    plot(x,y,'Color',color,'LineWidth',abs(thickness));

s = 0.3;
R = [x(1);y(1)] - s*(1-2*NEG)*(minr + abs(magnitude))*arm(   pi/4);
L = [x(1);y(1)] - s*(1-2*NEG)*(minr + abs(magnitude))*arm(pi-pi/4);

plot([x(1) R(1)],[y(1) R(2)],'-','Color',color,'Marker','none','LineWidth',abs(thickness))
plot([x(1) L(1)],[y(1) L(2)],'-','Color',color,'Marker','none','LineWidth',abs(thickness))

end



function[help] = Help()

help.Arguments = 'The required arguments are: Animate(tsim,qsim,obj,config)';
help.tsim = 'A column vector containing a series of timestamps (for their respective state-values)';
help.qsim = 'A matix containing the states of the animation at each timestamp in tsim. qsim = [q(t(1)); q(t(2)); q(t(3)); ...]';
help.obj = 'A struct containing all objects to be displayed in the animation.';
help.config = 'A struct containing several OPTIONAL fields for configuring the animation.';
help.object.fields = 'An object must have the field: type, AND any relevant *data-field* for that type. All objects can be colored using the field: *color*';
help.object.types = 'line, ball, point, arrow, and a timed_**** for all the types';
help.object.line = 'Requires: a,b - the positions of the ends of the line. BOTH are function handles that take the state vector as input';
help.object.ball = 'Requires: c,r - *c* is the center of the circle, and is a function handle that takes the state-vector as input.  *r* is the raduis of the circle, and is a constant';
help.object.point = 'Requires: p - the position of the point. *p* is a function handle that takes the state vector as input.';
help.object.arrow = 'Requires: head,tail - *head* is the position of the arrow head, where the arrow is pointing to. *tail* is the position of the back of the arrow. BOTH are function handles that take the state-vector as input';
help.object.curve = 'Requires: x,y - arrays of x-values and corresponding y-values. Can also take ''linestyle'' agrument ';
help.object.timed_ = 'the timed_**** versions of all types take in the time t as their second argument for all their fields(not color).';
help.object.color = 'use standard MATLAB color-codes';
help.configure.fields = 'tf, simspeed, axis, grid, enterToStart';
help.configure.tf = 'tf - *final time* - the time t at which to stop the animation. Make sure: ts < tf <= t(end).  (Defaults to tf = t(end))';
help.configure.ts = 'ts - *start time* - the time t at which to start the animation from. Make sure: ts < tf <= t(end).  (Defaults to ts = 0)';
help.configure.axis = 'axis - standard MATLAB axis command to set the window for display(recommended format: width = height*5/4).  (Defaults to *equal*)';
help.configure.enterToStart = 'enterToStart - if =1 - then the commandwindow asks you to press enter before running the animation.  (Defaults to =0)';
help.configure.simspeed = 'sets the relative speed to play the animation(relative to the timeseries t.). If =1,=0.5,=2, then the simulation plays in real time, half speed(slow motion), double speed(fast motion) respectively.  (Defaults to =1/real time)';

end

function[example] = Example()

% SIMULATE SYSTEM:
w = 2;
t = 0:0.01:10;
x1 = [cos(w.*t)./(0.1.*t+1);
     sin(w.*t)./(0.2.*t+1)];
x2 = [cos(2*w.*t);
      sin(2*w.*t)];
  
% CONFIGURE ANIMATION:
config.scale = 1;  %scale all coordinates
formatRatio = 5/4*1.55; %Recommended ratio for split screen figure
formatRatio = 5/4*0.75; %REcommended ratio for full screen figure
formatRatio = 5/4;      %Recommended ratio for standard pop-up figure
lift = 0;          %use to move in y-axis
shift = 0;         %use to move in x-axis
height = 3;        %use to set the window height
width = height*formatRatio;
config.axis = [-width/2 width/2 -height/2 height/2] + [shift shift lift lift];
config.simspeed = 1; %Set the simulation speed
config.ts = 0;       %Set start time of simulation
config.tf = t(end);  %Set finish time
config.grid = 'on';  %Enable grid
config.enterToStart = 1;  %Wait until you press enter before the simulation starts

% X1 :
obj.X1.type = {'ball','line'};
obj.X1.c = @(q) [q(1); q(2)];
obj.X1.r = 0.1;
obj.X1.a = @(q) obj.X1.c(q);
obj.X1.b = @(q)[0; 0];
obj.X1.color = 'b';

% Arrow:
obj.tau.type = 'arrow';
obj.tau.tail = @(q) obj.X1.c(q);
obj.tau.head = @(q) obj.X1.c(q) + 0.5*[-q(2); q(1)];
obj.tau.color = 'r';

% X2 :
obj.X2.type = 'timed_ball';
obj.X2.c = @(q,t) [q(3); q(4)];
obj.X2.r = @(q,t) 0.1*sin(3*t)+0.2;
darkmagenta = [139,0,139]/255;
obj.X2.color = darkmagenta;

example = 'See Example in code';
disp(example)
Animate(t',[x1;x2]',obj,config);

end



%% Description of video profiles:

% 'Archival'
% 
% Motion JPEG 2000 file with lossless compression
% 
% 'Motion JPEG AVI'
% 
% AVI file using Motion JPEG encoding
% 
% 'Motion JPEG 2000'
% 
% Motion JPEG 2000 file
% 
% 'MPEG-4'
% 
% MPEG-4 file with H.264 encoding (systems with Windows 7 or later, or macOS 10.7 and later)
% 
% 'Uncompressed AVI'
% 
% Uncompressed AVI file with RGB24 video
% 
% 'Indexed AVI'
% 
% Uncompressed AVI file with indexed video
% 
% 'Grayscale AVI'
% 
% Uncompressed AVI file with grayscale video
% 
% 








