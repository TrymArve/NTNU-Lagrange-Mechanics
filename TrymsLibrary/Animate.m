function[Fig] = Animate(tsim,qsim,obj,config)

Fig = figure;

% rotation matrices:
theta = pi/4+pi-pi/12;
R1 = [cos(theta) -sin(theta);
     sin(theta) cos(theta)];
theta = -theta;
R2 = [cos(theta) -sin(theta);
     sin(theta) cos(theta)];
arrow_arm = 0.2;
 
 
%Confugitation:
if ~isfield(config,'axis')
config.axis = 'equal';
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

axis(config.axis)
SimSpeed = config.simspeed;
tf = config.tf;


% Animation:
t_disp = 0;
tic
start = 0;
 while t_disp < tf/SimSpeed
    clf
    
    q   = interp1(tsim,qsim,SimSpeed*t_disp)';
    
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
        for ii = 1:L
            if isfield(getfield(obj,f),'color')
                col = getfield(obj,f,'color');
            else
                col = 'k';
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
                    
                case 'timed_arrow'
                    
                    h = getfield(obj,f,'head');
                    head = h(q,t_disp*SimSpeed);
                    h = getfield(obj,f,'tail');
                    tail = h(q,t_disp*SimSpeed);
                    arrow(head,tail,col,R1,R2,arrow_arm);
                    
                otherwise
                    error('WARNING: No valid object type chosen')
            end
        end
        grid(config.grid)
        axis(config.axis)
    end
    figure(Fig)
 
    t_disp = toc;
    
    
    if ~start && config.enterToStart
        dummy = input('Press Enter to start the animation');
        disp('Starting animation--')
        start = 1;
        t_disp = 0;
        tic
    end
    
 end
 
 disp('--Animation finished')
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
