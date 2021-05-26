function[] = Animate(tsim,qsim,obj,config)

Fig = figure;

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

axis(config.axis)
SimSpeed = config.simspeed;
tf = config.tf;


% Animation:
t_disp = 0;
tic
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
                    ball(c,r,col,Fig)

                case 'line'
                    
                    h = getfield(obj,f,'a');
                    a = h(q);
                    h = getfield(obj,f,'b');
                    b = h(q);
                    line(a,b,col,Fig)

                case 'point'
                    h = getfield(obj,f,'p');
                    p = h(q);
                    point(p,col,Fig)
                    
                otherwise
                    disp('WARNING: No valid object type chosen')
            end
        end
        grid(config.grid)
        axis(config.axis)
    end
    figure(Fig)
 
    t_disp = toc;
 end
 
 
end


function[] = ball(c,r,color,Fig)
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

function[] = line(a,b,color,Fig)
%figure(Fig)

hold on
plot([a(1) b(1)],[a(2) b(2)],'-','Color',color,'Marker','.')
end

function[] = point(p,color,Fig)
%figure(Fig)
hold on
plot(p(1),p(2),'.','Color',color)
end
