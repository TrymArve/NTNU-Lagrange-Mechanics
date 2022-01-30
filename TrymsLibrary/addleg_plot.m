function[] = addleg_plot(x,y,style,mark,addleg)
% o	Circle
% +	Plus sign
% *	Asterisk
% .	Point
% x	Cross
% s	Square
% d	Diamond
% ^	Upward-pointing triangle
% v	Downward-pointing triangle
% >	Right-pointing triangle
% <	Left-pointing triangle
% p	Pentagram
% h	Hexagram

% "addleg" - can be given as a cell input on the form:
% {'Ball','Coconut','Perfume';
%     2    ,   4    ,   3     }
% This example requires 9 inputs to plot across x. (y has 9 sets of values)
% Legends will be: 

% Ball_1
% Ball_2
% Coconut_1
% Coconut_2
% Coconut_3
% coconut_4
% Perfume_1
% Perfume_2
% Perfume_3


%colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560], ...
          [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0.5 0 0.5], ...
          [1 1 1]*0.5,[0.9686 0.498 0.7451]};
code = ['b','o','y','p','g','c','r','m','e','i'];
% e - grey
% i - pink
% p - purple
ind = find(code == style(1));
if ~isempty(ind)
    color = colors{ind};
else
    color = style(1);
end

p = size(y,1);

hold on

if p == 1
L = legend();
s = L.String;
if mark == ';'
plot(x,y,'Color',color,'LineStyle',style(2:end))
else
plot(x,y,'Color',color,'LineStyle',style(2:end),'Marker',mark)
end
legend([s, {addleg}])

elseif p > 1
    
    if iscell(addleg)
        [P,xx] = size(addleg);
    else
        P = 1;
        addleg = {addleg; p};
    end 
    plotcounter = 0;
    for ii = 1:P
        addleg_extract = addleg{1,ii};

        for i = 1:addleg{2,ii}
            plotcounter = plotcounter + 1;
        L = legend();
        s = L.String;
        if mark == ';'
            plot(x,y(plotcounter,:),'Color',color,'LineStyle',style(2:end))
        else
            plot(x,y(plotcounter,:),'Color',color,'LineStyle',style(2:end),'Marker',mark)
        end
        color = mod(color + [0.23 0 -0.13], 1);
        number = num2str(i);
        ADD = [addleg_extract, '_', number];
        legend([s, {ADD}])
        end

    end

else
msg = ['WARNING: Could not plot' , addleg, '. Check you y array'];
disp(msg)
end

grid on
end
