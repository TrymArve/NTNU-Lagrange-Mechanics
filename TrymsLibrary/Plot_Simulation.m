function[Ref,States_figure] = Plot_Simulation(S,tsim,xsim,controller)
%%% PLOT SIMPULATED SYSTEMS

States_figure = figure('Position', [10+200 30 900 600]);

names = ["th","theta","lmb","lambda","w","omega","xi","eta","zeta","sigma","phi";
         "\theta","\theta","\lambda","\lambda","\omega","\omega","\xi","\eta","\zeta","\sigma","\phi"];

Tiles_States = tiledlayout(S.nq,4,'TileSpacing','Compact');
title(Tiles_States,'<- Coordinates/Velocities ->')
ylabel(Tiles_States,'Respective Units')
xlabel(Tiles_States,'Time [s]')

Names = S.GenCoords;
if nargin > 3
    Ref = controller.ref(tsim')';
end

for i = 1:S.nq
    
    color = GetColorCode(i);
    
    % correct names
    ind = find(names(1,:) == Names{i},1);
    if ~isempty(ind)
        Names{i} = names(2,ind);
    end
    
    % plot Coordinates
    tile = nexttile(i*4-3,[1 2]);
    hold on
    plot(tsim',xsim(:,i)','Color',color)
    Leg = char(['$$' Names{i} '$$']);
    Leg = reshape(Leg,1,[]);
        Title_temp = title(Leg);
        set(Title_temp,'Interpreter','latex');
    if nargin > 3
        plot(tsim',Ref(:,i)','Color',GetColorCode('pink'),'Linestyle','--')
        
        Leg2 = char(['$$' ,Names{i}, '_{ref}$$']);
        Leg2 = reshape(Leg2,1,[]);
        Leg_temp = legend(Leg,Leg2);
        set(Leg_temp,'Interpreter','latex');        
    end
    grid on


    % plot Velocities
    tile = nexttile(i*4-1,[1 2]);
    plot(tsim',xsim(:,S.nq+i)','Color',color)
    grid on
    Leg = char(['$$\dot{' Names{i} '}$$']);
    Leg = reshape(Leg,1,[]);
    Title_temp = title(Leg);
    set(Title_temp,'Interpreter','latex');

end
end
