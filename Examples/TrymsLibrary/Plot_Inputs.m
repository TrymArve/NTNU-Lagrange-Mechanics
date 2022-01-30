function[inputs,Inputs_figure] = Plot_Inputs(tsim,usim)
%%% PLOT CONTROL INPUTS FROM SIMULATION

Inputs_figure = figure('Position', [10 30 900 600]);


%Clean Input Series;
inputs = CleanSeries(usim,tsim); %Removes duplicate timestamps and interpolates to fit tsim.

nu = length(inputs(1,:));

if nu > 4
    Tiles_Inputs = tiledlayout(nu,2,'TileSpacing','Compact');
else
    Tiles_Inputs = tiledlayout(nu,1,'TileSpacing','Compact');
end

title(Tiles_Inputs,'Control Inputs')
ylabel(Tiles_Inputs,'Respective Units')
xlabel(Tiles_Inputs,'Time [s]')

for i = 1:nu
    
    % plot
    nexttile
    plot(tsim',inputs(:,i)','Color',GetColorCode(i))
    grid on
    title(['u_', num2str(i)])
end
    
end








