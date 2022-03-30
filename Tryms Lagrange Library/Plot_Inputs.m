function[inputs,Inputs_figure] = Plot_Inputs(tsim,usim,mpc)
%%% PLOT CONTROL INPUTS FROM SIMULATION

Inputs_figure = figure('Position', [10 30 900 600]);

if (nargin == 3) && (mpc == 1)
    inputs = usim;
else
%Clean Input Series;
inputs = CleanSeries(usim,tsim); %Removes duplicate timestamps and interpolates to fit tsim.
end

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
    if (nargin == 3) && (mpc == 1)
        % Add padding:
        M = max(inputs(:,i)');
        N = min(inputs(:,i)');
        A = (M+N)/2;
        W = (M-A);
        plot([0 0],[N-W*0.05 M+W*0.05],'.','Color','white')
        hold on
        % Plot control inputs:
        stairs(tsim',inputs(:,i)','-','Linewidth',2,'Color',GetColorCode(i,1.2))
    else
        plot(tsim',inputs(:,i)','Color',GetColorCode(i))
    end
    grid on
    title(['u_', num2str(i)])
    
end
    





end








