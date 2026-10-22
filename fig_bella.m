function fig_bella(y_lim)
global x
    xlabel('$x$ [m]','Interpreter','latex')
    ylabel('$y$ [m]','Interpreter','latex')
    % set(gcf,'Color','white')
    if ~exist('y_lim','var')
        ylim([-1.2 1.2]*inf)
    else
        ylim(y_lim)
    end
    xlim(x([1 end]))
    grid on
    box on
end