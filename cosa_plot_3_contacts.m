function cosa_plot_3_contacts(W_n,t,DATA,y_lims)
    global y_lim x y velocities rho u v sigma sigma_c F11 F12 F21 F22 B1 B2
    global e_c dV glob rho0 a1 a2 a3
    clf
    set(gcf,'Position',get(0,'ScreenSize'))
    set(gcf,'Color','w')
    tt = tiledlayout(3,3,'TileSpacing','tight');
    % title(tt,['$t = ' num2str(t*1e3,3) '$ ms'],'Interpreter','latex')
    
    if ~exist('y_lims','var')
        y_lim_rho = 6000 + [-inf 10];
        y_lim_u = [-10 10]/100;
        y_lim_Yp = [0 1];
        y_lim_p = y_lim;
        y_lim_lmb = 1+[-1 1]/1000;
    else
        y_lim_rho = y_lims.rho;
        y_lim_u = y_lims.u;
        y_lim_Yp = y_lims.Yp;
        y_lim_p = y_lims.p;
        y_lim_lmb = y_lims.lmb;
    end
    
    [xx,yy] = meshgrid([x-dV.x/2, x(end)+dV.x(end)/2], [y-dV.y/2, y(end)+dV.y(end)/2]); % Creo la griglia che usa surf, cioè tutti i vertici in basso a sx delle celle
    n_g = 2; % Number of ghost cells
    x_box = [ xx([1+n_g end-n_g],1+n_g)' xx([1+n_g end-n_g],end-n_g)']; x_box(end+1) = x_box(1);
    y_box = [ yy(1+n_g,[1+n_g end-n_g])  yy(end-n_g,[1+n_g end-n_g]) ]; y_box = y_box([end 1:end]);
%%% Density
    nexttile
    zz = squeeze(rho(W_n))'; zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
    surf(xx, yy, zz, 'EdgeColor', 'none'); view(2)
    hold on
    plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
    colorbar;
    axis equal tight
    title(['Density [kg/m$^3$] \quad $t = ' num2str(t*1e3,3) '$ ms'],'Interpreter','latex');
    % clim(y_lim_rho)
    % Triangle and baricenter
    % bari_x = sum(zz(:).*xx(:)) / sum(zz(:));
    % bari_y = sum(zz(:).*yy(:)) / sum(zz(:));
    % plot3(bari_x,bari_y,2*max(zz,[],'all'),'ro','MarkerSize',2,'MarkerFaceColor','auto')
    % if exist('angolo','var')
    %     Lx = 0.2; Ly = 0.2/1.5;
    %     LL = 0.2;
    %     dd = LL/7;  ee = dd/2 *0;
    %     if angolo<pi/2
    %         % x_tri = [ dd, Lx/2, Lx-dd ];
    %         % y_tri = [ 0   tan(angolo)*(Lx/2-dd)    0];
    %         x_tri = [ dd, LL/2 ];
    %         y_tri = [ 0   tan(angolo)*(LL/2-dd) ]+ee;
    %     else
    %         % x_tri = [dd dd Lx-dd Lx-dd];
    %         % y_tri = [0 Ly/2 Ly/2 0];
    %         x_tri = [dd  dd  LL/2   ];
    %         y_tri = [0 Ly*2/3 Ly*2/3]+ee;
    %     end
    %     plot3(x_tri, y_tri, 2*max(zz,[],'all')*(x_tri*0+1), 'r')
    % end
    
%%% Speeds
    nexttile
    dovex = 3:round((length(x)-4)/20):length(x)-2;   dovey = 3:round((length(y)-4)/20):length(y)-2;
    uu = squeeze(u(W_n(:,dovex,dovey)))';  vv = squeeze(v(W_n(:,dovex,dovey)))';
    quiver(x(dovex), y(dovey), uu, vv, 'AutoScale', 'on', 'Color', 'r');
    hold on
    plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
    title(['Velocity Field [m/s] - max $' num2str(max(sqrt(uu.^2+vv.^2),[],'all'),3) '$ m/s'],'Interpreter','latex');
    fig_bella
    ylabel('$y$ [m]','Interpreter','latex');
    axis equal tight
    xlim([min(xx(:)),max(xx(:))]);   ylim([min(yy(:)),max(yy(:))])
    grid on


% %%% Density comparison
%     nexttile
%     rrrho = squeeze(rho(W_n))';
%     zG2 = squeeze(IG2(W_n))';
%     zJ = 1./sqrt(zG2);
% 
%     zz = rrrho.*zJ/rho0; zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
%     surf(xx, yy, zz, 'EdgeColor', 'none'); view(2)
%     hold on
%     plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
%     colorbar; axis equal tight
%     title('Density comparison [-]','Interpreter','latex');


% %%% Curl
%     nexttile
%     zz = squeeze(curl_b(W_n));  zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
%                                  zz(end-2,:) = zz(end-4,:); zz(:,end-2) = zz(:,end-4);
%                                  zz(end-3,:) = zz(end-4,:); zz(:,end-3) = zz(:,end-4);
%     surf(xx(3:end-3,3:end-3)+dV.x(1)/2, yy(3:end-3,3:end-3)+dV.y(1)/2, abs(zz(3:end-3,3:end-3)'), 'EdgeColor', 'none'); view(2)
%     hold on
%     max_z = 2*abs(max(zz(3:end-3,3:end-3),[],'all'));
%     plot3(x_box,y_box,max_z*ones(1,5),'k')
%     cb = colorbar;
%     fig_bella;
%     title(['$|\nabla\times a|$ -- $\oint\nabla\times a$ ='  num2str(glob.circu_a(end),2) ' -- $\oint\nabla\times b$ ='  num2str(glob.circu_b(end),2)],'Interpreter','latex');
%     plot3(x(round([end/3 end/2 end/2 end/3 end/3])),y(round([end*2/3 end*2/3 end/3 end/3 end*2/3 ])),max_z*ones(1,5),'k')
%     axis equal tight
%     xlim([xx(1,1) xx(end,end)]);  ylim([yy(1,1) yy(end,end)])
%     set(gca, 'ColorScale', 'log');


% %%% Deformation diagram
%     rrrho = squeeze(rho(W_n))';
%     % energia = squeeze(e_c(W_n))';
%             % IF1_tilde = @(W) sqrt(IG1_tilde(W)+2);
%             % test = squeeze(IG2(W_n).^(1/4) .* sqrt( ff(IF1_tilde(W_n)) ));
%     zG1_tilde = squeeze(IG1_tilde(W_n))'; % Queste sono invarianti di G
%     zG2 = squeeze(IG2(W_n))';
%     % zG2 = squeeze((a1(W_n).*b2(W_n)-a2(W_n).*b1(W_n)).^2)';
%     zF1_tilde = sqrt(zG1_tilde+2);
%     zJ = 1./sqrt(zG2);
% 
%     % xmax = max(zF1_tilde,[],'all');
%     % xmin = min(zF1_tilde,[],'all');
%     % ymax = max(zJ,[],'all');
%     % ymin = min(zJ,[],'all');
%     IIF1_tilde = [2+linspace(1e-3,4e-1,200), linspace(2+4e-1,23/3+3, 200)];
%     ff_val = ff(IIF1_tilde');
% 
%     nexttile
%     % Trasformo le invarianti di G = B^-1 date come input in invarianti di B
%     plot(IIF1_tilde-2,ff_val-1)
%     hold on
%     [rrrho,ord] = sort(rrrho(:));
%     scatter(zF1_tilde(ord)-2,  zJ(ord)-1   ,1,rrrho(:),'.')
%     % scatter(zF1_tilde(:)-2,  zJ(:)-1   ,1,energia(:),'.')
%     colorbar
%     xlabel('$Tr(F)\,/\,det(F)^{0.5}-2$','interpreter','latex');
%     ylabel('$det(F)-1$','interpreter','latex')
%     % axis square tight
%     % xlim([-1.2*(2-xmin), 1.2*(xmax-2)]+10*eps*[-1 1])
%     % ylim([-1.2*(1-ymin), 1.2*(ymax-1)]+10*eps*[-1 1])
%     set(gca,'XScale','log')
%     xlim([1e-3 23/3-2+3])
%     ylim([-0.05 0.25])
%     grid on


%%% sigma xx
    nexttile
    zz = squeeze(sigma.xx(W_n))'; zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
    % z = squeeze(p(W_n))';
    % z = squeeze(sigma_c.xx(W_n))';
    surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
    hold on
    plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
    fig_bella; title('$\sigma_{xx}$ [Pa]','Interpreter','latex')
    axis equal tight; 
    cb = colorbar;
    % clim(y_lim_p) % Occhio che sigma è negativo
    % set(cc, 'YTickLabel', cellstr(num2str(reshape(str2double(cc.TickLabels),[],1),'%0.1e')) )
    cb.Ruler.TickLabelFormat = '%.1f';
    % cb.Ruler.TickLabelFormat = '%.0e';
    % % Extract exponent from one of the tick labels
    % ticks = cb.Ticks;
    % expStr = regexp(sprintf('%.1e', max(abs(ticks))), 'e([+-]\d+)', 'tokens', 'once');
    % exponent = str2double(expStr{1});
    % % Remove the exponent part from tick labels
    % cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x / 10^exponent), ticks, 'UniformOutput', false);
    % % Add text '×10^n' on top of the colorbar
    % text(1.1, 0.9, ...
    %      sprintf('\\times10^{%d}', exponent), ...
    %      'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %      'FontSize', 10,'Interpreter', 'tex');

%%% sigma yy
    nexttile
    zz = squeeze(sigma.yy(W_n))'; zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
    surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
    hold on
    plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
    fig_bella; title('$\sigma_{yy}$ [Pa]','Interpreter','latex')
    axis equal tight; 
    cb = colorbar;
    cb.Ruler.TickLabelFormat = '%.1f';
    % cb.Ruler.TickLabelFormat = '%.0e';
    % Extract exponent from one of the tick labels
    % ticks = cb.Ticks;
    % expStr = regexp(sprintf('%.1e', max(abs(ticks))), 'e([+-]\d+)', 'tokens', 'once');
    % exponent = str2double(expStr{1});
    % % Remove the exponent part from tick labels
    % cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x / 10^exponent), ticks, 'UniformOutput', false);
    % % Add text '×10^n' on top of the colorbar
    % text(1.1, 0.9, ...
    %      sprintf('\\times10^{%d}', exponent), ...
    %      'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %      'FontSize', 10,'Interpreter', 'tex');

% %%% J
%     nexttile
%     zz = sqrt(squeeze((F11(W_n).*F22(W_n)-F12(W_n).*F21(W_n)).^2))';
%     % zz = zJ;
%     zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
%     surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
%     hold on
%     plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
%     colorbar
%     % hold on; plot(x,chi(1,:)/rho0,'o','MarkerSize',5)
%     fig_bella; title('$J$ [-]','Interpreter','latex')
%     axis equal tight

% %%% sigma xy
%     nexttile
%     % zz = squeeze(I1(W_n))';  zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
%     zz = squeeze(sigma.xy(W_n))'; zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
%     surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
%     hold on
%     plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
%     colorbar
%     % hold on; plot(x,chi(1,:)/rho0,'o','MarkerSize',5)
%     fig_bella; title('$\sigma_{xy}$ [Pa]','Interpreter','latex')
%     axis equal tight
%     clim([-1 1]*(max(abs(zz),[],'all')+2*eps))

%%% Parte immaginaria della velocità
    nexttile
    zz = squeeze(velocities.imag(W_n))';
    zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
    surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
    hold on
    plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
    colorbar
    % hold on; plot(x,chi(1,:)/rho0,'o','MarkerSize',5)
    fig_bella; title('max imag [m/s]','Interpreter','latex')
    axis equal tight
    clim([-1 1]*(max(abs(zz(3:end-3,3:end-3)),[],'all')+2*eps))

% %%% beta
%     nexttile
%     zz = squeeze(B1(W_n))'; zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
%     surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
%     hold on
%     plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
%     colorbar
%     % hold on; plot(x,chi(1,:)/rho0,'o','MarkerSize',5)
%     fig_bella; title('$\beta_{1}$ [-]','Interpreter','latex')
%     axis equal tight
%     clim([-1 1]*(max(abs(zz),[],'all')+2*eps))

%%% Kinetic energy
    nexttile
    uu = squeeze(u(W_n))';  vv = squeeze(v(W_n))';
    uu([1:2 end-1:end],:) = 0;  uu(:,[1:2 end-1:end]) = 0;
    vv([1:2 end-1:end],:) = 0;  vv(:,[1:2 end-1:end]) = 0;
    zz = 1/2 * squeeze(rho(W_n))' .* ( uu.^2 + vv.^2 ); zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
    surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
    hold on
    plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
    colorbar
    % hold on; plot(x,chi(1,:)/rho0,'o','MarkerSize',5)
    fig_bella; title('$\mathcal K$ [J/m$^3$]','Interpreter','latex')
    axis equal tight

%%% Characteristic speeds x
    nexttile
    zz = squeeze(velocities.max_x(W_n))';
    zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
    surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
    hold on
    plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
    colorbar
    fig_bella; title('Velocities x [m/s]','Interpreter','latex')
    axis equal tight
%%% Characteristic speeds y
    nexttile
    zz = squeeze(velocities.max_y(W_n))';
    zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
    surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
    hold on
    plot3(x_box,y_box,2*abs(max(zz,[],'all'))*ones(1,5),'k')
    colorbar
    fig_bella; title('Velocities y [m/s]','Interpreter','latex')
    axis equal tight
    
%%% Numero di contatti
    axx = nexttile;
    [contact1,b1_1,b1_2] = trova_conta(a1);
    [contact2,b2_1,b2_2] = trova_conta(a2);
    [contact3,b3_1,b3_2] = trova_conta(a3);
    zz = contact1 + contact2 + contact3;
    zz(end+1,:) = zz(end,:); zz(:,end+1) = zz(:,end);
    surf(xx,yy,zz, 'EdgeColor', 'none'); view(2)
    hold on
    max_z = 2*abs(max(zz,[],'all'));
    plot3(x_box,y_box,max_z*ones(1,5),'k')
    clim([-0.5 3.5])
    % Discrete colormap (4 colors)
    cmap = parula(4);
    colormap(axx,cmap)
    
    % Create colorbar
    cb = colorbar;
    
    % Ticks at class centers
    cb.Ticks = 0:3;
    
    % % Optional labels
    % cb.TickLabels = {'A','B','C','D','E'};
    % hold on; plot(x,chi(1,:)/rho0,'o','MarkerSize',5)
    fig_bella; title('$N$ contacts [-]','Interpreter','latex')
    axis equal tight


    ogni = floor(max([size(zz,1),size(zz,2)])/12);
    s = ogni/2*dV.x(1) / DATA.R;   % size of the local symbol
    for iy = floor(ogni/2) : ogni : size(zz,2)-floor(ogni/2)
        for ix = floor(ogni/2) : ogni : size(zz,1)-floor(ogni/2)
            xp = xx(ix,iy);
            yp = yy(ix,iy);
    
            % horizontal "--"
            if contact1(ix,iy)
                % plot3(xp+s*a1(1)*[-1 1]/2,yp+s*a1(2)*[-1 1]/2,max_z*[1 1],'k','LineWidth',1)
                plot3(xp +s*b1_1(ix,iy)*[-1 1]/2,  yp +s*b1_2(ix,iy)*[-1 1]/2,max_z*[1 1],'r','LineWidth',.8)
            end
            % "/"
            if contact2(ix,iy)
                % plot3(xp+s*a2(1)*[-1 1]/2,yp+s*a2(2)*[-1 1]/2,max_z*[1 1],'k','LineWidth',1)
                plot3(xp +s*b2_1(ix,iy)*[-1 1]/2,  yp +s*b2_2(ix,iy)*[-1 1]/2,max_z*[1 1],'r','LineWidth',.8)
            end
            % "\"
            if contact3(ix,iy)
                % plot3(xp+s*a3(1)*[-1 1]/2,yp+s*a3(2)*[-1 1]/2,max_z*[1 1],'k','LineWidth',1)
                plot3(xp +s*b3_1(ix,iy)*[-1 1]/2,  yp +s*b3_2(ix,iy)*[-1 1]/2,max_z*[1 1],'r','LineWidth',.8)
            end    
        end
    end
    

    function [out,b1,b2] = trova_conta(a)
        corrector = rho0./rho(W_n) ./( F11(W_n).*F22(W_n)-F12(W_n).*F21(W_n) );
        b1 = ( a(1).*F11(W_n)+a(2).*F12(W_n)  ) .*corrector; % I estimate the determinant from the density
        b2 = ( a(1).*F21(W_n)+a(2).*F22(W_n)  ) .*corrector;
        bb = sqrt( b1.^2 +  b2.^2);
        xi = sqrt(a(1)^2 + a(2)^2) - bb;
        out = squeeze(xi >=0)';
        b1 = squeeze(b1)';
        b2 = squeeze(b2)';
    end
end
