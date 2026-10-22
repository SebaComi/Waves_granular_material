% We consider the system obtained when the energy comes from three contacts
% rho_t + div(rho u) = 0
% (rho u)_t + div(rho u^2 - sigma ) = - g rho - gamma rho u
% (rho F)_t + div(rho F \otimes u - rho u F^-T) = 0
% We do splitting with respect to [rho, rho u] and [rho F]. In this way, we
% apply Bale's scheme to the first part (that has a static RHS and has a simple eigenstructure) and Rusanov
% or equivalents to the second part (which has no RHS and a complex eigenstructure)

clc; clear; close all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaulttextinterpreter','latex')
% addpath('Schemes\','Figure\','Funs\','Comi 3 contacts\','Comi 3 contacts\Data\')
tic
global velocities x y
global rho0 F11 F12 F21 F22
global rho u v sigma e_c a1 a2 a3 
global glob gg dV
global DATA

importa = false;
% importa = true;
if importa
    % load('Melin 1994 – 256x256 grains - 84 - precopr=2.20e-04.mat');
    % load('Melin 1994 – 256x256 grains - 84 - precopr=2.40e-04.mat');
    % load('Melin 1994 – 256x256 grains - 44 - precopr=2.40e-04.mat');
    % load('Melin 1994 – 256x256 grains - 104 - precopr=2.40e-04.mat');
    load('Melin 1994 – 256x256 grains - 204 - precopr=2.15e-04_T=0.05.mat');
    % load('Melin 1994 – 256x256 grains - 204 - precopr=2.15e-04_T=0.1.mat');
    % Here the data that I change
    DATA.T = 0.0125; % [s]    simulation time
    DATA.fric = @(t) 0;  % [Ns/m] arbitrary dissipation

    DATA.U0 = 2*0.05 *2*DATA.R; % [m] Spostamento imposto
    DATA.sigg = 5e-4;      % [s] Gaussian width
    DATA.t0 = 3*DATA.sigg; % [s] Gaussian center
    DATA.vv = @(t) DATA.U0/DATA.sigg/sqrt(2*pi) * exp(-(t-DATA.t0).^2/2/DATA.sigg^2);
    plots = false;
    % plots = true;

    DATA.time_to_save = 0;
    DATA.time_saved = [];
    snap1.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_rho_t='      num2str(t,2) 's'];  snap1.plot = @ plot_rho_3_contacts;
    snap2.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_kinetic_t='  num2str(t,2) 's'];  snap2.plot = @(W,t) plot_kinetic_3_contacts(W,t);
    snap3.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_sigma.xx_t=' num2str(t,2) 's'];  snap3.plot = @(W,t) plot_sigma_3_contacts(W,t,'xx');
    snap4.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_sigma.xy_t=' num2str(t,2) 's'];  snap4.plot = @(W,t) plot_sigma_3_contacts(W,t,'xy');
    snap5.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_sigma.yy_t=' num2str(t,2) 's'];  snap5.plot = @(W,t) plot_sigma_3_contacts(W,t,'yy');
    snap6.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_F11_t=' num2str(t,2) 's'];  snap6.plot = @(W,t) plot_F_3_contacts(W,t,'11');
    snap7.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_F12_t=' num2str(t,2) 's'];  snap7.plot = @(W,t) plot_F_3_contacts(W,t,'12');
    snap8.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_F21_t=' num2str(t,2) 's'];  snap8.plot = @(W,t) plot_F_3_contacts(W,t,'21');
    snap9.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_F22_t=' num2str(t,2) 's'];  snap9.plot = @(W,t) plot_F_3_contacts(W,t,'22');
    snap10.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_vx_t=' num2str(t,2) 's'];  snap10.plot = @(W,t) plot_v_3_contacts(W,t,'max_x');
    snap11.name = @(t) ['Comi 3 contacts/Snapshots/256x256-' num2str(DATA.Nx) '_vy_t=' num2str(t,2) 's'];  snap11.plot = @(W,t) plot_v_3_contacts(W,t,'max_y');
    DATA.snaps = {snap1, snap2, snap3, snap4, snap5, snap6, snap7, snap8, snap9, snap10, snap11};
    DATA.T = 0.0001; % [s]    simulation time

    DATA.time_to_save = (1.5:1.5:12)/1000;
    DATA.T = DATA.time_to_save(end)+0.1/1000; % [s]    simulation time
    snap2.name = @(t) ['Comi 3 contacts/Snapshots/256x256-104_kinetic_t='  num2str(t,3) 's'];  snap2.plot = @(W,t) plot_kinetic_3_contacts(W,t,'time');
    DATA.snaps = {snap2};
else
    DATA.gg = 9.81; % [m/s] gravity
    DATA.ggtheta = -pi/2*0-pi/2; % [rad] gravity angle with semipos x-axis
    DATA.R    = 0.005; %; % [m]   Particles radius by Melin
    DATA.rot = pi/4*0;
    DATA.a1 = [cos(DATA.rot), -sin(DATA.rot); sin(DATA.rot),cos(DATA.rot) ] * [1 0]' * DATA.R;
    DATA.a2 = [cos(pi/3), -sin(pi/3); sin(pi/3),cos(pi/3) ] * DATA.a1;
    DATA.a3 = [cos(pi/3), -sin(pi/3); sin(pi/3),cos(pi/3) ] * DATA.a2;
    DATA.rho_p = 7650; % [kg/m3] Density of the composing material
    DATA.rho0 = DATA.rho_p*sqrt(3)*pi/9; % [kg/m3] Reference density
    DATA.E_p  = 2.26e11;   % [Pa]  Young's modulus
    DATA.nu_p = 0.3;   % [-]   Poisson's ration
    DATA.n    = 3/2;   % [-]   Hertz's exponent
    DATA.kH   = sqrt(2)*DATA.E_p/(4*pi*DATA.rho_p*(1-DATA.nu_p^2)*DATA.R^2.5); % [Pa m^0.5 / kg] Hertz's constant

    % DATA.CFL = 0.96/2;
    DATA.CFL = 0.6;
    DATA.Nx = 50 + 4;

    DATA.precompr = 1 + 2.15e-4;

    DATA.vv = @(t) 0*t;

    DATA.Lx = 256*DATA.R*2     ; % [m]
    % DATA.Ly = DATA.Lx;
    DATA.Ly = DATA.Lx/2*sqrt(3);
    DATA.dx = DATA.Lx/(DATA.Nx-4);
    DATA.Ny = round(DATA.Ly/DATA.dx) + 4;
    DATA.dy = DATA.Ly/(DATA.Ny-4);

    DATA.T = 0.05; % [s]
    % DATA.T = 0.1; 
    DATA.fric = @(t) 800*(t<0.05); % [Ns/m] arbitrary dissipation

    % DATA.scheme_type = {'Bale','Ru2D'};
    DATA.scheme_type = {'Bale','MUSCL_H2D'};
    % DATA.scheme_type = 'MUSCL_H2D_curl';
    % DATA.scheme_type = 'MUSCL_H2D_curl_imp';
    plots = false;
    % plots = true;
end

salva = false;
% salva = true;

%% Importo i dati
T = DATA.T;
Nx = DATA.Nx;   Ny = DATA.Ny;
Lx = DATA.Lx;   Ly = DATA.Ly;
dx = DATA.dx;   dy = DATA.dy;
rho0 = DATA.rho0;
gg = DATA.gg;  theta = DATA.ggtheta;
a1 = DATA.a1;  a2 = DATA.a2;  a3 = DATA.a3;
CFL = DATA.CFL; % CFL
scheme_type = DATA.scheme_type;
fric = DATA.fric; % Arbitrary dissipation 

%% FEM Data
[x,y] = deal( (0:Nx-1)*dx+dx/2-2*dx,  (0:Ny-1)*dy+dy/2-2*dy );
% [x,y] = deal( (0:Nx-1)*dx+dx/2,  [0:floor(Ny/2)  (floor(Ny/2)+(Ny/2)^0.5*(0:floor(Ny/2)).^0.5)] *dy+dy/2 );
dV = [];
dV.x = [dx diff(x)]; % Questi sono i volumi delle celle
dV.y = [dy diff(y)];
[xx,yy] = meshgrid(x,y); xx = xx'; yy = yy';
dV.xy = zeros(2,Nx,Ny);   dV.xy(1,:,:) = xx;   dV.xy(2,:,:) = yy;
%% Variables and equations

% The vector W = [rho, rho u, rho v, rho F11, rho F12, rho F21, rho F22]
rho = @(W) W(1,:,:);
u   = @(W) W(2,:,:)./W(1,:,:);
v   = @(W) W(3,:,:)./W(1,:,:);
F11 = @(W) W(4,:,:)./W(1,:,:);
F12 = @(W) W(5,:,:)./W(1,:,:);
F21 = @(W) W(6,:,:)./W(1,:,:);
F22 = @(W) W(7,:,:)./W(1,:,:);

[sigma,e_c] = fun_cont_hertz_homo(rho,F11,F12,F21,F22,DATA);

restr = { 1:3 ,  4:7 };


%%% Funzioni per la prima parte dello splitting (Bale)
    %       f1 =    [rho u,         rho u^2 - sigma_xx,             rho u v - sigma_yx,             ]
    f1.x  = @(W,xy) [rho(W).*u(W);  rho(W).*u(W).^2-sigma.xx(W);    rho(W).*u(W).*v(W)-sigma.yx(W)  ];
    %       f1 =    [rho u,         rho u v - sigma_xy,             rho v^2 -sigma_yy,              ]
    f1.y  = @(W,xy) [rho(W).*v(W);  rho(W).*v(W).*u(W)-sigma.xy(W); rho(W).*v(W).^2-sigma.yy(W);    ];
    f1.u = u; f1.v = v; f1.rho = rho;
    
    % f1.A = @(W,xy) A_fun(W,a1,a2,a3,rho,u,v,F11,F12,F21,F22,DATA);
    % f1.B = @(W,xy) B_fun(W,a1,a2,a3,rho,u,v,F11,F12,F21,F22,DATA);
    f1.eigA = @(W,xy) eigA_fun_spl1(W,a1,a2,a3,rho,u,v,F11,F12,F21,F22,DATA,restr{1});
    f1.eigB = @(W,xy) eigB_fun_spl1(W,a1,a2,a3,rho,u,v,F11,F12,F21,F22,DATA,restr{1});
    % Divido il RHS in due componenti in modo che RHS = phix+phiy.
    % Quella che mi aspetto sia bilanciata dal flusso lungo x è phix, l'altra è phiy. 
    %   g  =            [ 0                       rho     (g cos(theta) - fric u )        0    0    0    0   0   ];
    f1.phix = @(W,xy,t) [zeros(size(W(1,:,:)));   rho(W).*(gg*cos(theta)-fric(t)*u(W));   zeros(size(W(1:5,:,:))) ];
    f1.dphix_dq = @(W,xy,t) dphix_dq_fun(W,xy,t,gg*cos(theta),fric);
    %   g  =            [ 0               0       rho     (g sin(theta) - fric v )        0     0     0     0   ];
    f1.phiy = @(W,xy,t) [zeros(size(W(1:2,:,:))); rho(W).*(gg*sin(theta)-fric(t)*v(W));   zeros(size(W(1:4,:,:))) ];
    f1.dphiy_dq = @(W,xy,t) dphiy_dq_fun(W,xy,t,gg*sin(theta),fric);

    g1 = [];
    if false % set true to run the scheme that is not well-balanced
        f1 = rmfield(f1,{'phix','phiy','dphix_dq','dphiy_dq'});
        g1 = @(t,W) [zeros(size(W(1,:,:)));  rho(W).*(gg*cos(theta)-fric(t)*u(W)); rho(W).*(gg*sin(theta)-fric(t)*v(W)) ];
    end

%%% Funzioni per la seconda parte dello splitting (Rusanov)
    %       f2 =    [rho(F11 u1 -F11 u1),  rho(F12 u1-F12 u1),  rho     (F21 u1    -  F11 u2),        rho     (F22 u1      -F12 u2),       ]
    f2.x  = @(W,xy) [zeros(size(W(1:2,:,:)));               rho(W).*(F21(W).*u(W)-F11(W).*v(W));  rho(W).*(F22(W).*u(W)-F12(W).*v(W))  ];
    %       f2 =    [rho     (F11 u2      -F21 u1),         rho     (F12 u2      -F22 u1),       rho(F21 u2 -F21 u2),  rho(F22 u2 -F22 u2) ]
    f2.y  = @(W,xy) [rho(W).*(F11(W).*v(W)-F21(W).*u(W));   rho(W).*(F12(W).*v(W)-F22(W).*u(W));      zeros(size(W(1:2,:,:)));             ];
    f2.u = u; f2.v = v; f2.rho = rho;
    
    %       g =  [   -u1 (rho Fj1)_,j      ,        -u1 (rho Fj2)_,j      ,       -u2 (rho Fj1)_,j      ,        -u2 (rho Fj2)_,j  ]
    g2 = @(t,W) [ -u(W).*num_div(rho(W).*F11(W),rho(W).*F21(W));   -u(W).*num_div(rho(W).*F12(W),rho(W).*F22(W));   -v(W).*num_div(rho(W).*F11(W),rho(W).*F21(W));   -v(W).*num_div(rho(W).*F12(W),rho(W).*F22(W)) ];

ftot = {f1, f2};
gtot = {g1, g2};

velocities = [];
velocities.x   = @(W) c_fun_Comi_3_contacts(W,rho,u,v,F11,F12,F21,F22,DATA,a1,a2,a3,'x') ;
velocities.y   = @(W) c_fun_Comi_3_contacts(W,rho,u,v,F11,F12,F21,F22,DATA,a1,a2,a3,'y') ;
velocities.max_x = @(W) c_fun_Comi_3_contacts(W,rho,u,v,F11,F12,F21,F22,DATA,a1,a2,a3,'max_x') ;
velocities.max_y = @(W) c_fun_Comi_3_contacts(W,rho,u,v,F11,F12,F21,F22,DATA,a1,a2,a3,'max_y') ;
velocities.max = @(W) c_fun_Comi_3_contacts(W,rho,u,v,F11,F12,F21,F22,DATA,a1,a2,a3,'max') ;
velocities.imag = @(W) c_fun_Comi_3_contacts(W,rho,u,v,F11,F12,F21,F22,DATA,a1,a2,a3,'imag') ;

veloc1 = velocities;  % this is not important because Bale's uses the velocities from the eigenvalues directly
veloc2 = velocities;  % This seems critical as the velocities are used by the Rusanov

% veloc1 = @(W) veloc_fun(W,a1,a2,a3,rho,u,v,F11,F12,F21,F22,DATA);

% veloc2.x = @(W) [u(W); -u(W)];  % If I give this to the Rusanov, it makes no sense and explodes
% veloc2.x = @(W) [v(W); -v(W)];
% veloc2.max_x = @(W) abs(u(W));
% veloc2.max_y = @(W) abs(v(W));

velocities_tot = {veloc1, veloc2};

%% Initial conditions
if importa
    W0 = W_n1;
else
    dove = true(size(yy));                                              % All grains
    rho_0 = zeros(Nx,Ny);  u_0 = rho_0; v_0 = u_0; F11_0 = u_0; F12_0 = u_0; F21_0 = u_0;  F22_0 = u_0;
    % Scelgo le densità
    rho_0(dove) = rho0* DATA.precompr;
    % rho1 = 200;
    % rho_0(~dove) = rho1;
       
    F11_0(:) = 1;   F12_0(:) = 0;
    F21_0(:) = 0;   F22_0(:) = 1;
    
    % Dove ho poca densità, i grani sono molto più sparsi, più distanti
    determ = F11_0.^2+F12_0.^2;
    F11_0 = F11_0./ sqrt(determ)  ./ (rho_0/rho0).^0.1;
    F12_0 = F12_0./ sqrt(determ)  ./ (rho_0/rho0).^0;
    F21_0 = F21_0./ sqrt(determ)  ./ (rho_0/rho0).^0;
    F22_0 = F22_0./ sqrt(determ)  ./ (rho_0/rho0).^0.9;
    
    % "Conservative form"
    W0 = zeros(7,Nx,Ny);
    W0(1,:,:) = rho_0;
    W0(2,:,:) = u_0.*rho_0;
    W0(3,:,:) = (0.01*0 + v_0).*rho_0;
    W0(4,:,:) = F11_0.*rho_0;
    W0(5,:,:) = F12_0.*rho_0;
    W0(6,:,:) = F21_0.*rho_0;
    W0(7,:,:) = F22_0.*rho_0;
end

%% Boundary conditions
BCs = [];
% BCs.L.type = 'periodic';
vv = DATA.vv;
BCs.L.type = 'wall';
BCs.L.v = @(W,t) [W(1,:);  2*vv(t)-W(2,:);   W(3,:);   W(4,:);-W(5,:);-W(6,:);W(7,:)];
% BCs.R.type = 'periodic';
BCs.R.type = 'wall';
BCs.R.v = @(W,t) [W(1,:);         -W(2,:);   W(3,:);   W(4,:);-W(5,:);-W(6,:);W(7,:)];
BCs.B.type = 'wall';
BCs.B.v = @(W,t) [W(1,:);          W(2,:);  -W(3,:);   W(4,:);-W(5,:);-W(6,:);W(7,:)];
BCs.T.type = 'wall';
BCs.T.v = @(W,t) [W(1,:);          W(2,:);  -W(3,:);   W(4,:);-W(5,:);-W(6,:);W(7,:)];

BCs.L.even = logical([ 1 1 0 1 0 0 1 ]);
BCs.R.even = logical([ 1 1 0 1 0 0 1 ]);
BCs.B.even = logical([ 1 0 1 1 0 0 1 ]);
BCs.T.even = logical([ 1 0 1 1 0 0 1 ]);


% These are the BCs for the RHS of the Bale's scheme
% Essentially, I invert the direction of gravity
BCs.BCsPHIx.L.type = 'wall';
BCs.BCsPHIx.L.v = @(PHI,t) [ PHI(1,:);  -PHI(2,:);   PHI(3:end,:) ];
BCs.BCsPHIx.R.type = 'wall';
BCs.BCsPHIx.R.v = @(PHI,t) [ PHI(1,:);  -PHI(2,:);   PHI(3:end,:) ];
BCs.BCsPHIx.B.type = 'wall';
BCs.BCsPHIx.B.v = @(PHI,t) [ PHI(1,:);   PHI(2,:);   PHI(3:end,:) ];
BCs.BCsPHIx.T.type = 'wall';
BCs.BCsPHIx.T.v = @(PHI,t) [ PHI(1,:);   PHI(2,:);   PHI(3:end,:) ];

BCs.BCsPHIy.L.type = 'wall';
BCs.BCsPHIy.L.v = @(PHI,t) [ PHI(1:2,:);   PHI(3,:);    PHI(4:end,:) ];
BCs.BCsPHIy.R.type = 'wall';
BCs.BCsPHIy.R.v = @(PHI,t) [ PHI(1:2,:);   PHI(3,:);    PHI(4:end,:) ];
BCs.BCsPHIy.B.type = 'wall';
BCs.BCsPHIy.B.v = @(PHI,t) [ PHI(1:2,:);  -PHI(3,:);    PHI(4:end,:) ];
BCs.BCsPHIy.T.type = 'wall';
BCs.BCsPHIy.T.v = @(PHI,t) [ PHI(1:2,:);  -PHI(3,:);    PHI(4:end,:) ];
%% Name and scheme
filename = ['Comi_3_cont CFL=' num2str(CFL) ' Nx=' num2str(Nx) ' Melin' ' sch=' strjoin(scheme_type, '_') '.gif'];
scheme = set_scheme(scheme_type,[],[],ftot,[],[],gtot,velocities_tot,BCs,restr);

%% Solver
y_lims.rho = 6000 + [-inf 15];
y_lims.u = [-10 10]/100;
y_lims.Yp = [-inf inf];
y_lims.p = [0 inf];
y_lims.lmb = 1+[-1 1]*1e-6*inf;
y_lim = inf*[-1 1];

[W_n1,dV] = solver(W0,dV,T,CFL,scheme,@(W_n,t) cosa_plot_3_contacts(W_n,t,DATA,y_lims),salva,plots,filename,'Eul',[], @ global_fun);

% return
%% Figura dell'energia
figure
plot(glob.t,glob.ep_c,   glob.t,glob.ep_gg,   glob.t, glob.ek,   glob.t,glob.ek+glob.ep_c+glob.ep_gg, 'LineWidth',1.2)
hold on
% plot(glob.t,glob.ep_c+glob.ek,'LineWidth',1.2)
ylabel('Energy [J]')
% ylim([1e-5 1e3])
xlim(glob.t([1 end]))
set(gca,'YScale','log')
grid on

yyaxis right
plot(glob.t, glob.m,'LineWidth',1.21)
ylabel('Mass [kg / m]')
xlabel('Time [s]');  
legend('E_c','E_{gg}','E_K','E_{tot}','m')
if salva
    saveas(gcf,[filename(1:end-21) '_energy' filename(end-20:end-3) 'png'])
end
toc

return
%% Esporto i dati a fine simulazione
% save(['Comi 3 contacts\Data\Melin 1994 – 256x256 grains - ' num2str(Nx) ' - precopr=' sprintf('%.2e',DATA.precompr-1) '_T=' num2str(DATA.T,2) '.mat'],"W_n1",'x','y','dV','DATA','glob')
% save(['Comi 3 contacts\Data\Melin 1994 – 256x256 grains - ' num2str(Nx) ' - precopr=' sprintf('%.2e',DATA.precompr-1) '_T=' num2str(DATA.T,2) '_Bad.mat'],"W_n1",'x','y','dV','DATA','glob')


%% Nested
function [] = global_fun(W,dV,t)
    global glob u v rho y gg e_c DATA a1 a2 a3 F11 F12 F21 F22
    if t == 0
        glob = [];  glob.ek = [];  glob.ep_gp = [];  glob.ep_c = [];  glob.ep_gg = [];  glob.m = [];  glob.t = [];
    end
    % I do not count the ghost cells as they change without conserving anything
    glob.ep_c(end+1) = dV.x(3:end-2) * squeeze(...  Potential energy stored in the contacts
                     + rho(W(:,3:end-2,3:end-2)).*e_c(W(:,3:end-2,3:end-2)) ...
                        ) * dV.y(3:end-2)'; 
    glob.ep_gg(end+1) = dV.x(3:end-2) * squeeze(...  Gravitationsl potential energy
                    gg*rho(W(:,3:end-2,3:end-2)).*reshape(y(3:end-2),1,1,[]) ...
                        ) * dV.y(3:end-2)'; 
    glob.ek(end+1) = 0.5*dV.x(3:end-2) * squeeze( ... Kinetic energy
                     rho(W(:,3:end-2,3:end-2)).* (u(W(:,3:end-2,3:end-2)).^2 + v(W(:,3:end-2,3:end-2)).^2) ...
                        ) * dV.y(3:end-2)';
    glob.m(end+1)  = dV.x(3:end-2) * squeeze(  rho(W(:,3:end-2,3:end-2))  ) *dV.y(3:end-2)'; % Total mass
    glob.t(end+1)  = t;  % Time

    %% Test of linear dependency of eigenvectors
        [RRA,~] = eigs_x_3_contacts(W(:,3:end-2,3:end-2),a1,a2,a3,rho,u,v,F11,F12,F21,F22,DATA); % Eigensystem of A evaluated in the interface orthogonal to x
        % [RRB,~] = f.eigB_tot(W(:,3:end-2,3:end-2)); % Eigensystem of B evaluated in the interface orthogonal to x
                
        RA = reshape(RRA(:,1:3,:,:),7,3,[]); % Verifico solo i primi 3 perché è solo lì che posso avere il problema.
        % Non uso tutta la matrice 7x7 perché capita che F = G_6 = G_7 = 0 (vedi notazione del paper)
        % RB = reshape(RRB,7,7,[]);
        
        rankA = arrayfun(@(k) rank(RA(:,:,k)), 1:size(RA,3));
        independentA = reshape(rankA == 3,size(RRA,3),size(RRA,4));
        % rankB = arrayfun(@(k) rank(RB(:,:,k)), 1:size(RB,3));
        % independentB = reshape(rankB == 7,size(RRB,3),size(RRB,4));
        
        % detA = arrayfun(@(k) det(RA(:,:,k)), 1:size(RA,3));
        % detA = reshape(detA,size(RRA,3),size(RRA,4));
        
        % if any(abs(detA) < 10,'all') % || any(~independentB,'all')
        if any(~independentA,'all') % || any(~independentB,'all')
            warning('Linearly dependent eigenvector columns detected.');
        end
end


% function A = A_fun_tot(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA)
%     % See my Matematica code: Hyperbolicity_Grains_F_no_beta.nb
%     A = zeros(size(W,1),size(W,1),size(W,2),size(W,3));
%     [u1,u2,kH,n,F11,F12,F21,F22,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,AA,BB,CC,DD,EE,FF,AAA,BBB,CCC,DDD,EEE,FFF] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
% 
% end


% function B = B_fun(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA)
%     % See my Matematica code: Hyperbolicity_Grains_F_no_beta.nb
%     B = zeros(size(W,1),size(W,1),size(W,2),size(W,3));
%     [u1,u2,kH,n,F11,F12,F21,F22,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,AA,BB,CC,DD,EE,FF,AAA,BBB,CCC,DDD,EEE,FFF] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
% 
% end


function [RRx,DDx] = eigA_fun_spl1(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA,restr)
    if exist('restr','var')
        sz = [length(restr), size(W,2), size(W,3)];
    else 
        sz = size(W);
    end

    RRx = zeros(sz(1),sz(1),sz(2),sz(3));
    DDx = zeros(sz(1),  1,  sz(2),sz(3));
%   [u1,u2,kH,n,F11,F12,F21,F22,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,xi1,xi2,xi3,AA,BB,CC,DD,EE,FF,AAA,BBB,CCC,DDD,EEE,FFF] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
    [u1,u2,~, ~,~,  ~,  ~,  ~,  ~, ~, b1,b2,~, ~, d1,d2,~, ~, f1,f2, ~,  ~,  ~, AA,~, CC, ~,EE, ~,AAA,~,  CCC,~,  EEE,~  ] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
        
    DDx(1,1,:,:) = u1;
    DDx(2,1,:,:) = u1 + sqrt(AAA + CCC + EEE);
    DDx(3,1,:,:) = u1 - sqrt(AAA + CCC + EEE);

    % RRx(1,1,:,:) = 0;
    % RRx(2,1,:,:) = 0;
    RRx(3,1,:,:) = 1;

    RRx(1,2,:,:) = 1;
    RRx(2,2,:,:) = u1 + sqrt(AAA + CCC + EEE);
    RRx(3,2,:,:) = u2 + (AA .*b1.*b2)./sqrt(AAA + CCC + EEE) + (CC .*d1.*d2)./sqrt(AAA + CCC + EEE) + (EE .*f1.*f2)./sqrt(AAA + CCC + EEE);

    RRx(1,3,:,:) = 1;
    RRx(2,3,:,:) = u1 - sqrt(AAA + CCC + EEE);
    RRx(3,3,:,:) = u2 - (AA .*b1.*b2)./sqrt(AAA + CCC + EEE) - (CC .*d1.*d2)./sqrt(AAA + CCC + EEE) - (EE .*f1.*f2)./sqrt(AAA + CCC + EEE);

end

function [RRy,DDy] = eigB_fun_spl1(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA,restr)
    if exist('restr','var')
        sz = [length(restr), size(W,2), size(W,3)];
    else 
        sz = size(W);
    end

    RRy = zeros(sz(1),sz(1),sz(2),sz(3));
    DDy = zeros(sz(1),1,sz(2),sz(3));

%   [u1,u2,kH,n,F11,F12,F21,F22,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,xi1,xi2,xi3,AA,BB,CC,DD,EE,FF,AAA,BBB,CCC,DDD,EEE,FFF,AAA2,BBB2,CCC2,DDD2,EEE2,FFF2] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
    [u1,u2,~, ~, ~,  ~,  ~,  ~,  ~, ~, b1,b2,~, ~,d1,d2,~, ~, f1,f2,~,  ~,  ~,  AA, ~,CC, ~,EE, ~, ~ , ~ , ~ , ~ , ~ , ~ ,AAA2, ~  ,CCC2, ~ , EEE2, ~  ] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
        
    
    DDy(1,1,:,:) = u2;
    DDy(2,1,:,:) = u2 + sqrt(AAA2 + CCC2 + EEE2);
    DDy(3,1,:,:) = u2 - sqrt(AAA2 + CCC2 + EEE2);
    % 
    % RRy(1,1,:,:) = 0;
    RRy(2,1,:,:) = 1;
    % RRy(3,1,:,:) = 0;

    RRy(1,2,:,:) = 1;
    RRy(2,2,:,:) = u1 + (AA .*b1.*b2)./sqrt(AAA2 + CCC2 + EEE2) + (CC .*d1.*d2)./sqrt(AAA2 + CCC2 + EEE2) + (EE .*f1.*f2)./sqrt(AAA2 + CCC2 + EEE2);
    RRy(3,2,:,:) = u2 + sqrt(AAA2 + CCC2 + EEE2);

    RRy(1,3,:,:) = 1;
    RRy(2,3,:,:) = u1 - (AA .*b1.*b2)./sqrt(AAA2 + CCC2 + EEE2) - (CC .*d1.*d2)./sqrt(AAA2 + CCC2 + EEE2) - (EE .*f1.*f2)./sqrt(AAA2 + CCC2 + EEE2);
    RRy(3,3,:,:) = u2 - sqrt(AAA2 + CCC2 + EEE2);

end

function [veloc] = veloc_fun(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA)

    veloc = zeros(6, size(W,2),size(W,3));
%   [u1,u2,kH,n,F11,F12,F21,F22,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,xi1,xi2,xi3,AA,BB,CC,DD,EE,FF,AAA,BBB,CCC,DDD,EEE,FFF,AAA2,BBB2,CCC2,DDD2,EEE2,FFF2] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
    [u1,u2,~, ~,~,  ~,  ~,  ~,  ~, ~,  ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~,  ~,  ~ , ~, ~, ~, ~, ~, ~,AAA,~,  CCC,~,  EEE,~  ,AAA2, ~  ,CCC2, ~ , EEE2, ~  ] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
        
    veloc(1,:,:) = u1;
    veloc(2,:,:) = u1 + sqrt(AAA + CCC + EEE);
    veloc(3,:,:) = u1 - sqrt(AAA + CCC + EEE);
    
    veloc(1,:,:) = u2;
    veloc(2,:,:) = u2 + sqrt(AAA2 + CCC2 + EEE2);
    veloc(3,:,:) = u2 - sqrt(AAA2 + CCC2 + EEE2);

end

function c = num_div(vec1,vec2)
% Questo calcola la divergenza calcolata nei centri cella dato il vettore [a1; a2]
global dV
    c = zeros(size(vec1(1,:,:)));
    dx = dV.x(1);   dy = dV.y(1);
    i = 2:length(dV.x)-1;  j = 2:length(dV.y)-1;
    c(1,i,j) = (vec1(1,i+1,j) - vec1(1,i-1,j))/2/dx  + (vec2(1,i,j+1) - vec2(1,i,j-1))/2/dy;

end

function out = dphix_dq_fun(W,xy,t,gx,fric)
    out = zeros(3,3,size(W,2),size(W,3));
    % There are many zeros
    out(2,1,:,:) = gx;
    out(2,2,:,:) = - fric(t);
end

function out = dphiy_dq_fun(W,xy,t,gy,fric)
    out = zeros(3,3,size(W,2),size(W,3));
    % There are many zeros
    out(3,1,:,:) = gy;
    out(3,3,:,:) = - fric(t);
end

