function [u_n1] = Rusanov(u_n,dV,dt,f,veloc,BCs,t,x_sweep,noncons,restr)
% restr (restriction) is an index vector that contains the indexes of the state to
% propagate with the Rusanov. This is useful for the splitting of the
% hyperbolic system with respect to the state.
% Accordingly, the flux function must be provided as a vector long
% length(restriction).
if ~exist('x_sweep','var')
    x_sweep = true;
end

if ~exist('noncons','var')
    noncons = false;
end
if ~exist('restr','var') || isempty(restr)
    restr = 1:size(u_n,1);
end

% Boundary conditions
u_n = app_BCs(u_n,BCs,t,x_sweep);
u_n1 = u_n;

% u_n1 = u_n + dt/dx (F_i-1/2 - F_i+1/2)
% F_i+1/2 = 1/2 (F_i + F_i+1) + 1/2 S^+ (U_i - U_i+1)
% Le celle interne, le uniche che posso davvero propagare, sono numerate 2:Nc-1
% Le interfacce interne, le uniche che posso davvero usare, sono numerate 1:Nc-1

% S = squeeze(zeros(size(u_n(1,:,:))));
% S = zeros(1,Nx+1);  % Questo è l'OLD per 1D
% S(1:end-1) = max([vel(:,1:end-1); vel(:,2:end)],[],1);
% S(i) è la S attraverso l'interfaccia i+1/2
% F(i) è il flusso attraverso l'interfaccia i+1/2
% dxx(i) è il dV medio attraverso l'interfaccia i+1/2 (cioè fra le celle i e i+1)
if x_sweep
    Nc = length(dV.x);
    dxx = (dV.x(1:Nc-1) + dV.x(2:Nc) )/2; % è definito solo a cavallo delle interfaces 1:Nc-1
    vel = abs(veloc.x(u_n));
    S = max([vel(:,1:Nc-1,:); vel(:,2:Nc,:)],[],1);
    % F_j is the flux across the interfaces 1:Nc-1
    %   that is the interfs. at the L of the cells 2:Nc or at the R of the cells 1:Nc-1
    % F_j = 0*u_n(:,1:Nc-1,:);
    F = f.x(u_n,dV.xy);
    F_i = 1/2 * (F(:,1:Nc-1,:) + F(:,2:Nc,:)) - 1/2*S(1,1:Nc-1,:) .* (u_n(restr,2:Nc,:)-u_n(restr,1:Nc-1,:));
    if noncons % If the are non conservative terms, I add the path integral to the flux
        i = 2:Nc-1; % Quando valuto i flussi in i, in realtà li sto valutando in i+1/2
        [Bnc,FncL,FncR] = f.ncx(u_n,u_n(:,1:Nc-1,:),u_n(:,2:Nc,:));  % L e R sono intesi come L e R rispetto alle interfacce
        Fnc_i = 1/2 * (FncL + FncR); % - 1/2*S(1,1:Nc-1,:) .* (??u_n(restr,2:Nc,:)-??u_n(restr,1:Nc-1,:));  % Forse non devo introdurre dissipazione perché ce'è già nel flusso conservativo
        % Bnc dovrebbe essere una matrice che si applica a tutto lo stato, ma per ora mi basta che sia un
        % vettore quindi faccio il prodotto Bnc.*Fnc e non tensorprod(Bnc,Fnc)
        u_n1(restr,i,:) = u_n(restr,i,:) - dt*(  ( F_i(:,i,:)    + Bnc(:,i,:).*Fnc_i(:,i,:)  ) ./reshape(dxx(i)  , 1,Nc-2,[]) ...
                                                -( F_i(:,i-1,:)  + Bnc(:,i,:).*Fnc_i(:,i-1,:)) ./reshape(dxx(i-1), 1,Nc-2,[])  );
    else
        i = 2:Nc-1;
        u_n1(restr,i,:) = u_n(restr,i,:) - dt*(  F_i(:,i,:)  ./reshape(dxx(i)  , 1,Nc-2,[]) ...
                                                -F_i(:,i-1,:)./reshape(dxx(i-1), 1,Nc-2,[]) );
    end

    if isfield(f,'Sx') % Questo se tratto alcuni termini di source dentro lo schema
        warning('Non ho implementato restriction')
        Sx = f.Sx(u_n);
        u_n1(:,i,:) = u_n1(:,i,:) + dt*Sx(:,i,:);  % Aggiungo il source locale
        % u_n1(:,i,:) = u_n1(:,i,:) + dt*(Sx(:,i-1,:) +2*Sx(:,i,:) +Sx(:,i+1,:))/4;  % Aggiungo il source locale
    end

else
    Nc = length(dV.y);
    dxx = (dV.y(1:Nc-1) + dV.y(2:Nc) )/2; % è definito solo a cavallo delle interfaces 1:Nc-1
    vel = abs(veloc.y(u_n));
    S = max([vel(:,:,1:Nc-1); vel(:,:,2:Nc)],[],1);
    F = f.y(u_n,dV.xy);
    F_i = 1/2 * (F(:,:,1:Nc-1) + F(:,:,2:Nc)) - 1/2*S(1,:,1:Nc-1) .* (u_n(restr,:,2:Nc)-u_n(restr,:,1:Nc-1));
    % F_im1_2 = 1/2 * (f(u_n(restr,:,2:Nc-1-1)) + f(u_n(restr,:,2:Nc-1))) - reshape(1/2*S(:,2:Nc-1-1), 1,[],Nx-1) .* (u_n(restr,:,2:Nc-1)-u_n(restr,:,2:Nc-1-1));
    if noncons % If the are non conservative terms, I add the path integral to the flux
        i = 2:Nc-1; % Quando valuto i flussi in i, in realtà li sto valutando in i+1/2
        [Bnc,FncL,FncR] = f.ncy(u_n,u_n(:,:,1:Nc-1),u_n(:,:,2:Nc));
        Fnc_i = 1/2 * (FncL + FncR); % - 1/2*S(1,1:Nc-1,:) .* (??u_n(restr,2:Nc,:)-??u_n(restr,1:Nc-1,:));  % Forse non devo introdurre dissipazione perché ce'è già nel flusso conservativo
        u_n1(restr,:,i) = u_n(restr,:,i) - dt*(  ( F_i(:,:,i)    + Bnc(:,:,i).*Fnc_i(:,:,i)  ) ./reshape(dxx(i)  , 1,[],Nc-2) ...
                                                -( F_i(:,:,i-1)  + Bnc(:,:,i).*Fnc_i(:,:,i-1)) ./reshape(dxx(i-1), 1,[],Nc-2)  );
    else
        i = 2:Nc-1;
        u_n1(restr,:,i) = u_n(restr,:,i) - dt*(  F_i(:,:,i)   ./reshape(dxx(i)  , 1,[],Nc-2) ...
                                                -F_i(:,:,i-1) ./reshape(dxx(i-1), 1,[],Nc-2) );
    end

    if isfield(f,'Sy') % Questo se tratto alcuni termini di source dentro lo schema
        warning('Non ho implementato restriction')
        Sy = f.Sy(u_n);
        u_n1(:,:,i) = u_n1(:,:,i) + dt*Sy(:,:,i) ;  % Aggiungo il source locale
    end
end

% Questo per il caso Euleriano dove la media sx e dx sono uguali
% dxx = (dV(2:Nc-1-1)/2 + dV(2:Nc-1) + dV(2:Nc-1+1)/2 )/2; % dxx è definito solo nelle celle 2:Nc-1
% u_n1(:,2:Nc-1) = u_n(restr,2:Nc-1) + dt*(F_im1_2 - F_ip1_2)./dxx;