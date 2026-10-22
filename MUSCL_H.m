function [u_out] = MUSCL_H(u_in,dV,dt,f,veloc,BCs,t,x_sweep,noncons,salta_ultimi,restr)
% This is the MUSCL-Hancock scheme
% We follow the Toro 1997, Sec. 14.4
% The curl-free is applied to the n_vec that are stored at the end of the
% states (they occupy the last 2*n_vec spaces as they are 2D vectors)

if ~exist('x_sweep','var')
    x_sweep = true;
end
if exist('salta_ultimi','var') && ~isempty(salta_ultimi)
    u_in_pre_BCs = u_in;
else
    salta_ultimi = 0;
end

if ~exist('restr','var')  || isempty(restr)
    restr = 1:size(u_in,1);
end


% Boundary conditions with 2 ghost cells per side
u_in = app_BCs(u_in,BCs,t,x_sweep,2);
u_out = u_in;
if x_sweep
    Nc = length(dV.x);
    % ii = 2:Nx; % Queste solo le celle interne, le uniche che posso davvero propagare
    % dxx(i) è il dV medio attraverso l'interfaccia i+1/2 (cioè fra le celle i-1 e i)
    % dxx = (dV(jj) + dV(jj+1) )/2; % è definito solo a cavallo delle interfaces jj
    dxx = (dV.x(1:Nc-1) + dV.x(2:Nc) )/2; % è definito solo a cavallo delle interfaces jj
else
    Nc = length(dV.y);
    dxx = (dV.y(1:Nc-1) + dV.y(2:Nc) )/2; % è definito solo a cavallo delle interfaces jj
end


limiter = 2;
switch limiter 
    case 1   % Questo serve se seguo l'approccio del Toro
        w = 1/3;
    case 2
        beta = 2;
        % Nota che beta = 1 corrisponde a minmod
        %          beta = 2 corrisponde a superbee
    case 3
        % Parameters for the limiter
        % https://chimeracfd.com/programming/gryphon/muscl.html
        % phi = 0;  k = 1;   % first order, since phi = 0, the value of w does not count
        % phi = 1;  k = -1;  % Upwind 2nd order
        % phi = 1;  k = 0;   % Fromm  2nd order
        phi = 1;  k = 1/3; % 3rd order
        %%% Select a limiter, see     https://en.wikipedia.org/wiki/Flux_limiter
        % psi_lim = @(r) ( r + abs(r) )./(1+abs(r)); % Van Leer limiter
        % psi_lim = @(r) ( r.^2 + r )./(r.^2 + 1); % Van Albada limiter
        % psi_lim = @(r) max(0,min(1,r)); % minmod
        beta = 1.8; psi_lim = @(r) max(0,min(r,beta)); % Osher
        % psi_lim = @(r) ; % Another limiter
end


if x_sweep
    %%% Step I: Data reconstruction
    % I reconstruct the slopes inside each cells, only for the cells 2:Nc-1
    D_j = u_in(:,2:Nc,:) - u_in(:,1:Nc-1,:); % Questo è un vettore definito sulle Nc-1 internal interfaces
    D_i = zeros(size(u_in));  % Questo è un vettore definito nelle celle, i valori in 1 e in Nc non vengono usati
    switch limiter % Caso non limitato, tratto da Toro Eq. (14.37)
        case 1
            D_i(:,2:Nc-1,:) = 0.5*(1+w)*D_j(:,1:Nc-2,:) + 0.5*(1-w)*D_j(:,2:Nc-1,:);
            % Reconstruct the values on the intefaces for the cells 2:Nc-1
            u_L = u_in;   u_R = u_L; % I valori delle celle 1 e Nc non vengono utilizzati
            u_L(:,2:Nc-1,:) = u_in(:,2:Nc-1,:)-0.5*D_i(:,2:Nc-1,:);
            u_R(:,2:Nc-1,:) = u_in(:,2:Nc-1,:)+0.5*D_i(:,2:Nc-1,:);
        case 2 % with limted slope (vedi Toro, pag. 504 eq. (14.44) )
            D_i(:,2:Nc-1,:) = ...
                  max( max(...
                           min(beta*D_j(:,1:Nc-2,:),       D_j(:,2:Nc-1,:) )  , ...
                           min(     D_j(:,1:Nc-2,:),  beta*D_j(:,2:Nc-1,:) )) , ...
                           0   ) .*(D_j(:,2:Nc-1,:) > 0) ...
                + min( min(...
                           max(beta*D_j(:,1:Nc-2,:),       D_j(:,2:Nc-1,:) )  , ...
                           max(     D_j(:,1:Nc-2,:),  beta*D_j(:,2:Nc-1,:) )) , ...
                           0   ) .*(D_j(:,2:Nc-1,:) < 0);
            if isfield(BCs.L,'even')
                 % Il limiter alle cells a cavallo del boundary da problemi quando le ghost cells hanno lo stesso valore (sono even)
                 % In particolare, c'è una pendenza = 0 che forza due pendenze = 0 e il boundary diventa approssimato solo all'ordine 1
                 % Quindi in quelle celle scelgo la pendenza delle interfaccie lontane dal boundary
                D_i(BCs.L.even,3,:) =  D_j(BCs.L.even,3,:); % Ricorda che D_j vive sulle interfacce
                D_i(BCs.L.even,2,:) = -D_j(BCs.L.even,3,:); % = +D_j(BCs.B.even,1,:)
            end
            if isfield(BCs.R,'even')
                D_i(BCs.R.even,Nc-2,:) =  D_j(BCs.R.even,Nc-3,:);
                D_i(BCs.R.even,Nc-1,:) = -D_j(BCs.R.even,Nc-3,:); % = +D_j(BCs.B.even,Nc-1,:)
            end
            % Reconstruct the values on the intefaces for the cells 2:Nc-1
            u_L = u_in;   u_R = u_L; % I valori delle celle 1 e Nc non vengono utilizzati
            u_L(:,2:Nc-1,:) = u_in(:,2:Nc-1,:)-0.5*D_i(:,2:Nc-1,:);
            u_R(:,2:Nc-1,:) = u_in(:,2:Nc-1,:)+0.5*D_i(:,2:Nc-1,:);
        case 3
            % r is the ratio between consecutive slopes (from L to R = D_R/D_L) defined inside the cells 2:Nc-1
            r = zeros(size(u_in)); psi_LR = r; psi_RL = r;
            r(:,2:Nc-1,:) = D_j(:,2:Nc-1,:) ./ D_j(:,1:Nc-2,:);
            r(isnan(r)) = 1; % If any 0/0 occurred
            r(r>100) = 100; % To avoid NaN due to inf / inf
            r(r<-100) = -100;
            r(r==0) = 10*eps; % To avoid problems when inverting
            
            % Van Leer limiter defined inside the cells 2:Nc-1,
            % either for L-to-R or R-to-L (which means taking the inverse)
            psi_LR(:,2:Nc-1,:) = psi_lim(r(:,2:Nc-1,:));
            psi_RL(:,2:Nc-1,:) = psi_lim(1./r(:,2:Nc-1,:));

            % Reconstruct the values on the intefaces for the cells 2:Nc-1
            u_L = u_in;   u_R = u_L; % L'unico modo che ho per ricostruire i valori delle celle 1 e Nc è quello di copiarli (questi valori servono solo se c'è la parte nc che contiene delle derivate)
            u_L(:,2:Nc-1,:) = u_in(:,2:Nc-1,:) -phi/4*((1+k)*psi_LR(:,2:Nc-1,:).*D_j(:,1:Nc-2,:) +(1-k)*psi_RL(:,2:Nc-1,:).*D_j(:,2:Nc-1,:));
            u_R(:,2:Nc-1,:) = u_in(:,2:Nc-1,:) +phi/4*((1-k)*psi_LR(:,2:Nc-1,:).*D_j(:,1:Nc-2,:) +(1+k)*psi_RL(:,2:Nc-1,:).*D_j(:,2:Nc-1,:));
    end
    % % Dove la denstà è negativa, metto la pendenza = 0
    % dove_neg = u_L(1,:,:) < 200 | u_R(1,:,:) < 200;
    % u_L(dove_neg) = u_in(dove_neg);
    % u_R(dove_neg) = u_in(dove_neg);
    
    %%% Step II: Evolution inside the cells 2:Nc-1 for half time step  dt/2
    uu_L = u_L; uu_R = u_R; evol = zeros(size(uu_L(restr,:,:)));
    i = 2:Nc-1; % Ricorda che quando valuto i flussi in i, in realtà li sto valutando in i+1/2
    F_L = f.x(u_L(:,i,:),dV.xy(:,i,:));   % attento a non chiamare f con valori irrealistici (e.g. 0) perché funzioni tipo alpha_p sballerebbero
    F_R = f.x(u_R(:,i,:),dV.xy(:,i,:));
    evol(:,i,:) = 0.5*dt./dV.x(i) .*( F_R - F_L ); % Questo è l'evoluzione (senza stabilizzazione -> solo predizione)
    if noncons % Se esiste, aggiungo anche il contributo non conservativo
        [Bnc,FncL,FncR] = f.ncx(u_in,uu_L,uu_R);  % Qui, L e R sono intesi come L e R rispetto alle facce perché la propagazione è interna
        evol(:,i,:) = evol(:,i,:) + 0.5*dt./dV.x(i)  .* Bnc(:,i,:).*( FncR(:,i,:) - FncL(:,i,:) );
    end
    uu_L(restr,i,:) = u_L(restr,i,:) - evol(:,i,:); % Indicato con \bar U^L sul Toro
    uu_R(restr,i,:) = u_R(restr,i,:) - evol(:,i,:);


    
    
    %%% Step III: Riemann problem
    % I solve this only for the cells 3:Nc-2
    % Let's use a Rusanov solver
    vel = abs(veloc.x(u_in));
    S = max([vel(:,1:Nc-1,:); vel(:,2:Nc,:)],[],1);
    % F_j is the flux across the interfaces 2:Nc-2
    %   that is the interfs. at the L of the cells 3:Nc-1 or at the R of the cells 2:Nc-2 (I have no info on uu_L and uu_R in 1 and Nc)
    
    % F_j = zeros(size(u_in,1), Nc-1, size(u_in,3));
    F_j = zeros(length(restr), Nc-1, size(u_in,3));
    f_L = f.x(uu_L(:,3:Nc-1,:),dV.xy(:,3:Nc-1,:));
    f_R = f.x(uu_R(:,2:Nc-2,:),dV.xy(:,2:Nc-2,:));
    F_j(:,2:Nc-2,:) = 0.5* (f_L + f_R) - 0.5*S(1,2:Nc-2,:) .* (uu_L(restr,3:Nc-1,:)-uu_R(restr,2:Nc-2,:));
    if noncons % If there are non conservative terms, I add the path integral to the flux
        % In questo caso il riferimento sono le interfaces, quindi 
        % uu_L è a Destra dell'interfaccia e uu_R è a sinistra dell'interfaccia
        [Bnc,FncL,FncR] = f.ncx(u_in,uu_R(:,1:Nc-1,:),uu_L(:,2:Nc,:));  % Qui, L e R sono intesi come L e R rispetto alle interfacce
        Fnc_i = 1/2 * (FncL + FncR);
        i = 3:Nc-2;
        u_out(restr,i,:) = u_in(restr,i,:) - dt*(  ( F_j(:,i,:)   + Bnc(:,i,:).*Fnc_i(:,i,:)   )./reshape(dxx(i)  , 1,Nc-4,[]) ...
                                                  -( F_j(:,i-1,:) + Bnc(:,i,:).*Fnc_i(:,i-1,:)   )./reshape(dxx(i-1), 1,Nc-4,[])  );
    else
        i = 3:Nc-2;
        u_out(restr,i,:) = u_in(restr,i,:) - dt*(  F_j(:,i,:)  ./reshape(dxx(i)  , 1,Nc-4,[]) ...
                                                  -F_j(:,i-1,:)./reshape(dxx(i-1), 1,Nc-4,[])  );
    end

else
    %%% Step I: Data reconstruction
    % I reconstruct the slopes inside each cells, only for the cells 2:Nc-1
    D_j = u_in(:,:,2:Nc) - u_in(:,:,1:Nc-1); % Questo è un vettore definito sulle Nc-1 internal interfaces
    D_i = zeros(size(u_in));  % Questo è un vettore definito nelle celle, i valori in 1 e in Nc non vengono usati
    switch limiter  
        case 1% Caso non limitato
            D_i(:,:,2:Nc-1) = 0.5*(1+w)*D_j(:,:,1:Nc-2) + 0.5*(1-w)*D_j(:,:,2:Nc-1);
            % Reconstruct the values on the intefaces for the cells 2:end-1
            u_L = u_in;   u_R = u_L; % I valori delle celle 1 e Nc non vengono utilizzati
            u_L(:,:,2:Nc-1) = u_in(:,:,2:Nc-1)-0.5*D_i(:,:,2:Nc-1);
            u_R(:,:,2:Nc-1) = u_in(:,:,2:Nc-1)+0.5*D_i(:,:,2:Nc-1);
        case 2 % with limted slope
            D_i(:,:,2:Nc-1) = ...
                  max( max(...
                           min(beta*D_j(:,:,1:Nc-2),       D_j(:,:,2:Nc-1) )  , ...
                           min(     D_j(:,:,1:Nc-2),  beta*D_j(:,:,2:Nc-1) )) , ...
                           0   ) .*(D_j(:,:,2:Nc-1) > 0) ...
                + min( min(...
                           max(beta*D_j(:,:,1:Nc-2),       D_j(:,:,2:Nc-1) )  , ...
                           max(     D_j(:,:,1:Nc-2),  beta*D_j(:,:,2:Nc-1) )) , ...
                           0   ) .*(D_j(:,:,2:Nc-1) < 0);
            if isfield(BCs.B,'even')
                D_i(BCs.B.even,:,3) =  D_j(BCs.B.even,:,3); % Ricorda che D_j vive sulle interfacce
                D_i(BCs.B.even,:,2) = -D_j(BCs.B.even,:,3); % = +D_j(BCs.B.even,:,1)
            end
            if isfield(BCs.T,'even')
                D_i(BCs.T.even,:,Nc-2) =  D_j(BCs.T.even,:,Nc-3).*(D_j(BCs.T.even,:,Nc-3) < 0) ;
                D_i(BCs.T.even,:,Nc-1) = -D_i(BCs.T.even,:,Nc-2); % -D_j(BCs.T.even,:,Nc-3); % = +D_j(BCs.B.even,:,Nc-1)
            end
            % Reconstruct the values on the intefaces for the cells 2:end-1
            u_L = u_in;   u_R = u_L; % I valori delle celle 1 e Nc non vengono utilizzati
            u_L(:,:,2:Nc-1) = u_in(:,:,2:Nc-1)-0.5*D_i(:,:,2:Nc-1);
            u_R(:,:,2:Nc-1) = u_in(:,:,2:Nc-1)+0.5*D_i(:,:,2:Nc-1);
        case 3
            % r is the ratio between consecutive slopes (from L to R = D_R/D_L) defined inside the cells 2:Nc-1
            r = zeros(size(u_in)); psi_LR = r; psi_RL = r;
            r(:,:,2:Nc-1) = D_j(:,:,2:Nc-1) ./ D_j(:,:,1:Nc-2);
            r(isnan(r)) = 1; % If any 0/0 occurred
            r(r>100) = 100; % To avoid NaN due to inf / inf
            r(r<-100) = -100;
            r(r==0) = 10*eps; % To avoid problems when inverting
            
            % Van Leer limiter defined inside the cells 2:Nc-1,
            % either for L-to-R or R-to-L (which means taking the inverse)
            psi_LR(:,:,2:Nc-1) = psi_lim(    r(:,:,2:Nc-1));
            psi_RL(:,:,2:Nc-1) = psi_lim( 1./r(:,:,2:Nc-1));
            % Reconstruct the values on the intefaces for the cells 2:Nc-1
            u_L = u_in;   u_R = u_L; % L'unico modo che ho per ricostruire i valori delle celle 1 e Nc è quello di copiarli (questi valori servono solo se c'è la parte nc che contiene delle derivate)
            u_L(:,:,2:Nc-1) = u_in(:,:,2:Nc-1) -phi/4*((1+k)*psi_LR(:,:,2:Nc-1).*D_j(:,:,1:Nc-2) +(1-k)*psi_RL(:,:,2:Nc-1).*D_j(:,:,2:Nc-1));
            u_R(:,:,2:Nc-1) = u_in(:,:,2:Nc-1) +phi/4*((1-k)*psi_LR(:,:,2:Nc-1).*D_j(:,:,1:Nc-2) +(1+k)*psi_RL(:,:,2:Nc-1).*D_j(:,:,2:Nc-1));
    end
    % % Dove la denstà è negativa, metto la pendenza = 0
    % dove_neg = u_L(1,:,:) < 200 | u_R(1,:,:) < 200;
    % u_L(dove_neg) = u_in(dove_neg);
    % u_R(dove_neg) = u_in(dove_neg);
    
    %%% Step II: Evolution inside the cells 2:end-1
    uu_L = u_L; uu_R = u_R; zeros(size(uu_L));  evol = zeros(size(uu_L(restr,:,:)));
    i = 2:Nc-1; % Ricorda che quando valuto i flussi in i, in realtà li sto valutando in i+1/2
    F_L = f.y(u_L(:,:,i),dV.xy(:,:,i));
    F_R = f.y(u_R(:,:,i),dV.xy(:,:,i));
    evol(:,:,i) = 0.5*dt./reshape(dV.y(2:Nc-1),1,1,[]) .*( F_R - F_L );
    if noncons % Se esiste, aggiungo anche il contributo non conservativo
        [Bnc,FncL,FncR] = f.ncy(u_in,uu_L,uu_R);  % Qui, L e R sono intesi come L e R rispetto alle facce perché la propagazione è interna
        evol(:,:,i) = evol(:,:,i) + 0.5*dt./reshape(dV.y(i),1,1,[])  .* Bnc(:,:,i).*( FncR(:,:,i) - FncL(:,:,i) );
    end
    uu_L(restr,:,i) = u_L(restr,:,i) - evol(:,:,i); % Indicato con \bar U^L sul Toro
    uu_R(restr,:,i) = u_R(restr,:,i) - evol(:,:,i);

    
    %%% Step III: Riemann problem
    % I solve this only for the cells 3:end-2
    % Let's use a Rusanov solver
    vel = abs(veloc.x(u_in));
    S = max([vel(:,:,1:Nc-1); vel(:,:,2:Nc)],[],1);
    % F_j is the flux across the interfaces 2:Nc-2
    %   that is the interfs. at the L of the cells 3:Nc-1 nor at the R of the cells 2:Nc-2 (I have no info on uu_L and uu_R in 1 and Nc)
    F_j = 0*u_in(restr,:,1:Nc-1);
    f_L = f.y(uu_L(:,:,3:Nc-1),dV.xy(:,:,3:Nc-1));   f_R = f.y(uu_R(:,:,2:Nc-2),dV.xy(:,:,2:Nc-2));
    F_j(:,:,2:Nc-2) = 0.5* (f_L + f_R) - 0.5*S(1,:,2:Nc-2) .* (uu_L(restr,:,3:Nc-1)-uu_R(restr,:,2:Nc-2));
    if noncons % If the are non conservative terms, I add the path integral to the flux
        % In questo caso il riferimento sono le interfaces, quindi 
        % uu_L è a Destra dell'interfaccia e uu_R è a sinistra dell'interfaccia
        [Bnc,FncL,FncR] = f.ncy(u_in,uu_R(:,:,1:Nc-1),uu_L(:,:,2:Nc));  % Qui, L e R sono intesi come L e R rispetto alle interfacce
        Fnc_i = 1/2 * (FncL + FncR);
        i = 3:Nc-2;
        u_out(restr,:,i) = u_in(restr,:,i) - dt*(  ( F_j(:,:,i)   + Bnc(:,:,i).*Fnc_i(:,:,i)   )./reshape(dxx(i)  , 1,1,Nc-4) ...
                                                  -( F_j(:,:,i-1) + Bnc(:,:,i).*Fnc_i(:,:,i-1) )./reshape(dxx(i-1), 1,1,Nc-4)  );
    else
        i = 3:Nc-2;
        u_out(restr,:,i) = u_in(restr,:,i) - dt*(  F_j(:,:,i)  ./reshape(dxx(i-1), 1,1,Nc-4) ...
                                                  -F_j(:,:,i-1)./reshape(dxx(i-1), 1,1,Nc-4) );
    end
end

if salta_ultimi
    % I do not want to propagate the curl-free vectors now
    u_out(end-salta_ultimi+1:end,:,:) = u_in_pre_BCs(end-salta_ultimi+1:end,:,:);
end