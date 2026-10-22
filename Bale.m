function [u_n1] = Bale(u_n,dV,dt,f,veloc,BCs,t,restr)
% This implementation follows the PhD thesis of Bale 2002.
% The method is also explained in Leveque 2002, chapter 17.14
% In short, the procedure is:
% 1. At each interface, compute the right and left eigenvectors of the matrix A
%     OLD: I compute the eigenvalues of such a matrix in each cell and take the
%          left/right eigenvalues according to the interface. This should be more stable than
%          computing the eigenvalues of the "average matrix", computed on the average of the state vectors between two cells.
%     Unfortunately, I could end up in more eigs than I need... so I use the average of the state and cross the fingers
% 2. Compute the flux difference at the interface i-1/2   (3.12) (or (4.42) in case of source)
% 3. Decompose the flux in f-waves Z^p_{i-1/2}            (3.12)
% 4. Compute the fluctuations A^{+-} \Delta Q_{i-1/2} = sum_{p\in P^{+-}} Z^p_{i-1/2}           (3.22)
% 5. Propagate using     (2.31) (or better (3.26))



if ~exist('restr','var') || isempty(restr) % For the definition of this variable, see Rusanov.m
    restr = 1:size(u_n,1);
end

sz = size(u_n(restr,:,:));

if length(sz) == 2
%% Questo è il caso 1D

    % Boundary conditions with 2 ghost cells per side
    u_n = app_BCs(u_n,BCs,t,true,2);
    u_n1 = u_n;

    Nc = sz(2);
    i = 1:Nc-1;
    dxx = (dV.x(i) + dV.x(i+1) )/2; % è definito solo a cavallo delle interfaces 1:Nc-1
    dxx = reshape(dxx,1,1,Nc-1);
    u_interf = ( u_n(:,i) + u_n(:,i+1) )/2; % Questi sono i valori medi nelle interfacce i+1/2
    A_int_x = f.A(u_interf);
    % Oppure
    % A = (f.A(u_n(:,i)) + f.A(u_n(:,i+1)) ) / 2;

        % Calcolo gli autovalori 
        % RR = zeros(sz(1),sz(1),Nc-1);
        % LL = zeros(sz(1),sz(1),Nc-1);
        % DD = zeros(sz(1),Nc-1);
        % for k = 1:Nc-1 % Calcolo tutti gli eigenvectors nelle interfacce
        %     [V,D] = eig(A(:,:,k),'vector');  % V right eigenvectors, D eigenvalues, W left eigenvectors
        %     % the columns of W are the corresponding left eigenvectors, so that W'*A = D*W'
        %     RR(:,:,k) = V;
        %     LL(:,:,k) = inv(V); % Sono obbligato a fre l'inversa perché i left-eigs che mi restituirebbe eig(A) sono normalizzati diversamente
        %     DD(:,k) = D;
        % end
    % In compact form:
    [RR,DD] = pageeig(A_int_x,'vector');

    F = f.x(u_n,dV.xy);
    DF = F(:,i+1) - F(:,i) ;
    if isfield(f,'phi')
        PHI = f.phi(u_n,dV.xy,t) .* dV.x; % Questo è il source in ogni cella, concentrato in delta di Dirac
        PHI = app_BCs(PHI,BCs.BCsPHI,t,true,2);      % Visto che questo schema considera anche il RHS delle ghost cells, 
        %                                            %  devo applicare la simmetria anche lì!!!!!!!!!
        %                                            % Almeno per un hard wall...
        DF = DF - ( PHI(restr,i+1) + PHI(restr,i) )/2;   % Definisco il source all'interfaccia come la media dei source nelle celle
    end
    %%% Faccio la proiezione sulle f-waves
        % beta = zeros(sz(1),sz(2)-1);
        % Z = zeros(sz(1),sz(2)-1,sz(2)-1);
        % for k = 1:sz(2)-1
        %     beta(:,k) = LL(:,:,k)' * DF(:,k);   % Eq. (3.13): i left-autovettori, che sono le righe di LL, sono contratti con DF
        %     Z(:,:,k) = RR(:,:,kk)  .* beta(:,k)' ;  % Ogni colonna di RR è moltiplicata da il corrispondente coefficiente beta
        % end
    % In compact form:
    DF = reshape(DF, sz(1), 1, Nc-1);
        % Se volessi usare i left eigenvectors
        % beta = squeeze( pagemtimes(LL, DF3)  ); % Ricorda che LL ha i left-eigenvectors sulle righe
        % LL = pageinv(RR); % Questo produce LL con i left-eigenvectors sulle righe
    beta = pagemldivide(RR,DF);
    Z = RR .* pagetranspose(beta);

    %%% Costruisco le variazioni
        % ApDQ = zeros(sz(1),sz(2)-1); % This is A^+ \Delta Q_{i-1/2}
        % AmDQ = zeros(sz(1),sz(2)-1); % A^-
        % ApDQ(:,k) = Z(:,:,k) .* (D(:,k)>0)'; % Ogni colonna di Z è moltiplicata per 1 o 0, a seconda se corrisponde a un autovalore positivo o negativo, rispettivamente.
        % AmDQ(:,k) = Z(:,:,k) .* (D(:,k)<0)'; % idem ma prendo solo i negativi
    % In compact form:
    ApDQ = squeeze(pagemtimes(Z, double(DD>0) ));
    AmDQ = squeeze(pagemtimes(Z, double(DD<0) ));
    
    %%% Applico il limiter a Z_{i+-1/2}
    % Uso il metodo dopo Eq. (22), a fondo pagina, del paper 1997 Leveque - Wave propagation algorithms for hyperbolic systems
    i = 2:Nc-2;
    theta = zeros(1,sz(1),Nc-1);
    for p = 1:sz(1)   %  Z = sz(1) x sz(1) x Nc-1
        i_uw = i(DD(p,i) > 0); %  Mi confronto con le celle upwind
        i_dw = i(DD(p,i) < 0);
        if any(i_uw)
            theta(1,p,i_uw) =   pagemtimes(Z(:,p,i_uw),'transpose', Z(:,p,i_uw-1),'none') ...
                             ./ pagemtimes(Z(:,p,i_uw),'transpose', Z(:,p,i_uw)  ,'none');
        end
        if any(i_dw)
            theta(1,p,i_dw) =   pagemtimes(Z(:,p,i_dw),'transpose', Z(:,p,i_dw+1),'none') ...
                             ./ pagemtimes(Z(:,p,i_dw),'transpose', Z(:,p,i_dw)  ,'none');
        end
    end
    % pphi = max(0,min(1,theta));                     % minmod
    % pphi = max(0,max(min(1,2*theta),min(2,theta))); % Superbee
    pphi = max(0, min(min((1+theta)/2, 2), 2*theta));  % monotonized centered
    Z_lim = Z .* pphi;
    % Z_lim = Z;
    % Z_lim = Z*0;

    %% Calcolo la correzione dei flussi F_{i+-1/2}
    coeff = 0.5* sign(DD).*( 1 - dt * abs(DD)./reshape(dxx,1,1,Nc-1,1) ); %  sz(1) x 1 x Nc-1
    F_i = pagemtimes(  Z_lim, coeff );
    F_i = squeeze(F_i);
    %% Faccio la propagazione Eq. (3.26)
    i = 3:Nc-2;
    u_n1(restr,i) = u_n(restr,i) - dt*(   AmDQ(:,i)   ./dxx(1,i)  ...
                                        + ApDQ(:,i-1) ./dxx(1,i-1) ) ...
                                 - dt*(    F_i(:,i)   ./dxx(1,i)   ...
                                        -  F_i(:,i-1) ./dxx(1,i-1) );
    if isfield(f,'phi')
        dphi_dq = f.dphi_dq(u_n,dV.xy,t);
        dxx = reshape(dxx,1,1,Nc-1);
        u_n1(restr,i) = u_n1(restr,i) - dt^2/4* squeeze(  pagemtimes(dphi_dq(:,:,i), sum( ... terms to add that depends on dphi/dW. See Eq. (4.45) of Bale 2002
                                                          Z_lim(:,:,i)   ./dxx(1,1,i) ...
                                                        + Z_lim(:,:,i-1) ./dxx(1,1,i-1)...
                                                       ,2)  ))       ; 
    end

else
%%% Caso 2D
    % Boundary conditions with 2 ghost cells per side
    u_n = app_BCs(u_n,BCs,t,true,2);  % Along x
    u_n = app_BCs(u_n,BCs,t,false,2); % Along y
    u_n1 = u_n;

    %% Inizialisation along x
    Ncx = sz(2);
    i = 1:Ncx-1;
    dxx = (dV.x(i) + dV.x(i+1) )/2; % è definito solo a cavallo delle interfaces 1:Ncx-1
    dxx = reshape(dxx,1,Ncx-1,1);

    u_interf_x = ( u_n(:,i,:) + u_n(:,i+1,:) )/2; % Questi sono i valori medi nelle interfacce i+1/2
    if ~isfield(f,'eigA')
        A_int_x = f.A(u_interf_x);
        B_int_x = f.B(u_interf_x);
        % Calcolo gli autovalori
        [RA_int_x,DA_int_x] = pageeig(A_int_x,'balance','vector');
        [RB_int_x,DB_int_x] = pageeig(B_int_x,'balance','vector');
    else
        [RA_int_x,DA_int_x] = f.eigA(u_interf_x); % Eigensystem of A evaluated in the interface orthogonal to x
        [RB_int_x,DB_int_x] = f.eigB(u_interf_x); % Eigensystem of B evaluated in the interface orthogonal to x
    end

    %% Inizialisation along y
    Ncy = sz(3);
    j = 1:Ncy-1;
    dyy = (dV.y(j) + dV.y(j+1) )/2; % è definito solo a cavallo delle interfaces 1:Ncy-1
    dyy = reshape(dyy,1,1,Ncy-1);

    u_interf_y = ( u_n(:,:,j) + u_n(:,:,j+1) )/2; % Questi sono i valori medi nelle interfacce j+1/2
    if ~isfield(f,'eigB')
        B_int_y = f.B(u_interf_y);
        A_int_y = f.A(u_interf_y);
        % Calcolo gli autovalori
        [RB_int_y,DB_int_y] = pageeig(B_int_y,'balance','vector');
        [RA_int_y,DA_int_y] = pageeig(A_int_y,'balance','vector');
    else
        [RB_int_y,DB_int_y] = f.eigB(u_interf_y); % Eigensystem of B evaluated in the interface orthogonal to y
        [RA_int_y,DA_int_y] = f.eigA(u_interf_y); % Eigensystem of A evaluated in the interface orthogonal to y
    end

    %% purtroppo a volte non trova tutti gli autovettori, quindi faccio un check
    if ((~isfield(f,'eigA')) || (~isfield(f,'eigB')) ) && (~isfield(f,'fix_eigs'))
        for ii = 1:Ncx
            for jj = 1:Ncy
                switch 10  %%%%%%%% Alla fine il metodo 1 funzionicchia per le Euler equations, ma non molto oltre...
                    case 1
                    %%% Metodo 1: se il numero di condizionamente è alto, aggiungo eps e ci riprovo
                        if ii~=Ncx && (~isfield(f,'eigA'))
                            if cond(RA_int_x(:,:,ii,jj)) > 100
                                [RA_int_x(:,:,ii,jj),DA_int_x(:,:,ii,jj)] = eig(squeeze(A_int_x(:,:,ii,jj))+eps,'balance','vector');
                            end
                        end
                        if jj~=Ncy && (~isfield(f,'eigB'))
                            if cond(RB_int_y(:,:,ii,jj)) > 100
                                [RB_int_y(:,:,ii,jj),DB_int_y(:,:,ii,jj)] = eig(squeeze(B_int_y(:,:,ii,jj))+eps,'balance','vector');
                            end
                        end
                    case 2
                    %%% Metodo 2: ci provo tante volte finché il numero di cond. non si abbassa
                        % kk = 1e-4;
                        % while cond(RRy(:,:,ii,jj)) > 500 || ~isreal(RRy(:,:,ii,jj))
                        %     [RRytmp,DDytmp] = eig(squeeze(B(:,:,ii,jj)+eps*kk),'nobalance','vector');
                        %     RRy(:,:,ii,jj) = RRytmp;
                        %     DDy(:,:,ii,jj) = DDytmp;
                        %     kk = 10*kk;
                        % end
                        % if kk > 1e2
                        %     warning('ss')
                        % end
                    case 3
                    %%% Metodo 3: se cond. grande, cambio di eps lo stato, ricalcolo le matrici e gli autov.
                        error('Per questo vedere il fix_eigs qui sotto')
                    case 4
                    %%% Metodo 4: gli eigs che fanno casino sono quelli nulli. Questi non contribuiscono al flusso.
                    %             Numericamente trovo sempre gli autovalori, il problema sono gli autovettori.
                    %             Quindi individuo gli autovalori = 0 e ricostruisco i relativi autovettori in modo che completino il null-space delle matrici
                    %             costruite con i soli autovettori di cui sono sicuro
                        if ii~=Ncx && (~isfield(f,'eigA'))
                            if cond(real(RA_int_x(:,:,ii,jj))) > 1e4
                                dove = abs(DA_int_x(:,:,ii,jj)) < 1e-10; % Questi sono gli autovalori nulli
                                RA_int_x(:,dove,ii,jj) = 0;  % Uccido completamente gli autovettori relativi
                                RA_int_x(:,dove,ii,jj) = null(RA_int_x(:,:,ii,jj));  % li rimpiazzo con il null space
                            end
                        end
                        if jj~=Ncy && (~isfield(f,'eigB'))
                            if cond(real(RB_int_y(:,:,ii,jj))) > 1e4
                                dove = abs(DB_int_y(:,:,ii,jj)) < 1e-10; % Questi sono gli autovalori nulli
                                RB_int_y(:,dove,ii,jj) = 0;  % Uccido completamente gli autovettori relativi
                                RB_int_y(:,dove,ii,jj) = null(RB_int_y(:,:,ii,jj));  % li rimpiazzo con il null space
                            end
                        end
                    case 5
                        % Schur decomposition
                        if ii~=Ncx && (~isfield(f,'eigA'))
                            if cond(real(RA_int_x(:,:,ii,jj))) > 1e4
                                [RA_int_x(:,:,ii,jj),DA_int_x(:,1,ii,jj)] = my_decomp1(A_int_x(:,:,ii,jj));
                            end
                        end
                        if jj~=Ncy && (~isfield(f,'eigB'))
                            if cond(real(RB_int_y(:,:,ii,jj))) > 1e4
                                [RB_int_y(:,:,ii,jj),DB_int_y(:,1,ii,jj)] = my_decomp1(B_int_y(:,:,ii,jj));
                            end
                        end
                end
            end
        end
        % The eigenvectors are slightly complex. We take the real part. This is
        % not a problem because the eigenvectors are going to be approximated,
        % but the scheme is still OK
        RA_int_x = real(RA_int_x);
        RB_int_y = real(RB_int_y);
    elseif ((~isfield(f,'eigA')) || (~isfield(f,'eigB')) ) && isfield(f,'fix_eigs')
        [RA_int_x,DA_int_x,RB_int_y,DB_int_y] = f.fix_eigs(Ncx,Ncy,u_interf_x,RA_int_x,DA_int_x,u_interf_y,RB_int_y,DB_int_y,f);
    end
    %% First-order along x
    F = f.x(u_n,dV.xy);
    DF = F(:,i+1,:) - F(:,i,:) ;
    if isfield(f,'phix')
        PHIx = f.phix(u_n,dV.xy,t) .* reshape(dV.x,1,Ncx,1); % Questo è il source in ogni cella, concentrato in delta di Dirac
        PHIx = app_BCs(PHIx,BCs.BCsPHIx,t,true,2);
        DF = DF - ( PHIx(restr,i+1,:) + PHIx(restr,i,:) )/2;   % Definisco il source all'interfaccia come la media dei source nelle celle       
    end
    %%% Faccio la proiezione sulle f-waves
    DF = reshape(DF, sz(1), 1, Ncx-1,Ncy);
    betax = pagemldivide(RA_int_x,DF);
    Zx = RA_int_x .* pagetranspose(betax);

    %%% Costruisco le variazioni (o fluttuazioni)
    ApDQ = pagemtimes(Zx, double(DA_int_x>0) );
    AmDQ = pagemtimes(Zx, double(DA_int_x<0) );


    %% First-order along y
    G = f.y(u_n,dV.xy);
    DG = G(:,:,j+1) - G(:,:,j) ;
    if isfield(f,'phiy')
        PHIy = f.phiy(u_n,dV.xy,t) .* reshape(dV.y,1,1,Ncy); % Questo è il source in ogni cella, concentrato in delta di Dirac
        PHIy = app_BCs(PHIy,BCs.BCsPHIy,t,false,2);
        DG = DG - ( PHIy(restr,:,j+1) + PHIy(restr,:,j) )/2;   % Definisco il source all'interfaccia come la media dei source nelle celle       
    end
    %%% Faccio la proiezione sulle f-waves
    DG = reshape(DG, sz(1), 1, Ncx,Ncy-1);
    betay = pagemldivide(RB_int_y,DG);
    Zy = RB_int_y.* pagetranspose(betay);

    %%% Costruisco le variazioni
    BpDQ = pagemtimes(Zy, double(DB_int_y>0) );
    BmDQ = pagemtimes(Zy, double(DB_int_y<0) );

    if any(isnan(Zx(:))) || any(isnan(Zy(:)))
        error('The eigenvector decomposition failed')
    end
    
    %% (1a )Calcolo il correttore del flusso del secondo ordine lungo la direzione x e lo limito
    % Uso il metodo dopo Eq. (22), a fondo pagina, del paper 1997 Leveque - Wave propagation algorithms for hyperbolic systems
    i = 2:Ncx-2;
    thetax = zeros(1,sz(1),Ncx-1,Ncy);
    for p = 1:sz(1)   %  Z = sz(1) x sz(1) x Nc-1
    %     i_uw = i(DA_int_x(p,1,i,:) > 0); %  Mi confronto con le celle upwind
    %     i_dw = i(DA_int_x(p,1,i,:) < 0);
    %     i_uw = DA_int_x(p,1,i,:) > 0; %  Mi confronto con le celle upwind
    %     i_dw = DA_int_x(p,1,i,:) < 0;
    %     if any(i_uw)
    %         thetax(1,p,i_uw,:) =   pagemtimes(Zx(:,p,i_uw,:),'transpose', Zx(:,p,i_uw-1,:),'none') ...
    %                             ./ pagemtimes(Zx(:,p,i_uw,:),'transpose', Zx(:,p,i_uw,:)  ,'none');
    %     end
    %     if any(i_dw)
    %         thetax(1,p,i_dw,:) =   pagemtimes(Zx(:,p,i_dw,:),'transpose', Zx(:,p,i_dw+1,:),'none') ...
    %                             ./ pagemtimes(Zx(:,p,i_dw,:),'transpose', Zx(:,p,i_dw,:)  ,'none');
    %     end
        thetax(1,p,i,:) = (pagemtimes(Zx(:,p,i,:),'transpose', Zx(:,p,i-1,:),'none') .* (DA_int_x(p,1,i,:) > 0) + ...   Mi confronto con le celle upwind
                           pagemtimes(Zx(:,p,i,:),'transpose', Zx(:,p,i+1,:),'none') .* (DA_int_x(p,1,i,:) < 0)   ...
                          )./ pagemtimes(Zx(:,p,i,:),'transpose', Zx(:,p,i,:) ,'none');
    end
    % pphi = max(0,min(1,thetax));                     % minmod
    % pphi = max(0,max(min(1,2*thetax),min(2,thetax))); % Superbee
    pphi = max(0, min(min((1+thetax)/2, 2), 2*thetax));  % monotonized centered
    Zx_lim = Zx .* pphi;
    % Zx_lim = Zx ;      % No limiter
    
    %%% Calcolo la correzione dei flussi F_{i+-1/2}
    coeff = 0.5*sign(DA_int_x).*( 1 - dt * abs(DA_int_x)./reshape(dxx,1,1,Ncx-1,1) ); %  sz(1) x 1 x Ncx-1 x 1
    Flux_corrector_x = pagemtimes(  Zx_lim, coeff );
    F_i = squeeze(Flux_corrector_x); % Questo è il flusso principale limitato

    %% (1b) Calcolo il correttore del flusso del secondo ordine lungo la direzione y e lo limito
    j = 2:Ncy-2;
    thetay = zeros(1,sz(1),Ncx,Ncy-1);
    for p = 1:sz(1)   %  Z = sz(1) x sz(1) x Nc-1
        % j_uw = j(DB_int_y(p,:,j) > 0); %  Mi confronto con le celle upwind
        % j_dw = j(DB_int_y(p,:,j) < 0);
        % if any(j_uw)
        %     thetay(1,p,:,j_uw) =   pagemtimes(Zy(:,p,:,j_uw),'transpose', Zy(:,p,:,j_uw-1),'none') ...
        %                         ./ pagemtimes(Zy(:,p,:,j_uw),'transpose', Zy(:,p,:,j_uw)  ,'none');
        % end
        % if any(i_dw)
        %     thetay(1,p,:,i_dw) =   pagemtimes(Zy(:,p,:,j_dw),'transpose', Zy(:,p,:,j_dw+1),'none') ...
        %                         ./ pagemtimes(Zy(:,p,:,j_dw),'transpose', Zy(:,p,:,j_dw)  ,'none');
        % end
        thetay(1,p,:,j) = (pagemtimes(Zy(:,p,:,j),'transpose', Zy(:,p,:,j-1),'none') .* (DB_int_y(p,1,:,j) > 0) + ...   Mi confronto con le celle upwind
                           pagemtimes(Zy(:,p,:,j),'transpose', Zy(:,p,:,j+1),'none') .* (DB_int_y(p,1,:,j) < 0)   ...
                          )./ pagemtimes(Zy(:,p,:,j),'transpose', Zy(:,p,:,j) ,'none');
    end
    % pphi = max(0,min(1,thetax));                     % minmod
    % pphi = max(0,max(min(1,2*thetax),min(2,thetax))); % Superbee
    pphi = max(0, min(min((1+thetay)/2, 2), 2*thetay));  % monotonized centered
    Zy_lim = Zy .* pphi;
    % Zy_lim = Zy ;      % No limiter
    
    %%% Calcolo la correzione dei flussi F_{i+-1/2}
    coeff = 0.5*sign(DB_int_y).*( 1 - dt * abs(DB_int_y)./reshape(dyy,1,1,1,Ncy-1) ); %  sz(1) x 1 x 1 x Ncy-1
    Flux_corrector_y = pagemtimes(  Zy_lim, coeff );
    G_i = squeeze(Flux_corrector_y); % Questo è il correttore del flusso al secondo ordine limitato

%% Costruisco le transverse variations
% Splitto le fluttuazioni ApDQ e AmDQ in verticale
    % la dimensione dei seguenti è:
    %  n_stati  x   n_stati          x  Nx  x  Ny
    %   vett    x   f-split along y  x  Nx  x  Ny
    betax = pagemldivide(RB_int_x,ApDQ - 2*Flux_corrector_x); % Applico la correzione secondo la sez. 3.4 di 1997 Leveque
    Zxx = RB_int_x .* pagetranspose(betax);
    BpApDQ = squeeze( pagemtimes(Zxx, double(DB_int_x>0))  );
    BmApDQ = squeeze( pagemtimes(Zxx, double(DB_int_x<0))  );

    betax = pagemldivide(RB_int_x,AmDQ + 2*Flux_corrector_x);
    Zxx = RB_int_x .* pagetranspose(betax);
    BpAmDQ = squeeze( pagemtimes(Zxx, double(DB_int_x>0))  );
    BmAmDQ = squeeze( pagemtimes(Zxx, double(DB_int_x<0))  );
    
    ApDQ = squeeze(ApDQ);
    AmDQ = squeeze(AmDQ);

% Splitto le fluttuazioni BpDQ e BmDQ in orizzontale
    betay = pagemldivide(RA_int_y,BpDQ + 2*Flux_corrector_y);
    Zyy= RA_int_y .* pagetranspose(betay);
    ApBpDQ = squeeze( pagemtimes(Zyy, double(DA_int_y>0))  );
    AmBpDQ = squeeze( pagemtimes(Zyy, double(DA_int_y<0))  );

    betay = pagemldivide(RA_int_y,BmDQ - 2*Flux_corrector_y);
    Zyy = RA_int_y .* pagetranspose(betay);
    ApBmDQ = squeeze( pagemtimes(Zyy, double(DA_int_y>0))  );
    AmBmDQ = squeeze( pagemtimes(Zyy, double(DA_int_y<0))  );
    
    BpDQ = squeeze(BpDQ);
    BmDQ = squeeze(BmDQ);

    %% (2) Correggo i flussi usando le fluttuazioni trasversali
% Considero l'interfaccia verticale in (i+1/2, j) e distribuisco i suoi contributi sulle interf. orizzontali

    %%% Influenza che le interfacce verticali hanno su quelle orizzontali al di sopra
    i = 1:Ncx-1; % Tutte le interf. verticali lungo la direzione x
    j = 1:Ncy-1; % Tutte le interf. verticali lungo la direzione y tranne l'ultima perché non ha un'interf. orizz. sopra
    
    G_i(:,i+1,j)   = G_i(:,i+1,j) - 0.5*dt* BpApDQ(:,i,j) ./dxx; % pareti sopra a dx   (quando valuto g, sopra = j, sotto = j-1)
    G_i(:,i,j)     = G_i(:,i,j)   - 0.5*dt* BpAmDQ(:,i,j) ./dxx; % pareti sopra a sx   (                  dx = i+1, sx = i     )
    
    %%% Influenza che le interfacce verticali hanno su quelle orizzontali al di sopra
    j = 2:Ncy; % Tutte le interf. verticali lungo la direzione y tranne la prima perché non ha un'interf. orizz. sotto

    G_i(:,i+1,j-1) = G_i(:,i+1,j-1) - 0.5*dt* BmApDQ(:,i,j) ./dxx;  % pareti sotto a dx
    G_i(:,i,j-1)   = G_i(:,i,j-1)   - 0.5*dt* BmAmDQ(:,i,j) ./dxx;  % pareti sotto a sx

% Considero l'interfaccia orizzontale in (i, j+1/2) e distribuisco i suoi contributi sulle interf. verticali
    %%% Influenza che le interfacce oriz. hanno su quelle vert. a dx
    i = 1:Ncx-1; % Tutte le interf. oriz. lungo la direzione x tranne l'ultima perché non ha un'interf. vert. a dx
    j = 1:Ncy-1; % Tutte le interf. oriz. lungo la direzione y
    
    F_i(:,i,j+1)   = F_i(:,i,j+1) - 0.5*dt* ApBpDQ(:,i,j) ./dyy; % pareti sopra a dx   (quando valuto F,   dx = i  ,    sx = i-1 )
    F_i(:,i,j)     = F_i(:,i,j)   - 0.5*dt* ApBmDQ(:,i,j) ./dyy; % pareti sotto a dx   (                sopra = j+1, sotto = j   )
    
    %%% Influenza che le interfacce oriz. hanno su quelle vert. a sx
    i = 2:Ncx; % Tutte le interf. oriz. lungo la direzione x tranne la prima perché non ha un'interf. vert. a sx

    F_i(:,i-1,j+1) = F_i(:,i-1,j+1) - 0.5*dt* AmBpDQ(:,i,j) ./dyy;  % pareti sopra a sx
    F_i(:,i-1,j)   = F_i(:,i-1,j)   - 0.5*dt* AmBmDQ(:,i,j) ./dyy;  % pareti sotto a sx


    %% Faccio la propagazione Eq. (3.26)
    i = 3:Ncx-2;
    j = 3:Ncy-2;
    % Nota che ho aggiunto già un po' di ,: come se fossimo già in 2D
    u_n1(restr,i,j) = u_n(restr,i,j) - dt*(   AmDQ(:,i,j)   ./dxx(1,i,1)  ...
                                            + ApDQ(:,i-1,j) ./dxx(1,i-1,1) ) ...
                                     - dt*(   BmDQ(:,i,j)   ./dyy(1,1,j)  ...
                                            + BpDQ(:,i,j-1) ./dyy(1,1,j-1) ) ...
                                     - dt*(    F_i(:,i,j)   ./dxx(1,i,1)   ...
                                            -  F_i(:,i-1,j) ./dxx(1,i-1,1) ) ...
                                     - dt*(    G_i(:,i,j)   ./dyy(1,1,j)   ...
                                            -  G_i(:,i,j-1) ./dyy(1,1,j-1) );
    %% Metto le correzioni del secondo ordine di phi
    if isfield(f,'phix')
        dphi_dq = f.dphix_dq(u_n,dV.xy,t);
        dxx = reshape(dxx,1,1,Ncx-1,1);
        u_n1(restr,i,:) = u_n1(restr,i,:) - dt^2/4* squeeze(  pagemtimes(dphi_dq(:,:,i,:), sum( ... terms to add that depends on dphi/dW. See Eq. (4.45) of Bale 2002
                                                     Zx_lim(:,:,i,:)   ./dxx(1,1,i,1) ...
                                                   + Zx_lim(:,:,i-1,:) ./dxx(1,1,i-1,1)...
                                               ,2)  ))       ; 
    end
    if isfield(f,'phiy')
        dphi_dq = f.dphiy_dq(u_n,dV.xy,t);
        dyy = reshape(dyy,1,1,1,Ncy-1);
        u_n1(restr,:,j) = u_n1(restr,:,j) - dt^2/4* squeeze(  pagemtimes(dphi_dq(:,:,:,j), sum( ... terms to add that depends on dphi/dW. See Eq. (4.45) of Bale 2002
                                                     Zy_lim(:,:,:,j)   ./dyy(1,1,1,j) ...
                                                   + Zy_lim(:,:,:,j-1) ./dyy(1,1,1,j-1)...
                                               ,2)  ))       ; 
    end


    if ~isreal(u_n1)
        warning('Cerco di sistemare con la parte reale !!!')
        u_n1 = real(u_n1);
    end

end

end

% function [Vi,lambda] = my_decomp1(A)
%     % A è la matrice in una sola cella
%     % applico questa procedura solo alle celle critiche perché è più stabile
%     % ma più costosa
%     tol = 1e-10;
% 
%     [Q,T] = schur(A,'real');
% 
%     lambda = ordeig(T);
%     [~,idx] = sort(abs(lambda),'ascend');
%     lambda = lambda(idx);
% 
%     % Count repetitions of elements in lambda
%     [lmb_uniq, ~, ic] = uniquetol(lambda,tol);
%     n_gruppi = length(lmb_uniq);
%     counts = groupcounts(ic);
%     counts_cum = [0; cumsum(counts)];
% 
%     Vi = zeros(size(A));
% 
%     for gg = 1:n_gruppi
%         % Function identifying eigenvalues belonging to group 1
%         selez = ismembertol(lambda,lmb_uniq(gg),tol);
% 
%         [Qtmp,~] = ordschur(Q,T,selez);
%         Vi(:,counts_cum(gg)+1:counts_cum(gg+1)) = Qtmp(:,1:counts(gg));
%     end
% 
%     % % Right invariant subspaces
%     % Vr1 = V1;
%     % Vr2 = V2;
%     % Vr3 = V3;
% 
%     % % Left invariant subspaces
%     % Vl1 = inv(Vr1'*Vr1)*Vr1';
%     % Vl2 = inv(Vr2'*Vr2)*Vr2';
%     % Vl3 = inv(Vr3'*Vr3)*Vr3';
% end

% function my_decomp2(A)
%     Proj = zeros(sz(1),sz(1),sz(1),sz(2),sz(3));
% 
%     D = pageeig(A);
% 
%     % for gg = 1:
%     % for k = 1:sz(1)
%     %     if ~ismembertol(lambda(k),lambda1,1e-10)
%     %         P1 = P1 * (A-lambda(k)*eye(4)) ...
%     %              /(lambda(1)-lambda(k));
%     %     end
%     % end
% end
