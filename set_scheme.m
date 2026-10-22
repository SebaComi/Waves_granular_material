function scheme = set_scheme(scheme_type,dx,Nx,f,df,u_s,g,vel,BCs,restr)

if ~exist('BCs','var')
    BCs.L.type = 'abs'; % Left
    BCs.R.type = 'abs'; % Right
    BCs.T.type = 'abs'; % Top
    BCs.B.type = 'abs'; % Bottom
end
if isfield(f,'ncx')
    nc = true;
else
    nc = false;
end
if ~exist('restr','var')
    restr = [];
end
if iscell(scheme_type) % In this case I have to apply the various splitting
    % Alcuni input sono obsoleti (tipo Nx e df), quindi metto [] e via
    spl = 1;
    % spl = length(scheme_type);
    scheme = set_scheme(scheme_type{spl},dx,[],f{spl},[],u_s,g{spl},vel{spl},BCs,restr{spl});

    for spl = 2:length(scheme_type)
    % for spl = length(scheme_type)-1:-1:1
        scheme_tmp = set_scheme(scheme_type{spl},dx,[],f{spl},[],u_s,g{spl},vel{spl},BCs,restr{spl});
        scheme = @(u_n,t,dV,dt) scheme(scheme_tmp(u_n,t,dV,dt),t,dV,dt);
    end

else
    switch scheme_type
        case 'UW'
            scheme = @(u_n,t,dV,dt) upwind(u_n,dV,dt,Nx,f,BCs);
        case 'UWnc' % Upwind non conservative
            scheme = @(u_n,t,dV,dt) upwind_nc(u_n,dx,dt,Nx,dV,df);
        case 'DW'
            scheme = @(u_n,t,dV,dt) downwind(u_n,dV,dt,Nx,f,BCs);
        case 'DWnc' % Downwind non conservative
            scheme = @(u_n,t,dV,dt) downwind_nc(u_n,dx,dt,Nx,dV,df);
        case 'LF'
            scheme = @(u_n,t,dV,dt) Lax_Friedrichs(u_n,dV,dt,f,BCs,t);
        case 'LFnc'
            scheme = @(u_n,t,dV,dt) Lax_Friedrichs_nc(u_n,dx,dt,Nx,dV,df);
        case {'LW','LWnc'}
            scheme = @(u_n,t,dV,dt) Lax_Wendroff(u_n,dx,dt,Nx,dV,f,df);
        case 'RLW'
            scheme = @(u_n,t,dV,dt) Richtmyer_Lax_Wendroff(u_n,dx,dt,Nx,dV,f);
        case 'Go'
            scheme = @(u_n,t,dV,dt) Godunov(u_n,dx,dt,Nx,dV,f,df,u_s,BCs,t);
        case 'H_lin' % high order for linear scalar equations
            scheme = @(u_n,t,dV,dt) high_order_UW_LW(u_n,dx,dt,Nx,dV,df);
        case 'Ru'
            scheme = @(u_n,t,dV,dt) Rusanov(u_n,dV,dt,f,vel,BCs,t,true,nc,restr);
        case 'HLL'
            scheme = @(u_n,t,dV,dt) HLL(u_n,dV,dt,f,vel,BCs,t);
        case 'HLLC'
            error('This scheme is not set yet')
            scheme = @(u_n,t,dV,dt) HLLC(u_n,dV,dt,f,vel,BCs,t,stim);
        case 'MUSCL_H'
            scheme = @(u_n,t,dV,dt) MUSCL_H(u_n,dV,dt,f,vel,BCs,t,true,nc);
        case 'Bale'
            scheme = @(u_n,t,dV,dt) Bale(u_n,dV,dt,f,vel,BCs,t,restr);
        case 'LF2D'
            % dV in realtà non è il volume ma solo dx o dy (alias, è dV*dy o dV*dx, rispettivamente)
            %    visto che Nx e Ny possono essere diversi, salvo solo dV.x e dV.y
            scheme_x = @(u_n,t,dV,dt) Lax_Friedrichs(u_n,dV.x,dt,f.x,BCs,t, true);
            %%% Non posso semplicemente permutare le dimensioni 2 e 3 perché le funzioni interne (e.g. alpha, rho, u, ecc) non se ne accorgono...
            scheme_y = @(u_n,t,dV,dt) Lax_Friedrichs(u_n,dV.y,dt,f.y,BCs,t, false);
            scheme = @(u_n,t,dV,dt) scheme_y(scheme_x(u_n, t,dV,dt), t,dV,dt);
        case 'Ru2D'
            scheme_x = @(u_n,t,dV,dt) Rusanov(u_n,dV,dt,f,vel,BCs,t, true,nc,restr);
            scheme_y = @(u_n,t,dV,dt) Rusanov(u_n,dV,dt,f,vel,BCs,t, false,nc,restr);
            scheme = @(u_n,t,dV,dt) scheme_y(scheme_x(u_n, t,dV,dt), t,dV,dt);
        case 'HLL2D'
            scheme_x = @(u_n,t,dV,dt) HLL(u_n,dV.x,dt,f.x,vel,BCs,t, true);
            scheme_y = @(u_n,t,dV,dt) HLL(u_n,dV.y,dt,f.y,vel,BCs,t, false);
            scheme = @(u_n,t,dV,dt) scheme_y(scheme_x(u_n, t,dV,dt), t,dV,dt);
        case 'MUSCL_H2D'
            scheme_x = @(u_n,t,dV,dt) MUSCL_H(u_n,dV,dt,f,vel,BCs,t, true,nc,[],restr);
            scheme_y = @(u_n,t,dV,dt) MUSCL_H(u_n,dV,dt,f,vel,BCs,t, false,nc,[],restr);
            scheme = @(u_n,t,dV,dt) scheme_y(scheme_x(u_n, t,dV,dt), t,dV,dt);
        case 'MUSCL_H2D_curl'
            n_ghost = 2; n_vec = 2; nc = false;
            scheme_x    = @(u_n,t,dV,dt)           MUSCL_H(u_n,dV,dt,f,vel,BCs,t, true ,nc,2*n_vec);
            scheme_y    = @(u_n,t,dV,dt)           MUSCL_H(u_n,dV,dt,f,vel,BCs,t, false,nc,2*n_vec);
            scheme_curl = @(u_n,t,dV,dt)  curl_free_scheme(u_n,dV,dt,f,vel,BCs,t,n_ghost, n_vec);
            scheme      = @(u_n,t,dV,dt) scheme_curl(scheme_y(scheme_x(u_n, t,dV,dt), t,dV,dt), t,dV,dt);
        case 'MUSCL_H2D_curl_imp'
            n_ghost = 2; n_vec = 2;
            scheme_x    = @(u_n,t,dV,dt)              MUSCL_H(u_n,dV,dt,f,vel,BCs,t, true ,nc,n_vec);
            scheme_y    = @(u_n,t,dV,dt)              MUSCL_H(u_n,dV,dt,f,vel,BCs,t, false,nc,n_vec);
            scheme_curl = @(u_n,t,dV,dt) curl_free_scheme_imp(u_n,dV,dt,f,vel,BCs,t, n_ghost, n_vec);
            scheme      = @(u_n,t,dV,dt) scheme_curl(scheme_y(scheme_x(u_n, t,dV,dt), t,dV,dt), t,dV,dt);
        otherwise
            error('Scheme not defined')
    end

    if exist('g','var') % In this case I use the splitting method to handle the RHS g(t,x)
        
        
        if ~isempty(g)
            if isfield(g,'type')
                if ~isempty(restr)
                    error('I''ve not yet implemented the restriction for this type pf RHS')
                end
                switch g.type
                    case 'ana_exp' % Here I only have the part that I solve analytically
                        scheme = @(u_n,t,dV,dt)  scheme(analytic_exp(u_n,t,dt,g),t,dV,dt);
                    case 'exp_and_RK' % Here I have both: analytic and RK
                        scheme = @(u_n,t,dV,dt)  scheme(analytic_exp(u_n,t,dt,g),t,dV,dt);
                        scheme = @(u_n,t,dV,dt)  scheme(Runge_Kutta(u_n,t,dt,g.g),t,dV,dt);
                end
            else
                % I first propagate the solution due to g and then apply the scheme
                % if ~isfield(g,'solver')
                % Ho sempre fatto scheme(RK(u_n))
                    % scheme = @(u_n,t,dV,dt)  scheme(Runge_Kutta(u_n,t,dt,g),t,dV,dt);
                % Provo a fare RK(scheme(u_n))
                    % scheme = @(u_n,t,dV,dt)  Runge_Kutta(scheme(u_n,t,dV,dt),t,dt,g,restr);
                % else
                    scheme = @(u_n,t,dV,dt)  scheme(Runge_Kutta(u_n,t,dt,g,restr),t,dV,dt);  %% Questo sembra funzionare meglio!! (Cioè, se l'ultimo step è lo "scheme" la simulazione si ferma più vicina a un equilibrio statico
                % end
                % error('Questo è da controllare, soprattutto con i dV')
            end
        end
    end


end % end of the variable splitting


% if exist('auxiliary','var')
% % In questo caso chiamo la funzione ausiliaria per fare degli update
%     scheme = @(u_n,t,dV,dt) auxiliary_scheme(u_n,t,dV,dt, scheme);
% end
% 
% function out = auxiliary_scheme(u_n,t,dV,dt, scheme)
%     auxiliary(u_n); % Questo fa l'update della variabile globale alpha_p_old
%     out = scheme(u_n,t,dV,dt);
% end

end