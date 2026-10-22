function [u_n1,dV] = solver(u0,dV,T,CFL,scheme,cosa_plot,salva,plots,filename,frame,auxiliary,globali_update)
global u x velocities DATA
if ~exist('frame','var')
    frame = 'Eul';
end

if exist('globali_update','var') && ~isempty(globali_update)
    globali = true;
    globali_update(u0,dV,0);
else
    globali = false;
end

if isfield(dV,'y')
    dim = 2; % Dimension of the problem
else
    dim = 1;
    % warning('Recentemente ho aggiungo l''input dV');
end
once = false;
figure(1)
cosa_plot(u0,0)
drawnow
save_GIF(salva,1,filename)
u_n1 = u0;
time_plot = T/50;  % I display about 50 times
t = 0; n = 0; n_plot = 1;
while t <= T

    % Salvo gli istanti temporali che mi servono
    if isfield(DATA,'time_to_save')
        if ~isempty(DATA.time_to_save)
            if t >= DATA.time_to_save(1)
                for fun = DATA.snaps
                    ff = figure(7);
                    fig_bella
                    fun{1}.plot(u_n1,t);
                    salvataggio(true,fun{1}.name(t))
                    clf
                end    
                close(ff)
                DATA.time_saved(end+1) = t;
                DATA.time_to_save(1) = [];
            end
        end
    end


    if t > 6e-3
        1;
    end
    n = n + 1;
    CFL_ok = false;
    % Compute the new time step
    if dim == 1
        dt = CFL * abs(min(dV.x(3:end-2) ./ max(velocities.x(u_n1(:,3:end-2)),[],1))); 
    else % dim == 2 The /sqrt(2) comes from the condition of travelling across a square cell
      % dt = CFL/sqrt(2) *  abs(min([dV.x,dV.y],[],'all') ./ max([velocities.x(u_n1(:,3:end-2,3:end-2));velocities.y(u_n1(:,3:end-2,3:end-2))],[],'all'));
        % dt = CFL         *  abs(min([dV.x,dV.y],[],'all') ./ max([velocities.x(u_n1(:,3:end-2,3:end-2));velocities.y(u_n1(:,3:end-2,3:end-2))],[],'all'));
        dt = CFL         *  abs( min([dV.x,dV.y],[],'all') ./ max(velocities.max(u_n1(:,3:end-2,3:end-2)),[],'all')   );
    end
    dt = min(dt,T/10); % In case dt is inf, typical of the initial case where no compaction is there
    while ~CFL_ok % I do the check if the change of velocities is acceptable at the end of the time step
        u_n_test = scheme(u_n1,t,dV,dt);
        % auxiliary(u_n_test);
        if dim == 1
            CFL_new = dt * abs(max(max(velocities.x(u_n_test(:,3:end-2)),[],1) ./ dV.x(3:end-2) ));
        else
            % old = max([velocities.x(u_n_test(:,3:end-2,3:end-2));velocities.y(u_n_test(:,3:end-2,3:end-2))],[],'all');
            new = abs(max(velocities.max(u_n_test),[],'all'));
            % CFL_new = dt*sqrt(2) * new ./ min([dV.x,dV.y],[],'all');
            CFL_new = dt * new ./ min([dV.x,dV.y],[],'all');
            % if all(isreal(new(:)))
            %     fprintf('New = %.2f     old = %.2f \n', new,old);
            % else
            %     fprintf('New = %.2f + %.2fi      old = %.2f + %.2fi  \n', real(new),imag(new),real(old), imag(old));
            % end
        end

        if CFL_new > 5 || dt <0
            warning('!!')
        end
        if CFL_new <= CFL*1.03 && CFL_new < 0.99 % Consento che il CFL aumenti del 3% rispetto al desiderato, ma che non sia >=0.99
            CFL_ok = true;
            t = t + dt;
            u_n1 = u_n_test;
        else
            dt = dt * sqrt(CFL/CFL_new); % riprovo con un CFL più basso (CFL/CFL_new < 1/1.03 = 0.97). Faccio sqrt per stare un po' a metà fra le due iterazioni
            % auxiliary(u_n1);
            if CFL_new > 5
                warning('!!')
            end
        end
    end
    
    if strcmp(frame,'Lag') % update of the mesh position and volumes 
        if dim == 2
            error('Impossible to update lagrangian coordinates in 2D')
        end
        uu = u(u_n1);
        x = x + dt*uu;
        dV(2:end-1) = dV(2:end-1) + dt* (uu(3:end)-uu(1:end-2))/2;
        dV(1) = dV(1) + dt* (uu(2)-uu(1))/2;
        dV(end) = dV(end) + dt* (uu(end)-uu(end-1))/2;
    end
    
    if (plots && t > n_plot*time_plot ) || t >= T %|| 1
        n_plot = n_plot + 1;
    % if (plots && (mod(n,plot_every)==0 || abs(mod(t, time_plot)) < dt/2 ) ) || t >= T
        cosa_plot(u_n1,t)
        %%% Questo se voglio vedere come cambiano i volumi
        % figure(2)
        % plot(dV)
        % figure(1)
        fprintf("Iter %d   -   t = %.2e s   -  dt = %.2e s  -  CFL = %.2f - Disegno\n",n, t, dt, CFL_new)

        drawnow
        save_GIF(salva,n+1,filename)
    elseif mod(n,20) == 0
        perc = t/T*100;
        fprintf("Iter %d   -   t = %.2e s   -  dt = %.2e s  -  CFL = %.2f - t_tot %.2f min - avanz %.1f %% \n",n, t, dt, CFL_new,toc/60,perc)
    end

    if globali
        globali_update(u_n1,dV,t);
    end

    % if t >= 0.05 && ~once   % Questo è un trucco per fare un'azione magica
    %     u_n1(2:3,:,:) = 0;
    %     once = true;
    % end
end
