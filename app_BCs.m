function u_n = app_BCs(u_n,BCs,t,x_sweep,n_ghost)
% Apply boundary conditions based on the specified BCs
if ~exist('x_sweep','var')
    x_sweep = true;
end
if ~exist('n_ghost','var')
    n_ghost = 1;
end

for gg = 1:n_ghost % The cell gg is the ghost cell counted from the boundary (1 is the first from the boundary)
    if x_sweep
        switch BCs.L.type % Left boundary
            case 'abs'
                u_n(:,gg,:) =         u_n(:,2*n_ghost-gg+1,:);
            case 'wall'
                u_n(:,gg,:) = BCs.L.v(u_n(:,2*n_ghost-gg+1,:),t);
            case 'periodic'
                u_n(:,gg,:) =         u_n(:,end-n_ghost+gg,:);
                BCs.R.type = 'periodic';
            case 'some' % This applies the BC only to some variables
                some = BCs.L.some;
                u_n(some,gg,:) = BCs.L.v(u_n(some,2*n_ghost-gg+1,:),t);                
        end
        switch BCs.R.type  % Right boundary
            case 'abs'
                u_n(:,end-gg+1,:) =         u_n(:,end-2*n_ghost+gg,:);
            case 'wall'
                u_n(:,end-gg+1,:) = BCs.R.v(u_n(:,end-2*n_ghost+gg,:),t);
            case 'periodic'
                u_n(:,end-gg+1,:) =         u_n(:,n_ghost-gg+1,:);
            case 'some' % This applies the BC only to some variables
                some = BCs.R.some;
                u_n(some,end-gg+1,:) = BCs.R.v(u_n(some,end-2*n_ghost+gg,:),t);
        end
    else
        switch BCs.B.type % Bottom boundary
            case 'abs'
                u_n(:,:,gg) =         u_n(:,:,2*n_ghost-gg+1);
            case 'wall'
                u_n(:,:,gg) = BCs.B.v(u_n(:,:,2*n_ghost-gg+1),t);
            case 'periodic'
                u_n(:,:,gg) =         u_n(:,:,end-n_ghost+gg);
                BCs.T.type = 'periodic';
            case 'some' % This applies the BC only to some variables
                some = BCs.B.some;
                u_n(some,:,gg) = BCs.B.v(u_n(some,:,2*n_ghost-gg+1),t);
        end
        switch BCs.T.type % Top boundary
            case 'abs'
                u_n(:,:,end-gg+1) =         u_n(:,:,end-2*n_ghost+gg);
            case 'wall'
                u_n(:,:,end-gg+1) = BCs.T.v(u_n(:,:,end-2*n_ghost+gg),t);
            case 'periodic'
                u_n(:,:,end-gg+1) =         u_n(:,:,n_ghost-gg+1);
            case 'some' % This applies the BC only to some variables
                some = BCs.T.some;
                u_n(some,:,end-gg+1) = BCs.T.v(u_n(some,:,end-2*n_ghost+gg),t);  
        end
    end
end