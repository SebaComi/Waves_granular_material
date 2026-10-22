function [sigma_c,eps_c] = fun_cont_hertz_homo(rho,F11,F12,F21,F22,DATA)
    global rho0 a1 a2 a3
    
    % [~,Y_p0,rho0,E_p,nu_p] = data_gas;
    kH = DATA.kH;
    n = DATA.n;
    rho0 = DATA.rho0;

    % Questo è lo sforzo dovuto ai contatti di 3 direzioni: a1, a2 e a3
    sigma_c = [];
    sigma_c.xx = @(W) sigma_xx_fun(a1,W) + sigma_xx_fun(a2,W) + sigma_xx_fun(a3,W);
    sigma_c.xy = @(W) sigma_xy_fun(a1,W) + sigma_xy_fun(a2,W) + sigma_xy_fun(a3,W);
    sigma_c.yx = @(W) sigma_xy_fun(a1,W) + sigma_xy_fun(a2,W) + sigma_xy_fun(a3,W);
    sigma_c.yy = @(W) sigma_yy_fun(a1,W) + sigma_yy_fun(a2,W) + sigma_yy_fun(a3,W);

    % Energia dei 3 contatti (specifica sul volume)
    eps_c = @(W) e_c_fun(a1,W) + e_c_fun(a2,W) + e_c_fun(a3,W);

    % function out = eps_lmb_fun(W)
    %     %%% Da verificare se questa funzione costa meno
    %     lmb_val = lmb(W);
    %     out = 1/3*E_p/(1-nu_p^2)/Y_p0/rho0 *(  (lmb_val-1).*(lmb_val>1) ./(lmb_val.*(lmb_val>1)+1e-12)  ).^nn   ;
    %     % out = out + 1/200;  % Sto sommando all'energia (lmb-1)/1000
    % end
    % function out = eps_lmblmb_fun(W)
    %     lmb_val = lmb(W);
    %     out = nn/3*E_p/(1-nu_p^2)/Y_p0/rho0 *(lmb_val.*(lmb_val>1)+1e-12).^(-nn-1).*((lmb_val-1).*(lmb_val>1)).^(nn-1)  ;
    % end


    function out = sigma_xx_fun(a,W)
    % - rho kH (a1.*F11+a2.*F12).^2.*((a1.*F11+a2.*F12).^2+(a1.*F21+a2.*F22).^2).^(-1/2)  .*  (aa-((a1.*F11+a2.*F12).^2+(a1.*F21+a2.*F22).^2).^(1/2)).^n
        corrector = rho0./rho(W) ./( F11(W).*F22(W)-F12(W).*F21(W) );
        b1 = ( a(1).*F11(W)+a(2).*F12(W)  ) .*corrector; % I estimate the determinant from the density
        b2 = ( a(1).*F21(W)+a(2).*F22(W)  ) .*corrector;
        bb = sqrt( b1.^2 +  b2.^2);
        xi = sqrt(a(1)^2 + a(2)^2) - bb;
        out = -kH*rho(W).*b1.^2 ./bb .*  ( xi.*(xi>0) ).^n;
    end

    function out = sigma_xy_fun(a,W)
    % - rho kH (a1.*F11+a2.*F12).*(a1.*F21+a2.*F22).*((a1.*F11+a2.*F12).^2+(a1.*F21+a2.*F22).^2).^(-1/2).*(aa -((a1.*F11+a2.*F12).^2+(a1.*F21+a2.*F22).^2).^(1/2)).^n
        corrector = rho0./rho(W) ./( F11(W).*F22(W)-F12(W).*F21(W) );
        b1 = ( a(1).*F11(W)+a(2).*F12(W)  ) .*corrector; % I estimate the determinant from the density
        b2 = ( a(1).*F21(W)+a(2).*F22(W)  ) .*corrector;
        bb = sqrt( b1.^2 +  b2.^2);
        xi = sqrt(a(1)^2+a(2)^2) - bb;
        out = -kH*rho(W).*b1.*b2./bb .* ( xi.*(xi>0) ).^n;
    end

    function out = sigma_yy_fun(a,W)
    % - rho kH (a1.*F21+a2.*F22).^2   .*  (  (a1.*F11+a2.*F12).^2 + (a1.*F21+a2.*F22).^2).^(-1/2).*(aa -( (a1.*F11+a2.*F12).^2+(a1.*F21+a2.*F22).^2 ).^(1/2)).^n
        corrector = rho0./rho(W) ./( F11(W).*F22(W)-F12(W).*F21(W) );
        b1 = ( a(1).*F11(W)+a(2).*F12(W)  ) .*corrector; % I estimate the determinant from the density
        b2 = ( a(1).*F21(W)+a(2).*F22(W)  ) .*corrector;
        bb = sqrt( b1.^2 +  b2.^2);
        xi = sqrt(a(1)^2+a(2)^2) - bb;
        out = -kH*rho(W).*b2.^2./bb .* (xi.*(xi>0)).^n;
    end

    function out = e_c_fun(a,W)
        corrector = rho0./rho(W) ./( F11(W).*F22(W)-F12(W).*F21(W) );
        b1 = ( a(1).*F11(W)+a(2).*F12(W)  ) .*corrector; % I estimate the determinant from the density
        b2 = ( a(1).*F21(W)+a(2).*F22(W)  ) .*corrector;
        bb = sqrt( b1.^2 +  b2.^2);
        xi = sqrt(a(1)^2+a(2)^2) - bb;
        out = kH/(n+1) *(xi.*(xi>0)).^(n+1);
    end
end
