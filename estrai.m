function [u1,u2,kH,n,F11,F12,F21,F22,...
          a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,xi1,xi2,xi3,...
          AA,BB,CC,DD,EE,FF,AAA,BBB,CCC,DDD,EEE,FFF,AAA2,BBB2,CCC2,DDD2,EEE2,FFF2] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA)
    u1 = u(W);     u2 = v(W);
    kH = DATA.kH;   n = DATA.n;
    
    F11 = F11_fun(W);      F12 = F12_fun(W);
    F21 = F21_fun(W);      F22 = F22_fun(W);

    corrector = DATA.rho0./rho(W) ./( F11.*F22-F12.*F21 );

    a1 = a1_vec(1);
    a2 = a1_vec(2);
    b1 = ( a1*F11+a2*F12  ) .*corrector; % I estimate the determinant from the density
    b2 = ( a1*F21+a2*F22  ) .*corrector;
    b0 = sqrt( b1.^2 +  b2.^2);
    xi1 = sqrt(a1^2 + a2^2) - b0;
    detach1 = xi1<0;
    xi1(detach1) = 0;

    c1 = a2_vec(1);
    c2 = a2_vec(2);
    d1 = ( c1*F11+c2*F12  ) .*corrector; % I estimate the determinant from the density
    d2 = ( c1*F21+c2*F22  ) .*corrector;
    d0 = sqrt( d1.^2 +  d2.^2);
    xi2 = sqrt(c1^2 + c2^2) - d0;
    detach2 = xi2<0;
    xi2(detach2) = 0;

    e1 = a3_vec(1);
    e2 = a3_vec(2);
    f1 = ( e1*F11+e2*F12  ) .*corrector; % I estimate the determinant from the density
    f2 = ( e1*F21+e2*F22  ) .*corrector;
    f0 = sqrt( f1.^2 +  f2.^2);
    xi3 = sqrt(e1^2 + e2^2) - f0;
    detach3 = xi3<0;
    xi3(detach3) = 0;
    
    if any(detach1(:)) || any(detach2(:)) || any(detach3(:))
        warning('Detach detected!!')
    end

if nargout>23
    % Con la nuova convenzione, kH è sempre integrata in AAA e in AAA perché semplifica i calcoli in Matematica
    AA = kH*xi1.^(n-1)*n;
    BB = kH*xi1.^n./b0;
    CC = kH*xi2.^(n-1)*n;
    DD = kH*xi2.^n./d0;
    EE = kH*xi3.^(n-1)*n;
    FF = kH*xi3.^n./f0;

    AAA = AA.*b1.^2;
    BBB = BB.*b1.^2;
    CCC = CC.*d1.^2;
    DDD = DD.*d1.^2;
    EEE = EE.*f1.^2;
    FFF = FF.*f1.^2;
end
if nargout>35
    AA2 = kH*xi1.^(n-1)*n;
    BB2 = kH*xi1.^n./b0;
    CC2 = kH*xi2.^(n-1)*n;
    DD2 = kH*xi2.^n./d0;
    EE2 = kH*xi3.^(n-1)*n;
    FF2 = kH*xi3.^n./f0;

    AAA2 = AA2.*b2.^2;
    BBB2 = BB2.*b2.^2;
    CCC2 = CC2.*d2.^2;
    DDD2 = DD2.*d2.^2;
    EEE2 = EE2.*f2.^2;
    FFF2 = FF2.*f2.^2;
end
end