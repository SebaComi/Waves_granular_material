function [EIGS,DDx] = eigs_x_3_contacts(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA)
% See my Matematica code: Hyperbolicity_Grains_F_no_beta_CLEAN_15.nb

sz = size(W);
EIGS = zeros(sz(1),sz(1),sz(2),sz(3));
DDx = zeros(sz(1),1,sz(2),sz(3));
% [u1,u2,kH,n,F11,F12,F21,F22,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,xi1,xi2,xi3,AA,BB,CC,DD,EE,FF,AAA,BBB,CCC,DDD,EEE,FFF] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA)
  [u1,u2, ~,~,F11,F12,F21,F22,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2,  ~,  ~,  ~,AA,BB,CC,DD,EE,FF,AAA,BBB,CCC,DDD,EEE,FFF] = estrai(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);

b0 = sqrt( b1.^2 +  b2.^2);
d0 = sqrt( d1.^2 +  d2.^2);
f0 = sqrt( f1.^2 +  f2.^2);
a1t = a1./b0;
a2t = a2./b0;
b1t = b1./b0;
b2t = b2./b0;
c1t = c1./d0;
c2t = c2./d0;
d1t = d1./d0;
d2t = d2./d0;
e1t = e1./f0;
e2t = e2./f0;
f1t = f1./f0;
f2t = f2./f0;

    %% The following is the output by Wolfram and some modifications
b0_2 = b0.^2;   b1_2 = b1.^2;   b2_2 = b2.^2;
d0_2 = d0.^2;   d1_2 = d1.^2;   d2_2 = d2.^2;
f0_2 = f0.^2;   f1_2 = f1.^2;   f2_2 = f2.^2;
b1_3 = b1.^3;   b2_3 = b2.^3;
d1_3 = d1.^3;   d2_3 = d2.^3;
f1_3 = f1.^3;   f2_3 = f2.^3;

cc1a =  AAA-BBB+CCC-DDD+EEE-FFF;
cc2a = -4*(AAA+BBB).*(b2t.*d1t-b1t.*d2t).^2.*(CCC+DDD)...
       -4*(AAA+BBB).*(b2t.*f1t-b1t.*f2t).^2.*(EEE+FFF)...
       -4*(CCC+DDD).*(d2t.*f1t-d1t.*f2t).^2.*(EEE+FFF)...
       +(AAA+BBB+CCC+DDD+EEE+FFF).^2;

cc2a_1_2 = cc2a.^(1/2);


DDx(1,1,:,:) = u1;
DDx(2,1,:,:) = u1;
DDx(3,1,:,:) = u1;
DDx(4,1,:,:) = u1 + sqrt(cc1a + cc2a_1_2)/sqrt(2);
DDx(5,1,:,:) = u1 - sqrt(cc1a + cc2a_1_2)/sqrt(2);
DDx(6,1,:,:) = u1 + sqrt(cc1a - cc2a_1_2)/sqrt(2);
DDx(7,1,:,:) = u1 - sqrt(cc1a - cc2a_1_2)/sqrt(2);

A = A_fun(W,a1_vec,a2_vec,a3_vec,rho,u,v,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
A25 = reshape( A(2,5,:,:), size(u1) );
A27 = reshape( A(2,7,:,:), size(u1) );
A35 = reshape( A(3,5,:,:), size(u1) );
A37 = reshape( A(3,7,:,:), size(u1) );

E = A27 .* A35 - A25 .* A37;

EIGS(1,1,:,:) = E;
EIGS(2,1,:,:) = u1 .*E;
EIGS(3,1,:,:) = u2 .*E;
% EIGS(4,1,:,:) = 0;
% EIGS(5,1,:,:) = (a2t.*b1t.*d0_2.*f0_2.*(b0_2.*b1t.*BB.*(-b1t.*CC.*d1t.^2.*f0_2-b2t.*(CC.*d1t.*d2t.*f0_2+d0_2.*EE.*f1t.*f2t)+b1t.*d0_2.*EE.*((-1)+f2t.^2))-AA.*(b1t.^2.*BB.*d0_2.*f0_2+b0_2.*b2t.*(-b2t.*CC.*d1t.^2.*f0_2+b1t.*(CC.*d1t.*d2t.*f0_2+d0_2.*EE.*f1t.*f2t)+b2t.*d0_2.*EE.*((-1)+f2t.^2))))-b0.^4.*(c2t.*CC.*d1t.^3.*DD.*f0.^4+c2t.*d0_2.*d1t.*EE.*f0_2.*(-CC.*d2t.^2.*f1t.^2+d1t.*d2t.*(CC+DD).*f1t.*f2t+d1t.^2.*(DD-DD.*f2t.^2))+d0.^4.*e2t.*EE.*f1t.^3.*FF+CC.*d0_2.*d1t.*e2t.*f0_2.*f1t.*(d2t.*f1t.*f2t.*(EE+FF)+d1t.*(FF-f2t.^2.*(EE+FF))))+AA.*b0_2.*b1t.*d0_2.*f0_2.*(-b2t.*(c2t.*d1t.^2.*d2t.*(CC+DD).*f0_2+d0_2.*e2t.*f1t.^2.*f2t.*(EE+FF))+b1t.*(c2t.*d1t.*(-DD+d2t.^2.*(CC+DD)).*f0_2+d0_2.*e2t.*f1t.*(-FF+f2t.^2.*(EE+FF))))).*(a2t.^2.*((-1)+b2t.^2).*(AA-BB).*BB.*d0.^4.*f0.^4-a2t.*b0_2.*d0_2.*f0_2.*(-b1t.*b2t.^2.*(AA+BB).*c2t.*CC.*d1t.*f0_2-AA.*b1t.*c2t.*d1t.*d2t.^2.*(CC+DD).*f0_2+2*b1t.*b2t.^2.*(AA+BB).*c2t.*d1t.*d2t.^2.*(CC+DD).*f0_2+2*b2t.^3.*(AA+BB).*c2t.*d2t.^3.*(CC+DD).*f0_2-b2t.*(2.*AA+BB).*c2t.*d2t.^3.*(CC+DD).*f0_2-b2t.^3.*(AA+BB).*c2t.*d2t.*(2.*CC+DD).*f0_2+b1t.*c2t.*d1t.*(BB.*(CC-2*DD)+AA.*DD).*f0_2+b2t.*c2t.*d2t.*(BB.*CC+AA.*(2.*CC+DD)).*f0_2-b1t.*b2t.^2.*(AA+BB).*d0_2.*e2t.*EE.*f1t-AA.*b1t.*d0_2.*e2t.*f1t.*f2t.^2.*(EE+FF)+2*b1t.*b2t.^2.*(AA+BB).*d0_2.*e2t.*f1t.*f2t.^2.*(EE+FF)+2*b2t.^3.*(AA+BB).*d0_2.*e2t.*f2t.^3.*(EE+FF)-b2t.*(2.*AA+BB).*d0_2.*e2t.*f2t.^3.*(EE+FF)-b2t.^3.*(AA+BB).*d0_2.*e2t.*f2t.*(2.*EE+FF)+b1t.*d0_2.*e2t.*f1t.*(BB.*(EE-2*FF)+AA.*FF)+b2t.*d0_2.*e2t.*f2t.*(BB.*EE+AA.*(2.*EE+FF)))+b0.^4.*(c2t.^2.*((-1)+d2t.^2).*(CC-DD).*DD.*f0.^4+d0.^4.*e2t.^2.*((-1)+f2t.^2).*(EE-FF).*FF-c2t.*d0_2.*e2t.*f0_2.*(-d1t.*d2t.^2.*(CC+DD).*EE.*f1t-CC.*d1t.*f1t.*f2t.^2.*(EE+FF)+2*d1t.*d2t.^2.*(CC+DD).*f1t.*f2t.^2.*(EE+FF)+2*d2t.^3.*(CC+DD).*f2t.^3.*(EE+FF)-d2t.*(2.*CC+DD).*f2t.^3.*(EE+FF)-d2t.^3.*(CC+DD).*f2t.*(2.*EE+FF)+d1t.*f1t.*(DD.*(EE-2*FF)+CC.*FF)+d2t.*f2t.*(DD.*EE+CC.*(2.*EE+FF))))).^(-1);
EIGS(5,1,:,:) = A37 .*(AAA + CCC + EEE) -  A27 .*(AA.*b1.*b2 + CC.*d1.*d2 + EE.*f1.*f2);
% EIGS(6,1,:,:) = 0;
% EIGS(7,1,:,:) = b0.^4.*d0.^4.*f0.^4.*((a2t.*b1t.^2.*b2t.*(AA+BB).*d0_2.*f0_2+b0_2.*(c2t.*d1t.^2.*d2t.*(CC+DD).*f0_2+d0_2.*e2t.*f1t.^2.*f2t.*(EE+FF))).*(-a2t.*b2t.*(-AA+b2t.^2.*(AA+BB)).*d0_2.*f0_2-b0_2.*(c2t.*d2t.*(-CC+d2t.^2.*(CC+DD)).*f0_2+d0_2.*e2t.*f2t.*(-EE+f2t.^2.*(EE+FF))))+(a2t.*b1t.*(BB-b2t.^2.*(AA+BB)).*d0_2.*f0_2+b0_2.*(c2t.*d1t.*(DD-d2t.^2.*(CC+DD)).*f0_2+d0_2.*e2t.*f1t.*(FF-f2t.^2.*(EE+FF)))).*(-a2t.*b1t.*(-AA+BB+b2t.^2.*(AA+BB)).*d0_2.*f0_2+b0_2.*(-c2t.*d1t.*(-CC+DD+d2t.^2.*(CC+DD)).*f0_2-d0_2.*e2t.*f1t.*(-EE+FF+f2t.^2.*(EE+FF))))).^(-1).*((-AA.*b0.^(-2).*b1t.^2-CC.*d0.^(-2).*d1t.^2-EE.*f0.^(-2).*f1t.^2).*(b0.^(-2).*(-a2t.*AA.*b2t+a2t.*b2t.^3.*(AA+BB))+c2t.*d0.^(-2).*d2t.*(-CC+d2t.^2.*(CC+DD))+e2t.*f0.^(-2).*f2t.*(-EE+f2t.^2.*(EE+FF)))+b0.^(-4).*d0.^(-4).*f0.^(-4).*(AA.*b1t.*b2t.*d0_2.*f0_2+b0_2.*(CC.*d1t.*d2t.*f0_2+d0_2.*EE.*f1t.*f2t)).*(a2t.*b1t.*(-AA+BB+b2t.^2.*(AA+BB)).*d0_2.*f0_2+b0_2.*(c2t.*d1t.*(-CC+DD+d2t.^2.*(CC+DD)).*f0_2+d0_2.*e2t.*f1t.*(-EE+FF+f2t.^2.*(EE+FF)))));
EIGS(7,1,:,:) = -A35.*(AAA + CCC + EEE) + A25 .*(AA.*b1.*b2 + CC.*d1.*d2 + EEE.*f1.*f2);

% EIGS(1,2,:,:) = 0;
% EIGS(2,2,:,:) = 0;
% EIGS(3,2,:,:) = 0;
EIGS(4,2,:,:) =  a1t.*b0.^(-2).*b1t.*(-BB+b2t.^2.*(AA+BB))   +c1t.*d0.^(-2).*d1t.*(-DD+d2t.^2.*(CC+DD))   +e1t.*f0.^(-2).*f1t.*(-FF+f2t.^2.*(EE+FF));
% EIGS(5,2,:,:) = 0;
EIGS(6,2,:,:) = -a1t.*b0.^(-2).*b1t.*(-AA+BB+b2t.^2.*(AA+BB))-c1t.*d0.^(-2).*d1t.*(-CC+DD+d2t.^2.*(CC+DD))-e1t.*f0.^(-2).*f1t.*(-EE+FF+f2t.^2.*(EE+FF));
% EIGS(7,2,:,:) = 0;

% EIGS(1,3,:,:) = 0;
% EIGS(2,3,:,:) = 0;
% EIGS(3,3,:,:) = 0;
EIGS(4,3,:,:) = -b0.^(-1).*d0.^(-1).*f0.^(-1).*( a2t.*b0_2.*b2t.*(-AA+b2t.^2.*(AA+BB))    +c2t.*d0_2.*d2t.*(-CC+d2t.^2.*(CC+DD))    +e2t.*f0_2.*f2t.*(-EE+f2t.^2.*(EE+FF))    );
EIGS(5,3,:,:) = -b0.^(-1).*d0.^(-1).*f0.^(-1).*( a1t.*b0_2.*b1t.*(-AA+BB+b2t.^2.*(AA+BB)) +c1t.*d0_2.*d1t.*(-CC+DD+d2t.^2.*(CC+DD)) +e1t.*f0_2.*f1t.*(-EE+FF+f2t.^2.*(EE+FF)) );
% EIGS(6,3,:,:) = 0;
% EIGS(7,3,:,:) = 0;

%%% Occhio che ci sono dei *2 in giro per colpa di cc1a e cc2a che non hanno 1/sqrt(2) inclusi
coso1 = cc1a+cc2a_1_2 +2.*BBB-2*b2t.^2.*(AAA+BBB) +2*DDD-2*d2t.^2.*(CCC+DDD) +2*FFF-2*f2t.^2.*(EEE+FFF);
coso2 = cc1a-cc2a_1_2 +2.*BBB-2*b2t.^2.*(AAA+BBB) +2*DDD-2*d2t.^2.*(CCC+DDD) +2*FFF-2*f2t.^2.*(EEE+FFF);
coso3 = b1t.*b2t.*(AAA+BBB)+d1t.*d2t.*(CCC+DDD)+f1t.*f2t.*(EEE+FFF);

EIGS(1,4,:,:) = coso1;
EIGS(2,4,:,:) = ( 2.^(-1/2).*(cc1a+cc2a_1_2).^(1/2)+u1 ) .*coso1;
EIGS(3,4,:,:) = 2.^(1/2).*(cc1a+cc2a_1_2).^(1/2).*coso3 + u2.*coso1;
% EIGS(4,4,:,:) = 0;
% EIGS(5,4,:,:) = 0;
EIGS(6,4,:,:) = F21.*coso1 - 2*F11.*coso3;
EIGS(7,4,:,:) = F22.*coso1 - 2*F12.*coso3;

EIGS(1,5,:,:) = coso1;
EIGS(2,5,:,:) = ( -2.^(-1/2).*(cc1a+cc2a_1_2).^(1/2)+u1 ) .*coso1;
EIGS(3,5,:,:) = -sqrt(2)*(cc1a+cc2a_1_2).^(1/2).*coso3 + u2.*coso1;
% EIGS(4,5,:,:) = 0;
% EIGS(5,5,:,:) = 0;
EIGS(6,5,:,:) = F21.*coso1 - 2*F11.*coso3;
EIGS(7,5,:,:) = F22.*coso1 - 2*F12.*coso3;

EIGS(1,6,:,:) = coso2;
EIGS(2,6,:,:) = coso2.* ( 2.^(-1/2).*(cc1a-cc2a_1_2).^(1/2)+u1 );
EIGS(3,6,:,:) = 2.^(1/2).*(cc1a-cc2a_1_2).^(1/2).*coso3 + u2.*coso2;
% EIGS(4,6,:,:) = 0;
% EIGS(5,6,:,:) = 0;
EIGS(6,6,:,:) = F21 .*coso2 - 2*F11.*coso3;
EIGS(7,6,:,:) = F22 .*coso2 - 2*F12.*coso3;

EIGS(1,7,:,:) = coso2;
EIGS(2,7,:,:) = coso2.* ( -2.^(-1/2).*(cc1a-cc2a_1_2).^(1/2)+u1 );
EIGS(3,7,:,:) = -sqrt(2)*(cc1a-cc2a_1_2).^(1/2).*coso3 + u2.*coso2;
% EIGS(4,7,:,:) = 0;
% EIGS(5,7,:,:) = 0;
EIGS(6,7,:,:) = F21.*coso2 - 2*F11.*coso3;
EIGS(7,7,:,:) = F22.*coso2 - 2*F12.*coso3;

end