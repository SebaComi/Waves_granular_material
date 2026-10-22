function vel_out = c_fun_Comi_3_contacts(W,rho_fun,u_fun,v_fun,F11_fun,F12_fun,F21_fun,F22_fun, DATA, a1_vec,a2_vec,a3_vec, quale)

% Le velocità lungo la direzione x sono [0, 0, u, u, u,  u +/- sqrt(c1 +/- sqrt(c2)) ]

[u1,u2,kH,n,~,~,~,~,~,~,b1,b2,~,~,d1,d2,~,~,f1,f2,xi1,xi2,xi3,] = estrai(W,a1_vec,a2_vec,a3_vec,rho_fun,u_fun,v_fun,F11_fun,F12_fun,F21_fun,F22_fun,DATA);
                                                                        
% rho = rho_fun(W);

% corrector = rho0./rho ./( F11.*F22-F12.*F21 );
% b1 = ( a1(1).*F11+a1(2).*F12  ) .*corrector; % I estimate the determinant from the density
% b2 = ( a1(1).*F21+a1(2).*F22  ) .*corrector;
% d1 = ( a2(1).*F11+a2(2).*F12  ) .*corrector; % I estimate the determinant from the density
% d2 = ( a2(1).*F21+a2(2).*F22  ) .*corrector;
% f1 = ( a3(1).*F11+a3(2).*F12  ) .*corrector; % I estimate the determinant from the density
% f2 = ( a3(1).*F21+a3(2).*F22  ) .*corrector;
% 
% b0 = sqrt(b1.^2+ b2.^2);
% xi1 = norm(a1) - b0;
% xi1 = xi1.*(xi1>0);
% d0 = sqrt(d1.^2+ d2.^2);
% xi2 = norm(a2) - d0;
% xi2 = xi2.*(xi2>0);
% f0 = sqrt(f1.^2+ f2.^2);
% xi3 = norm(a3) - f0;
% xi3 = xi3.*(xi3>0);
b0 = sqrt(b1.^2+b2.^2);
d0 = sqrt(d1.^2+d2.^2);
f0 = sqrt(f1.^2+f2.^2);

if strcmp(quale,'x') || strcmp(quale,'max')  || strcmp(quale,'min')  || strcmp(quale,'max_x') || strcmp(quale,'imag')
    A1 = n* b1.^2      .*( xi1 ).^(n-1); % Qui scelgo b_1 perché sto guardando lungo l'asse x
    B1 =    b1.^2./b0 .*( xi1 ).^n;
    A2 = n* d1.^2      .*( xi2 ).^(n-1); % Qui scelgo b_1 perché sto guardando lungo l'asse x
    B2 =    d1.^2./d0 .*( xi2 ).^n;
    A3 = n* f1.^2      .*( xi3 ).^(n-1); % Qui scelgo b_1 perché sto guardando lungo l'asse x
    B3 =    f1.^2./f0 .*( xi3 ).^n;
    
    c1x = A1-B1 + A2-B2 + A3-B3;
    c2x = (A1+B1 + A2+B2 + A3+B3).^2 ...
          - 4*(A1+B1).*(A2+B2) .*(b1.*d2-b2.*d1).^2 ./b0.^2./d0.^2 ...
          - 4*(A2+B2).*(A3+B3) .*(d1.*f2-d2.*f1).^2 ./d0.^2./f0.^2 ...
          - 4*(A3+B3).*(A1+B1) .*(f1.*b2-f2.*b1).^2 ./f0.^2./b0.^2;
    c1x = c1x*kH;
    c2x = c2x*kH^2;
end
if strcmp(quale,'y') || strcmp(quale,'max')  || strcmp(quale,'min') || strcmp(quale,'max_y') || strcmp(quale,'imag')
    A1 = n* b2.^2     .*( xi1 ).^(n-1); % Qui scelgo b_2 perché sto guardando lungo l'asse y
    B1 =    b2.^2./b0 .*( xi1 ).^n;
    A2 = n* d2.^2     .*( xi2 ).^(n-1); % Qui scelgo b_2 perché sto guardando lungo l'asse y
    B2 =    d2.^2./d0 .*( xi2 ).^n;
    A3 = n* f2.^2     .*( xi3 ).^(n-1); % Qui scelgo b_2 perché sto guardando lungo l'asse y
    B3 =    f2.^2./f0 .*( xi3 ).^n;

    c1y = A1-B1 + A2-B2 + A3-B3;
    c2y = (A1+B1 + A2+B2 + A3+B3).^2 ...
          - 4*(A1+B1).*(A2+B2) .*(b1.*d2-b2.*d1).^2 ./b0.^2./d0.^2 ...
          - 4*(A2+B2).*(A3+B3) .*(d1.*f2-d2.*f1).^2 ./d0.^2./f0.^2 ...
          - 4*(A3+B3).*(A1+B1) .*(f1.*b2-f2.*b1).^2 ./f0.^2./b0.^2;
    c1y = c1y*kH;
    c2y = c2y*kH^2;
end

switch quale
    case 'x'
        c2xs = sqrt(c2x);
        vel_out = [u1 + sqrt( c1x + c2xs)/sqrt(2); 
                   u1 - sqrt( c1x + c2xs)/sqrt(2); 
                   u1 + sqrt( c1x - c2xs)/sqrt(2); 
                   u1 - sqrt( c1x - c2xs)/sqrt(2) ];
    case 'y'
        c2ys = sqrt(c2y);
        vel_out = [u2 + sqrt( c1y + c2ys)/sqrt(2); 
                   u2 - sqrt( c1y + c2ys)/sqrt(2); 
                   u2 + sqrt( c1y - c2ys)/sqrt(2);
                   u2 - sqrt( c1y - c2ys)/sqrt(2) ];
    case 'max'
        vel_out = max([abs(u1) + sqrt( c1x + sqrt(c2x) )/sqrt(2);  ...
                       abs(u2) + sqrt( c1y + sqrt(c2y) )/sqrt(2);  ...
                      ], [], 1);
    case 'max_x'
        vel_out = abs(u1) + sqrt( c1x + sqrt(c2x) )/sqrt(2);
    case 'max_y'
        vel_out = abs(u2) + sqrt( c1y + sqrt(c2y) )/sqrt(2);
    case 'imag'
        vel_out = max(abs(imag( sqrt( c1x - sqrt(c2x))/sqrt(2) )),  abs(imag( sqrt( c1y - sqrt(c2y))/sqrt(2) )) );
end


if ~all(isreal(vel_out(:)))
    % warning('Imaginary velocities')
    vel_out = real(vel_out);
end
end