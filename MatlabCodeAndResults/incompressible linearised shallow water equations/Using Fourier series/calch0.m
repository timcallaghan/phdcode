function y = calch0(w,Fr,Ro,arad)

% Zonal flow component that depends on w
A = w*Fr^2/2*(1/Ro+w);
% The cubic equation coefficients
a3 = 4*pi/3;
a2 = 8*pi*A/3 + 4*pi*arad;
a1 = 32*pi*A^2/15 + 16*pi*A*arad/3 + 4*pi*arad^2;
a0 = 64*pi*A^3/105 + 32*pi*A^2/15*arad + 8*pi*A*arad^2/3 - 4*pi/3 - 4*pi*arad*(1+arad);

% Now construct the real root of the cubic equation
term1 = -2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2;
term2 = (term1 + sqrt(4*(-a2^2 + 3*a1*a3)^3 + term1^2))^(1/3);

y = -a2/(3*a3) - 2^(1/3)*(-a2^2 + 3*a1*a3)/(3*a3*term2) + term2/(3*a3*2^(1/3));