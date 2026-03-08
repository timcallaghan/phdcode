function masserr = masserror(h0,w,a,gamma,Fr,Ma,Ro,Mb)
%MASSERROR Calculates the error in the mass specification
% equation for given values of w and h0.

% Calculate some constants for the problem
B=(gamma-1)/gamma*(Ma/Fr)^2;
Aw=w*Fr^2/2.0*(1.0/Ro+w);
% Calculate the total mass of the compressible sytem.
% We use adaptive recursive quadrature here (Gauss-Lobatto quadrature
% with Kronrod extension)
Mz=4.0*pi*(gamma-1.0)*a^2/((B^3)*gamma*(2.0*gamma-1)*(3.0*gamma-2.0))*adaptlob(@integrnd,0,pi/2,1.0*10^(-1),[],h0,B,Aw,a,gamma);

% Now we form the error expression for the function to be zeroed.

masserr = Mz/Mb - 1.0;
