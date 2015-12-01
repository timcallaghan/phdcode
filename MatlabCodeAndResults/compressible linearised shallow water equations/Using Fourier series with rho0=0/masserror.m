function masserr = masserror(h0,w,a,gamma,Fr,Ma,Ro,Mbp)
%MASSERROR Calculates the error in the mass specification
% equation for given values of w and h0.

% Calculate some constants for the problem
B=(gamma-1)/gamma*(Ma/Fr)^2;
Aw=w*Fr^2/2.0*(1.0/Ro+w);
% Calculate the total mass of the compressible system.
% We use adaptive recursive quadrature here (Gauss-Lobatto quadrature
% with Kronrod extension)
MzonMb=(1.0/Mbp)*adaptlob(@integrnd,0,pi/2,[],[],h0,B,Aw,a,gamma);

% Now we form the error expression for the function to be zeroed.

masserr = MzonMb - 1.0;
