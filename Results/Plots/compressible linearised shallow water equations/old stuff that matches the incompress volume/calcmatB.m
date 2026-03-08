% Script to calculate the R.H.S matrix B in the 
% generalized eigenvalue problem for the Linearized
% compressible shallow water equations.

% Define all global variables
global N kappa w MATB h0 Sr Fr Ro

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the entries for the lambda-momentum equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first equation
MATB(1,1)=-kappa*Sr/2;
% The second equation
MATB(2,1)=-kappa*Sr/2;
MATB(2,2)=-kappa*Sr/2;
% Now do the other N-2 equations using the formula
for i=2:N-1
    MATB(i+1,i)=-kappa*Sr/2;
    MATB(i+1,i+1)=-kappa*Sr/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the entries for the phi-momentum equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We can do all N equations in one go, so off we go!
for i=1:N
    MATB(N+i,N+i)=kappa*Sr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the entries for the mass equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first equation
MATB(2*N+1,2*N+1)=3*kappa*Sr/2;
MATB(2*N+1,2*N+2)=-kappa*Sr/2;
% The second equation
MATB(2*N+2,2*N+1)=kappa*Sr/2;
MATB(2*N+2,2*N+2)=-kappa*Sr;
MATB(2*N+2,2*N+3)=kappa*Sr/2;
% Now do the next N-3 equations using the formula
for i=3:N-1
    MATB(2*N+i,2*N+i-1)=kappa*Sr/2*(-1)^i;
    MATB(2*N+i,2*N+i)=-kappa*Sr*(-1)^i;
    MATB(2*N+i,2*N+i+1)=kappa*Sr/2*(-1)^i;   
end
% The last equation
MATB(3*N,3*N-1)=kappa*Sr/2*(-1)^N;
MATB(3*N,3*N)=-kappa*Sr*(-1)^N;