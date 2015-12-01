% Script to calculate the L.H.S matrix A in the 
% generalized eigenvalue problem for the Linearized
% Incompressible shallow water equations.

% Define all global variables
global N kappa w MATA h0 Sr Fr Ro

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the entries for the lambda-momentum equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first equation
MATA(1,1)=-kappa*w/2;
MATA(1,N+1)=-1/4*(1/Ro+2*w);
MATA(1,2*N+1)=kappa/Fr^2;
% The second equation
MATA(2,1)=-kappa*w/2;
MATA(2,2)=-kappa*w/2;
MATA(2,N+2)=-1/4*(1/Ro+2*w);
MATA(2,2*N+1)=kappa/Fr^2;
MATA(2,2*N+2)=-kappa/Fr^2;
% Now do the other N-2 equations using the formula
for i=2:N-1
    MATA(i+1,i)=-kappa*w/2;
    MATA(i+1,i+1)=-kappa*w/2;
    MATA(i+1,N+i-1)=1/4*(1/Ro+2*w);
    MATA(i+1,N+i+1)=-1/4*(1/Ro+2*w);
    MATA(i+1,2*N+i)=-kappa/(Fr^2)*(-1)^i;
    MATA(i+1,2*N+i+1)=kappa/(Fr^2)*(-1)^i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the entries for the phi-momentum equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first N-1 equations
for i=1:N-1
    MATA(N+i,i)=1/2*(1/Ro+2*w);
    MATA(N+i,i+1)=-1/2*(1/Ro+2*w);
    MATA(N+i,N+i)=kappa*w;
    MATA(N+i,2*N+i)=2*i/(Fr^2)*(-1)^(i+1);
    MATA(N+i,2*N+i+1)=-2*i/(Fr^2)*(-1)^(i+1);
end
% The last equation
MATA(2*N,N)=1/2*(1/Ro+2*w);
MATA(2*N,2*N)=kappa*w;
MATA(2*N,3*N)=2*N/(Fr^2)*(-1)^(N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the entries for the mass equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first equation
MATA(2*N+1,1)=-kappa*h0-3/8*kappa*w*Fr^2*(1/Ro+w);
MATA(2*N+1,2)=-1/8*kappa*w*Fr^2*(1/Ro+w);
MATA(2*N+1,N+1)=h0/2+1/8*w*Fr^2*(1/Ro+w);
MATA(2*N+1,N+2)=1/16*w*Fr^2*(1/Ro+w);
MATA(2*N+1,2*N+1)=3*kappa*w/2;
MATA(2*N+1,2*N+2)=-kappa*w/2;
% The second equation
MATA(2*N+2,1)=-1/8*kappa*w*Fr^2*(1/Ro+w);
MATA(2*N+2,2)=-kappa*h0-1/4*kappa*w*Fr^2*(1/Ro+w);
MATA(2*N+2,3)=-1/8*kappa*w*Fr^2*(1/Ro+w);
MATA(2*N+2,N+1)=3*h0/2+9/16*w*Fr^2*(1/Ro+w);
MATA(2*N+2,N+2)=3*h0/2+9/16*w*Fr^2*(1/Ro+w);
MATA(2*N+2,N+3)=3/16*w*Fr^2*(1/Ro+w);
MATA(2*N+2,2*N+1)=kappa*w/2;
MATA(2*N+2,2*N+2)=-kappa*w;
MATA(2*N+2,2*N+3)=kappa*w/2;
% Now do the next N-3 equations using the formula
for i=3:N-1
    MATA(2*N+i,i-1)=-1/8*kappa*w*Fr^2*(1/Ro+w);
    MATA(2*N+i,i)=-kappa*h0-1/4*kappa*w*Fr^2*(1/Ro+w);
    MATA(2*N+i,i+1)=-1/8*kappa*w*Fr^2*(1/Ro+w);
    MATA(2*N+i,N+i-2)=1/16*w*Fr^2*(1/Ro+w)*(2*i-1);
    MATA(2*N+i,N+i-1)=(2*i-1)*h0/2+3/16*w*Fr^2*(1/Ro+w)*(2*i-1);
    MATA(2*N+i,N+i)=(2*i-1)*h0/2+3/16*w*Fr^2*(1/Ro+w)*(2*i-1);
    MATA(2*N+i,N+i+1)=1/16*w*Fr^2*(1/Ro+w)*(2*i-1);
    MATA(2*N+i,2*N+i-1)=kappa*w/2*(-1)^i;
    MATA(2*N+i,2*N+i)=-kappa*w*(-1)^i;
    MATA(2*N+i,2*N+i+1)=kappa*w/2*(-1)^i;
end
% The last equation
MATA(3*N,N-1)=-1/8*kappa*w*Fr^2*(1/Ro+w);
MATA(3*N,N)=-kappa*h0-1/4*kappa*w*Fr^2*(1/Ro+w);
MATA(3*N,2*N-2)=1/16*w*Fr^2*(1/Ro+w)*(2*N-1);
MATA(3*N,2*N-1)=(2*N-1)*h0/2+3/16*w*Fr^2*(1/Ro+w)*(2*N-1);
MATA(3*N,2*N)=(2*N-1)*h0/2+3/16*w*Fr^2*(1/Ro+w)*(2*N-1);
MATA(3*N,3*N-1)=kappa*w/2*(-1)^N;
MATA(3*N,3*N)=-kappa*w*(-1)^N;
