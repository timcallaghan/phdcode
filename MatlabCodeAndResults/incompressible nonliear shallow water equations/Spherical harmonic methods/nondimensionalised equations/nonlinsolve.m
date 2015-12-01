% Script nonlinsolve...used to solve for the nonliear incompressible
% shallow water coefficients

% Clear all constants from memory
clear all

% Declare all constants as global variables

global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro

% Read in the parameters from the linearised solution
params=dlmread('linparams.txt','\r');

N=params(1,1); % The linearised truncation level
kappa=params(1,2);% Number of wavelengths around a latitude circle
h0=params(1,3);% Nondimensional Polar free surface height
w=params(1,4);% Nondimensional Base zonal flow angular velocity (user defined)
Amplin=params(1,5);% The Nondimensional amplitude of the R-H wave at latitude 45N
Sr=params(1,6); % Strouhal number
Fr=params(1,7); % Froude number
Ro=params(1,8); % Rossby number
href=params(1,9); % Characteristic length scale (m)
cref=params(1,10);% Characterisitc angular velocity scale of the Rossby wave (s^-1)
vref=params(1,11);% Characterisitc velocity scale(m*s^-1)
wref=params(1,12);%Characterisitc angular velocity scale of the zonal flow (s^-1)

% Declare all constants
M = N; % The nonlinear series truncation level
numunknowns = 3*M^2+1; % The number of unknowns in the nonlinear model
Amp = Amplin; % The nonliear Nondimensional amplitude condition
g=9.80616; % Gravitational acceleration
Omega=2*pi/(24*60*60); % Earth's angular velocity
a=6.37122*10^6; % Radius of the Earth

% Read in the coefficients calculated from the linear model
lincoeffs=dlmread('lincoeffs.txt','\r');

A=zeros(N,1);
B=zeros(N,1);
C=zeros(N,1);
for i=1:N/2
    A(2*i-1,1)=lincoeffs(1,i);
    B(2*i,1)=lincoeffs(1,N+i);
    C(2*i-1,1)=lincoeffs(1,2*N+i);
end

% Define the initial guess for the coeffs from the linear model
coeffs=zeros(numunknowns,1);
if (M<N)
    for i=1:M
		coeffs(i,1)=lincoeffs(1,i);
		coeffs(M*M+i,1)=lincoeffs(1,i+N);
		coeffs(2*M*M+i,1)=lincoeffs(1,i+2*N);
    end
else
	for i=1:N
		coeffs(i,1)=lincoeffs(1,i);
		coeffs(M*M+i,1)=lincoeffs(1,i+N);
		coeffs(2*M*M+i,1)=lincoeffs(1,i+2*N);
    end
end
%define the wavespeed c
coeffs(numunknowns,1)=lincoeffs(1,3*N+1);

% Declare the collocation points in mu and phi
mu = gauleg(-1,1,M);
phi = asin(mu);

% Cache the basis functions
C=cacheC(M,kappa);
S=cacheS(M,kappa);
P=cacheL(mu,M,kappa);
Pamp=cacheLamp(M,kappa);
%Note that Lp accepts phi (not mu) as argument!!!
Pp=cacheLp(phi,M,kappa);


% Parameters for the optimisation process
opts = struct('Display','testing','LargeScale','on', ...
   'TolX',1e-8,'TolFun',1e-8,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','LevenbergMarquardt','off'); 

%X=FSOLVE('residual',coeffs,optimset('Display','off'),phi,P,Pp,C,S,Pamp,M,numunknowns);
%X=FSOLVE('residual',coeffs,optimset('Display','testing'));
%[X,FVAL,EXITFLAG,OUTPUT,JACOB]=FSOLVE('residual',coeffs,optimset(opts));

%tic
%X=residual(coeffs);
%jac=fdjac(coeffs,X);
%toc

ntrial=5;
tolx=1e-4;
tolf=1e-4;
[X,errf,errx] = mnewt(ntrial,coeffs,tolx,tolf)