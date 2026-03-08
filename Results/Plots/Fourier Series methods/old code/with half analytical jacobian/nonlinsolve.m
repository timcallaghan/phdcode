% Script nonlinsolve...used to solve for the nonliear incompressible
% shallow water coefficients using Fourier Series

% Clear all constants from memory
clear all

% Declare all constants as global variables
global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed Vzon intol a maxhalvings

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
Vzon=params(1,13); % Hemispherical volume for the zonal flow used in the Linearized model
g=params(1,14); % Gravitational acceleration
Omega=params(1,15); % Earth's angular velocity
aref=params(1,16); % Radius of the Earth
a=params(1,17); % Dimensionless value of the Earth's radius...relative to href...

% Declare all constants
M = N; % The nonlinear longitude series truncation level
numunknowns = 3*M*N+1; % The number of unknowns in the nonlinear model
Amp = Amplin; % The nonlinear Nondimensional amplitude condition

% The double integration tolerance
intol=10^(-14);

% Read in the coefficients calculated from the linear model
lincoeffs=dlmread('lincoeffs.txt','\r');

% Define the initial guess for the coeffs from the linear model
coeffs=zeros(numunknowns+1,1);
% Assign all the linearized starting guess coeffs
for i=1:N
    coeffs(N+i,1)=lincoeffs(1,i);
	coeffs(M*N+i,1)=lincoeffs(1,i+N);
	coeffs(2*M*N+N+1+i,1)=lincoeffs(1,i+2*N);
end
% Assign the zonal flow structure to h
coeffs(2*M*N+1,1)=h0+w*Fr^2*(1/Ro+w)/4;
coeffs(2*M*N+2,1)=w*Fr^2*(1/Ro+w)/4;
%define the wavespeed c
coeffs(numunknowns+1,1)=lincoeffs(1,3*N+1);

% Declare which coeff/parameter we will hold fixed
fixed=3*M*N+2;
% If we are holding the wavespeed fixed we should make it slightly larger
% than the linearized value
coeffs(fixed,1)=1.0*coeffs(fixed,1);
%coeffs(fixed,1)=0.0;

% Declare the collocation points in phi
delphi=pi/(2*N);
epsphi=delphi/2;
phi=zeros(N,1);
for i=1:N
    phi(i,1)=(i-1)*delphi+epsphi;
end

% Cache the basis functions
Ceta=cacheCeta(N,kappa);
Seta=cacheSeta(N,kappa);
Cphi1=cacheCphi1(N,phi);
Sphi1=cacheSphi1(N,phi);
Cphi2=cacheCphi2(N,phi);
Sphi2=cacheSphi2(N,phi);

% Calculate the volume of the fluid.
%intol=10^(-10);
%intol=[];
%Volume2=dblquad(@intgrand,0,pi/2,0,2*pi,intol,@adaptlob,coeffs(2*M*N+1:numunknowns,1),M,N,kappa,a);


%tic
%X=residual(coeffs);
%epsilon=10^(-5);
%coeffs(1,1)=coeffs(1,1)+epsilon;
%X2=residual(coeffs);
%(X2(1,1)-X1(1,1))/epsilon
%tic
%jac=fdjac(coeffs,X);
%toc
%toc
%dlmwrite('fdjacob.txt',jac,'\t')

%tic
df = jacobian(coeffs);
[U,S,V] = svd(df);
%cond(df)
%toc
%dlmwrite('analjacob.txt',df,'\t')

%X=residual(coeffs);
%jac=fdjac(coeffs,X);
%dlmwrite('fdjacob.txt',jac,'\t')

%jacdiff=abs(df-jac);
%dlmwrite('jacdiff.txt',jacdiff,'\t')

maxhalvings=30;
ntrial=1;
tolx=1e-10;
tolf=1e-10;
%[X,errf,errx] = mnewt(ntrial,coeffs,tolx,tolf);