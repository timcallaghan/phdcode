% Script nonlinsolve...used to solve for the nonliear incompressible
% shallow water coefficients using Fourier Series

% Clear all constants from memory
clear all

% Declare all constants as global variables
global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed

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
M = N; % The nonlinear longitude series truncation level
numunknowns = 3*M*N; % The number of unknowns in the nonlinear model
Amp = Amplin; % The nonliear Nondimensional amplitude condition
g=9.80616; % Gravitational acceleration
Omega=2*pi/(24*60*60); % Earth's angular velocity
a=6.37122*10^6; % Radius of the Earth

% Read in the coefficients calculated from the linear model
lincoeffs=dlmread('lincoeffs.txt','\r');

% Define the initial guess for the coeffs from the linear model
coeffs=zeros(numunknowns+1,1);
% Assign all the linearized starting guess coeffs
for i=1:N
    coeffs(N+i,1)=lincoeffs(1,i);
	coeffs(M*N+i,1)=lincoeffs(1,i+N);
	coeffs(2*M*N+N+i,1)=lincoeffs(1,i+2*N);
end
%define the wavespeed c
coeffs(numunknowns+1,1)=lincoeffs(1,3*N+1);

% Declare which coeff/parameter we will hold fixed
fixed=numunknowns+1;
% If we are holding the wavespeed fixed we should make it slightly larger
% than the linearized value
coeffs(numunknowns+1,1)=0.95*coeffs(numunknowns+1,1);

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


%tic
%X=residual(coeffs);
%epsilon=10^(-5);
%coeffs(1,1)=coeffs(1,1)+epsilon;
%X2=residual(coeffs);
%(X2(1,1)-X1(1,1))/epsilon
%jac=fdjac(coeffs,X);
%toc
%dlmwrite('fdjacob.txt',jac,'\t')

%df = jacobian(coeffs);
%dlmwrite('analjacob.txt',df,'\t')

%X=residual(coeffs);
%jac=fdjac(coeffs,X);
%dlmwrite('fdjacob.txt',jac,'\t')

%jacdiff=abs(df-jac);
%dlmwrite('jacdiff.txt',jacdiff,'\t')

ntrial=10;
tolx=1e-4;
tolf=1e-4;
[X,errf,errx] = mnewt(ntrial,coeffs,tolx,tolf);