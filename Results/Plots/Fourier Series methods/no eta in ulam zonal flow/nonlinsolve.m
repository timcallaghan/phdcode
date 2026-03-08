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

% The double integration tolerance
intol=10^(-14);

% Read in the coefficients calculated from the linear model
lincoeffs=dlmread('lincoeffs.txt','\r');

% Define the initial guess for the coeffs from the linear model
coeffs=zeros(numunknowns+1,1);
% Assign all the linearized starting guess coeffs
for i=1:N
    coeffs(i,1)=lincoeffs(1,i);
	coeffs(M*N+i,1)=lincoeffs(1,i+N);
	coeffs(2*M*N+N+1+i,1)=lincoeffs(1,i+2*N);
end
% Assign the zonal flow structure to h
coeffs(2*M*N+1,1)=h0+w*Fr^2*(1/Ro+w)/4;
coeffs(2*M*N+2,1)=w*Fr^2*(1/Ro+w)/4;
%define the wavespeed c
coeffs(numunknowns+1,1)=lincoeffs(1,3*N+1);

% Declare which coeff/parameter we will hold fixed
fixed=2*M*N+N+2;
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

% Used for testing...
%X=residual(coeffs);
%jac=jacobian(coeffs);
%dlmwrite('jacobian.txt',jac,'\t');
%return

% Do the Newton iterations
maxhalvings=30;
ntrial=10;
tolx=1e-12;
tolf=1e-12;
[X,errf,errx] = mnewt(ntrial,coeffs,tolx,tolf);

% Write the coefficients and the wavespeed c to file for use in the plotting routine
% First clear the file contents by opening and closing without appending (ie writing)
fid=fopen('nlincoeffs.txt','wt');
fclose(fid);
% Now open again and append all the constants and variables
fid=fopen('nlincoeffs.txt','at');
for i=1:numunknowns+1
    fprintf(fid,'%.16e',X(i,1));
    fprintf(fid,'\r');
end
fclose(fid);

% Now write all the parameters to file for use in data analysis
% First clear the file contents by opening and closing without appending (ie writing)
fid=fopen('nlinparams.txt','wt');
fclose(fid);
% Now open again and append all the constants and variables
fid=fopen('nlinparams.txt','at');
fprintf(fid,'%d',N);
fprintf(fid,'\r');
fprintf(fid,'%d',kappa);
fprintf(fid,'\r');
fprintf(fid,'%.16e',h0);
fprintf(fid,'\r');
fprintf(fid,'%.16e',w);
fprintf(fid,'\r');
fprintf(fid,'%d',M);
fprintf(fid,'\r');
fprintf(fid,'%.16e',Sr);
fprintf(fid,'\r');
fprintf(fid,'%.16e',Fr);
fprintf(fid,'\r');
fprintf(fid,'%.16e',Ro);
fprintf(fid,'\r');
fprintf(fid,'%.16e',href);
fprintf(fid,'\r');
fprintf(fid,'%.16e',cref);
fprintf(fid,'\r');
fprintf(fid,'%.16e',vref);
fprintf(fid,'\r');
fprintf(fid,'%.16e',wref);
fprintf(fid,'\r');
fprintf(fid,'%.16e',intol);
fprintf(fid,'\r');
fprintf(fid,'%.16e',g);
fprintf(fid,'\r');
fprintf(fid,'%.16e',Omega);
fprintf(fid,'\r');
fprintf(fid,'%.16e',aref);
fprintf(fid,'\r');
fprintf(fid,'%.16e',a);
fprintf(fid,'\r');
fprintf(fid,'%d',fixed);
fprintf(fid,'\r');
fprintf(fid,'%.16e',tolx);
fprintf(fid,'\r');
fprintf(fid,'%.16e',tolf);
fprintf(fid,'\r');
fprintf(fid,'%.16e',errx);
fprintf(fid,'\r');
fprintf(fid,'%.16e',errf);
fprintf(fid,'\r');
fclose(fid);