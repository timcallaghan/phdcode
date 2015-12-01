% Script nonlinsolve...used to solve for the nonliear incompressible
% shallow water coefficients. We include all possible basis functions
% and the only synmmetry assumed is that of latitudinal symmetry for 
% ulam and h and latitudinal anti-symmetry for uphi. This ammounts
% to every second coefficient being zero in each series representation.

% Clear all constants from memory
clear all

% Declare all constants as global variables
global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

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
M = kappa+1; % The nonlinear series truncation level must be > kappa if we want to use the linear solution as a starting guess
numunknowns = 6*M^2; % The number of unknowns in the nonlinear model
Amp = Amplin; % The nonliear Nondimensional amplitude condition
g=9.80616; % Gravitational acceleration
Omega=2*pi/(24*60*60); % Earth's angular velocity
a=6.37122*10^6; % Radius of the Earth


% Read in the coefficients calculated from the linear model
lincoeffs=dlmread('lincoeffs.txt','\r');

% Define the initial guess for the coeffs from the linear model.
% Note that this depends on the value of kappa used so that's why it
% appears here.
coeffs=zeros(numunknowns+1,1);
if (M<N)
    % Check if we have more linear coeffs that we need.
    for i=1:M
        % Note that the linear solution sets up parts of A1
        % Also we have to allow for the m=0 terms in each series.
        coeffs(kappa*M+i,1)=lincoeffs(1,i);
        % Note that the linear solution sets up parts of B2 (not B1)
        coeffs(3*M*M+(kappa-1)*M+i,1)=lincoeffs(1,i+N);
        %coeffs(3*M*M+kappa*M+i,1)=lincoeffs(1,i+N);
        % Note that the linear solution sets up parts of C1
        coeffs(4*M*M+kappa*M+i,1)=lincoeffs(1,i+2*N);
    end
else
    % We don't have more linear coeffs than we need so just fill
    % out as per normal.
    for i=1:N
        % Note that the linear solution sets up parts of A1
        % Also we have to allow for the m=0 terms in each series.
        coeffs(kappa*M+i,1)=lincoeffs(1,i);
        % Note that the linear solution sets up parts of B2 (not B1)
        coeffs(3*M*M+(kappa-1)*M+i,1)=lincoeffs(1,i+N);
        %coeffs(3*M*M+kappa*M+i,1)=lincoeffs(1,i+N);
        % Note that the linear solution sets up parts of C1
        coeffs(4*M*M+kappa*M+i,1)=lincoeffs(1,i+2*N);
    end
end


% Declare which coeff/parameter we will hold fixed
% Hold the nonlinear wavespeed fixed
fixed=numunknowns+1
% Hold the first coeff of the uphi velocity fixed
%fixed=kappa*M+1;
% Adjust the nonlinear wavespeed so that it is slightly faster than the
% value calculated from the linearised theory (to make sure we are hopefully in a solution region)
waveeps=0.05;
coeffs(numunknowns+1,1)=(1.0+waveeps)*lincoeffs(1,3*N+1);

%%%%%%%%%%%%%TESTING%%%%%%%%%%%%%%%%%%%%%%%%
% For testing purposes we set all coeffs to 1
%coeffs=0.01*ones(numunknowns+1,1);
%coeffs=zeros(numunknowns+1,1);
%%%%%%%%%%%%%TESTING%%%%%%%%%%%%%%%%%%%%%%%%


% Declare the collocation points in mu and phi.
% Get the mu values using the Gauss-Legendre abcissas.
mu = gauleg(-1,1,M);
% Since mu=sin(phi) we can invert to get the phi points.
phi = asin(mu);

% Cache the basis functions
C=cacheC(M);
S=cacheS(M);
P=cacheL(mu,M,kappa);
%Note that cacheLp accepts phi (not mu) as argument!!!
Pp=cacheLp(phi,M,kappa);

% Used to test the initial value of the residual vector.
%X=residual(coeffs);
% Used to test the initial value of the jacobian matrix.
%jac=fdjac(coeffs,X);
% Used to write the jacobian to file for analysis.
%dlmwrite('fdjacob.txt',jac,'\t')

% The maximum number of interations to perform in the
% Newton root finding method
ntrial=10;
% The convergence error tolerance on the coefficients vector.
tolx=1e-8;
% The convergence error tolerance on the residual vector.
tolf=1e-8;
% Calls the Newton method and returns the final coefficient
% vector X as well as the exit error values.
[X,errf,errx] = mnewt(ntrial,coeffs,tolx,tolf);


% for the sake of testing we set X=coeffs
%X=coeffs;

%%%%%%%%%%PLOTTING ROUTINE%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start the plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the number of contours for plotting
numcontours=20;

% The plotting grid points
eta = 0:2*pi/150:2*pi;
phi = [0:29*pi/(180*5):29*pi/180 pi/6:40*pi/(180*30):70*pi/180 71*pi/180:19*pi/(180*5):pi/2];

% Set up the coefficient vectors
A1 = X(1:M*M,1);
A2 = X(M*M+1:2*M*M,1);
B1 = X(2*M*M+1:3*M*M,1);
B2 = X(3*M*M+1:4*M*M,1);
C1 = X(4*M*M+1:5*M*M,1);
C2 = X(5*M*M+1:6*M*M,1);
%define the wavespeed c
c=X(numunknowns+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the free surface height field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the c++ function to calculate the free surface height at all grid points
hfield=fsfield(eta,phi,C1,C2,w,Fr,Ro,h0,M);

eatmat=zeros(max(size(phi)),max(size(eta)));
phimat=zeros(max(size(phi)),max(size(eta)));

[rows,cols]=size(hfield);

% Set up the latitude-longitude mesh grid points
for i=1:cols;
phimat(:,i)=phi';
end

for i=1:rows
etamat(i,:)=eta;
end

xgrid=etamat;
ygrid=phimat;
for i=1:rows
    for j=1:cols
        xgrid(i,j)=cos(phimat(i,j))*cos(etamat(i,j))/(1+sin(phimat(i,j)));
        ygrid(i,j)=cos(phimat(i,j))*sin(etamat(i,j))/(1+sin(phimat(i,j)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot without grid lines first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
contour(xgrid,ygrid,hfield,numcontours)
COLORBAR('horiz')
axis off

ws = num2str(w);
kappas = num2str(kappa);
Ms = num2str(M);
h0s = num2str(h0);
cs = num2str(c);

t3 = ' with w=';
t2 = strcat(t3,ws);
t3 = ', kappa=';
t2 = strcat(t2,t3,kappas);
t3 = ', h0=';
t2 = strcat(t2,t3,h0s);
t3 = ', c=';
t2 = strcat(t2,t3,cs);
t3 = ', and M=';
t2 = strcat(t2,t3,Ms);
t1 = 'Polar stereographic height contours';
tit=strcat(t1,t2);
title(tit)