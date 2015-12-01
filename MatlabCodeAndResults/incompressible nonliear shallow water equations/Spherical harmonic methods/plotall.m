% Script nonlinsolve...used to solve for the nonliear incompressible
% shallow water coefficients

% Clear all constants from memory
clear all

% Declare all constants as global variables

global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp

% Read in the parameters from the linearised solution
params=dlmread('linparams.txt','\r');

N=params(1,1); % The linearised truncation level
kappa=params(1,2);% Number of wavelengths around a latitude circle
h0=params(1,3);% Polar free surface height
w=params(1,4);% Base zonal flow angular velocity (user defined)
Amplin=params(1,5);% The amplitude of the R-H wave at latitude 45N

% Declare all constants
M = N; % The nonlinear series truncation level
numunknowns = 3*M^2+1; % The number of unknowns in the nonlinear model
Amp = Amplin; % The nonliear amplitude condition
g=9.80616; % Gravitational acceleration
Omega=2*pi/(24*60*60); % Earth's angular velocity
a=6.37122*10^6; % Radius of the Earth

% Set the number of contours for plotting
numcontours=20;

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

typX=coeffs;
for i=1:numunknowns
    if typX(i,1)==0
        typX(i,1)=1;
    end
end

% Parameters for the optimisation process
opts = struct('Display','testing','LargeScale','on', ...
   'TolX',1e-8,'TolFun',1e-9,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX',typX,...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','LevenbergMarquardt','off'); 

%'ones(numberOfVariables,1)',...
%X=coeffs;

%X=FSOLVE('residual',coeffs,optimset('Display','off'),phi,P,Pp,C,S,Pamp,M,numunknowns);
%X=FSOLVE('residual',coeffs,optimset('Display','testing'));
%[X,FVAL,EXITFLAG,OUTPUT,JACOB]=FSOLVE('residual',coeffs,optimset(opts));


%ntrial=0;
%tolx=1e-8;
%tolf=1e-8;
%[X,errf,errx] = mnewt(ntrial,coeffs,tolx,tolf)

%define the wavespeed c
c=X(numunknowns,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start the plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The plotting grid points
eta = 0:2*pi/150:2*pi;
phi = [0:29*pi/(180*5):29*pi/180 pi/6:40*pi/(180*30):70*pi/180 71*pi/180:19*pi/(180*5):pi/2];

% Set up the coefficient vectors
A = X(1:M*M,1);
B = X(M*M+1:2*M*M,1);
C = X(2*M*M+1:3*M*M,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the pressure field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the c++ function to calculate the free surface height at all grid points
hfield=fsfield(eta,phi,C,kappa,w,a,g,Omega,h0,M);

eatmat=zeros(max(size(phi)),max(size(eta)));
phimat=zeros(max(size(phi)),max(size(eta)));

[rows,cols]=size(hfield);

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
Ns = num2str(N);
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
t3 = ', and N=';
t2 = strcat(t2,t3,Ns);
t1 = 'Polar stereographic height contours';
tit=strcat(t1,t2);
title(tit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot with grid lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
contour(xgrid,ygrid,hfield,numcontours)
COLORBAR('horiz')
axis off
title(tit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add some grid lines to the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add some grid lines...
res=200;
xcirc=zeros(1,res);
ycirc=zeros(res,1);
theta=0:2*pi/(res-1):2*pi;

for j=0:8
    phival=j*pi/18;
    for i=1:res
        xcirc(1,i)=cos(theta(1,i))*cos(phival)/(1+sin(phival));
        ycirc(i,1)=sin(theta(1,i))*cos(phival)/(1+sin(phival));
    end
    plot(xcirc,ycirc,'k')
end

% phi=pi/3
phival=pi/4;
for i=1:res
    xcirc(1,i)=cos(theta(1,i))*cos(phival)/(1+sin(phival));
    ycirc(i,1)=sin(theta(1,i))*cos(phival)/(1+sin(phival));
end
plot(xcirc,ycirc,'k--')

res=200;
xcirc=zeros(1,res);
ycirc=zeros(res,1);
phival=0:pi/(2*(res-1)):pi/2;

for j=0:18
    thetaval=j*2*pi/18;
    for i=1:res
        xcirc(1,i)=cos(thetaval)*cos(phival(1,i))/(1+sin(phival(1,i)));
        ycirc(i,1)=sin(thetaval)*cos(phival(1,i))/(1+sin(phival(1,i)));
    end
    plot(xcirc,ycirc,'k')
    plot(xcirc,-ycirc,'k')
end
for i=0:9;
    degrees=num2str(i*10);
    text(cos(i*pi/18)*cos(151*pi/90)/(1+sin(i*pi/18)),cos(i*pi/18)*sin(151*pi/90)/(1+sin(i*pi/18)),degrees);
end

for i=0:17
    degrees=num2str(i*20);
    if i<6
        text(1.06*cos(20*i*pi/180),1.06*sin(20*i*pi/180),degrees);
    elseif i==6
        text(1.1*cos(20*i*pi/180),1.1*sin(20*i*pi/180),degrees);        
    elseif (i<12 & i~=6)
        text(1.15*cos(20*i*pi/180),1.15*sin(20*i*pi/180),degrees);
    elseif i==12
        text(1.1*cos(20*i*pi/180),1.1*sin(20*i*pi/180),degrees);
    else
        text(1.06*cos(20*i*pi/180),1.06*sin(20*i*pi/180),degrees);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the velocity vector field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the c++ function to calculate the velocity components at all grid points
ulamfield=ulam(eta,phi,A,kappa,w,a,g,Omega,h0,M);
uphifield=uphi(eta,phi,B,kappa,w,a,g,Omega,h0,M);

eatmat=zeros(max(size(phi)),max(size(eta)));
phimat=zeros(max(size(phi)),max(size(eta)));

[rows,cols]=size(ulamfield);

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

U=ulamfield;
V=uphifield;
for i=1:rows
    for j=1:cols
        U(i,j)=-sin(etamat(i,j))*ulamfield(i,j)-cos(etamat(i,j))*uphifield(i,j);
        V(i,j)=cos(etamat(i,j))*ulamfield(i,j)-sin(etamat(i,j))*uphifield(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot without grid lines first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
quiver(xgrid,ygrid,U,V);

t1 = 'Polar stereographic horizontal velocity vector field ';
tit=strcat(t1,t2);
title(tit)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot with grid lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
quiver(xgrid,ygrid,U,V);
title(tit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add some grid lines to the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res=200;
xcirc=zeros(1,res);
ycirc=zeros(res,1);
theta=0:2*pi/(res-1):2*pi;

for j=0:8
    phival=j*pi/18;
    for i=1:res
        xcirc(1,i)=cos(theta(1,i))*cos(phival)/(1+sin(phival));
        ycirc(i,1)=sin(theta(1,i))*cos(phival)/(1+sin(phival));
    end
    plot(xcirc,ycirc,'k')
end

% phi=pi/3
phival=pi/4;
for i=1:res
    xcirc(1,i)=cos(theta(1,i))*cos(phival)/(1+sin(phival));
    ycirc(i,1)=sin(theta(1,i))*cos(phival)/(1+sin(phival));
end
plot(xcirc,ycirc,'k--')

res=200;
xcirc=zeros(1,res);
ycirc=zeros(res,1);
phival=0:pi/(2*(res-1)):pi/2;

for j=0:18
    thetaval=j*2*pi/18;
    for i=1:res
        xcirc(1,i)=cos(thetaval)*cos(phival(1,i))/(1+sin(phival(1,i)));
        ycirc(i,1)=sin(thetaval)*cos(phival(1,i))/(1+sin(phival(1,i)));
    end
    plot(xcirc,ycirc,'k')
    plot(xcirc,-ycirc,'k')
end
for i=0:9;
    degrees=num2str(i*10);
    text(cos(i*pi/18)*cos(151*pi/90)/(1+sin(i*pi/18)),cos(i*pi/18)*sin(151*pi/90)/(1+sin(i*pi/18)),degrees);
end

for i=0:17
    degrees=num2str(i*20);
    if i<6
        text(1.06*cos(20*i*pi/180),1.06*sin(20*i*pi/180),degrees);
    elseif i==6
        text(1.1*cos(20*i*pi/180),1.1*sin(20*i*pi/180),degrees);        
    elseif (i<12 & i~=6)
        text(1.15*cos(20*i*pi/180),1.15*sin(20*i*pi/180),degrees);
    elseif i==12
        text(1.1*cos(20*i*pi/180),1.1*sin(20*i*pi/180),degrees);
    else
        text(1.06*cos(20*i*pi/180),1.06*sin(20*i*pi/180),degrees);
    end
end

%%%%%%%%%%%%%%%%%%%%%
% Plot both pressure and velocity field
%%%%%%%%%%%%%%%%%%%%%
figure
axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
contour(xgrid,ygrid,hfield,numcontours)
COLORBAR('horiz')
axis off
quiver(xgrid,ygrid,U,V);
t1 = 'Polar stereographic height and horizontal velocity vector field ';
tit=strcat(t1,t2);
title(tit)