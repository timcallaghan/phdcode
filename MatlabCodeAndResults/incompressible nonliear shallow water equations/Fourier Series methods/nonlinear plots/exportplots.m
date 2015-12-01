% Script to plot up the nonlinear solutions found
% using the c++ code

% Clear all constants from memory
clear all

% Declare and initialise all constants
g=9.80616; % Gravitational acceleration
Omega=2*pi/(24*60*60); % Earth's angular velocity
aref=6.37122*10^6; % Radius of the Earth

% Declare characteristic scales
href=8*10^3; % Characteristic length scale (m)
cref=Omega/30; % Characteristic angular velocity scale of the Rossby wave (s^-1)
vref=40; % Characteristic velocity scale(m*s^-1)
wref=7.848*10^(-6); %Characteristic angular velocity scale of the zonal flow (s^-1)

% Define the nondimensional numbers
Sr=aref*cref/vref; % Strouhal number
Fr=vref/sqrt(g*href); % Froude number
Ro=vref/(2*Omega*aref); % Rossby number

% Dimensionless zonal flow parameters
wbase=wref*aref/vref; % Dimensionless base zonal flow angular velocity
h0=1; % Dimensionless polar free suface height (h/href=1 is assumed at the poles)
a=aref/href; % Dimensionless value of the Earth's radius...relative to href...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Read in all the parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the truncation levels used
MN=dlmread('oldparams.txt','\r');
% Assign the parameters
M=MN(1,1); % The nonlinear longitudinal truncation.
N=MN(2,1); % The nonlinear latitudinal truncation.

% Get the other model parameters (only read first 5 rows)
params=dlmread('linparams.txt','\r',[0 0 4 0]);
% Assign all the parameters
kappa=params(1,1); % Number of wavelengths around a latitude circle
wpercent=params(3,1); % Nondimensional Base zonal flow angular velocity (as a percentage of w defined above)
% Get h0
h0 = params(4,1);

% Modify the value of w based on wpercent
w=wpercent*wbase;

% Read in the coefficients calculated from the nonlinear model
nlincoeffs=dlmread('current.txt',' ');
% Assign the coeffs to their respective vectors
P=-nlincoeffs(1,1:M*N)';
Q=nlincoeffs(1,M*N+1:2*M*N)';
H=nlincoeffs(1,2*M*N+1:3*M*N+1)';
% Assign the wavespeed c
c=nlincoeffs(1,3*M*N+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% We do all plotting here %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the number of contours for plotting
numcontours=20;

% The plotting grid points and meshes
eta = 0:2*pi/1440:2*pi;
phi = [0:29*pi/(180*60):29*pi/180 pi/6:40*pi/(180*120):70*pi/180 71*pi/180:19*pi/(180*60):pi/2];
etamat=zeros(max(size(phi)),max(size(eta)));
phimat=zeros(max(size(phi)),max(size(eta)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the height field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we can rotate the solution by a certain amount in radians
rotangle = pi/kappa;

% Call the c++ function to calculate the free surface height at all grid points
hfield=fsfield(eta,phi,H,kappa,M,N);

% Call the c++ function to calculate the free surface height at all grid points with rotation
%hfield=fsfieldrot(eta,phi,H,kappa,M,N,rotangle);

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
[conts,hands]=contour(xgrid,ygrid,hfield,numcontours);
%clabel(conts,hands,'manual');
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
[conts,hands]=contour(xgrid,ygrid,hfield,numcontours);
%COLORBAR('horiz')
axis off
title(tit)
clabel(conts,hands,'manual');

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
    plot(xcirc,ycirc,'k:','LineWidth',0.25)
end

% phi=pi/3
phival=pi/4;
for i=1:res
    xcirc(1,i)=cos(theta(1,i))*cos(phival)/(1+sin(phival));
    ycirc(i,1)=sin(theta(1,i))*cos(phival)/(1+sin(phival));
end
plot(xcirc,ycirc,'k--','LineWidth',0.25)

% phi=0;
phival=0;
for i=1:res
    xcirc(1,i)=cos(theta(1,i))*cos(phival)/(1+sin(phival));
    ycirc(i,1)=sin(theta(1,i))*cos(phival)/(1+sin(phival));
end
plot(xcirc,ycirc,'k')

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
    plot(xcirc,ycirc,'k:','LineWidth',0.25)
    plot(xcirc,-ycirc,'k:','LineWidth',0.25)
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

%return

% Call the c++ function to calculate the velocity components at all grid points
% Back up old xgrid and ygrid for later use
xgridold=xgrid;
ygridold=ygrid;

% Redefine the mesh resolution for better plots...
eta = 0:2*pi/100:2*pi;
phi = [0:29*pi/(180*30):29*pi/180 pi/6:40*pi/(180*50):70*pi/180 70*pi/180:19*pi/(180*30):pi/2];
etamat=zeros(max(size(phi)),max(size(eta)));
phimat=zeros(max(size(phi)),max(size(eta)));

% Call the c++ function to calculate the velocity components at all grid points
ulamfield=ulam(eta,phi,P,kappa,w,M,N);
uphifield=uphi(eta,phi,Q,kappa,M,N);

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
contour(xgridold,ygridold,hfield,numcontours)
%COLORBAR('horiz')
axis off
quiver(xgrid,ygrid,U,V);
t1 = 'Polar stereographic height and horizontal velocity vector field ';
tit=strcat(t1,t2);
title(tit)