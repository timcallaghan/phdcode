% Clear all variables in current workspace
clear all
tic
% Declare all constants as global variables
global N kappa w MATA MATB h0 Sr Fr Ro Ma gamma

% Declare and initialise all constants
N=100; % Truncation level (user defined)
kappa=4; % Base wavenumber (user defined)
g=9.80616; % Gravitational acceleration
Omega=2*pi/(24*60*60); % Earth's angular velocity
aref=6.37122*10^6; % Radius of the Earth
gamma=1.4; % Ratio of specific heats for dry air

% Declare characteristic scales
cref=Omega/30; % Characteristic angular velocity scale of the Rossby wave (s^-1)
vref=40; % Characteristic velocity scale(m*s^-1)
wref=7.848*10^(-6); %Characteristic angular velocity scale of the zonal flow (s^-1)
pref=1.01325*10^5; % Characteristic pressure value. (1 atmosphere)
rhoref=1.29; % Characteristic density value.

% Compute the value for href that makes the density and pressure at
% sea level equal to their characteristic values.
href=gamma*pref/((gamma-1)*g*rhoref);

% Define the nondimensional numbers
Sr=aref*cref/vref; % Strouhal number
Fr=vref/sqrt(g*href); % Froude number
Ro=vref/(2*Omega*aref); % Rossby number
Ma=vref/sqrt(pref/rhoref); % Mach number

% Dimensionless zonal flow parameters
w=wref*aref/vref; % Dimensionless base zonal flow angular velocity
h0=1; % Dimensionless polar free suface height (h/href=1 is assumed at the poles)
a=aref/href; % Dimensionless value of the Earth's radius...relative to href...

% Now compute the base mass of the system and also the base pseudo mass for later calculation.
Mb=basemass(gamma,a,Ma,Fr);
Mbp=basepmass(gamma,a,Ma,Fr);

% We really just want the mass of the zonal flow because this is the true
% mass that will get perturbed. The Linearized mass will be incorrect for
% the true nonlinear case since to conserve mass the mean height must drop if there are waves
% and the linerized solution can't do this with its series terms.


% Now set the new value of w to what we want it to be...if required...
% else comment out to default to Williamson's w value (as initialised above)
%w = 0.80*w;


% Now find the new value of h0 that will make the new mass equal in size
% to the mass from the base compressible calculations. We use the in-built
% root finding procedure to do this.
[h0,errorval,exitflag,output] = fzero(@masserror,h0,[],w,a,gamma,Fr,Ma,Ro,Mbp);

% This uses a minimisation process to solve it...
%h0 = fsolve(@masserror,h0,[],w,a,gamma,Fr,Ma,Ro,Mbp);

% Now we should have new values for w and h0 that still have the same mass
% as that for the previous values. Thus we can continue to find the eigenvalues now.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Finished Initialisation     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% First initialise the two matrices
MATA=zeros(3*N);
MATB=zeros(3*N);
% Now calculate the elements of MATA and MATB using our routines
calcmatA;
calcmatB;

% Now solve the generalized eigenvalue problem.
[eigenvectors,eigvalues] = eig(MATA,MATB);

% Set up eigenvalues so that we can sort it and still
% keep track of the original index so we can find its
% corresponding eigenvector. 
eigenvalues=zeros(3*N,2);
for i=1:3*N
    eigenvalues(i,1)=eigvalues(i,i);
    % Store the current eigenvalue's original position
    eigenvalues(i,2)=i;
end

% Before we sort we find the minimum eigenvalue (absolute value)
% This will hold the min eig index
mineigindex=1;
mineigenvalue=abs(eigenvalues(1,1));
for i=2:3*N
    if (abs(eigenvalues(i,1)) < mineigenvalue)
        mineigenvalue=abs(eigenvalues(i,1));
        mineigindex=i;
    end
end
% mineigindex holds the index of the min abs val eigenvalue
% so we use this to assign mineigenvalue and mineigenvector 
mineigenvalue=eigenvalues(mineigindex,1);
mineigenvector=eigenvectors(:,mineigindex);
% Now we sort eigenvalues from smallest to largest
eigenvalues=sortrows(eigenvalues,1);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Finished eigenvalue/eigenvector calculation     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% We do all plotting here %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the number of contours for plotting
numcontours=20;

% Define the wavespeed c
c=mineigenvalue;

P=mineigenvector(1:N,1);
Q=mineigenvector(N+1:2*N,1);
H=mineigenvector(2*N+1:3*N,1);

% Find the base height for the free surface at latitude 45N, longitude 0.
baseheight=0;
for n=1:N
    baseheight=baseheight+H(n,1)*(cos(2*n*pi/4)*(-1)^n-cos(2*(n-1)*pi/4)*(-1)^(n-1));
end

% Some parameters necessary for aligning the amplitude
K=7.848*10^(-6); % Parameter used by Phillips to fix the amplitude
% Find the maximum value of the free surface at latitude pi/4
temp1=a^2/g*(wref*(2*Omega+wref)/4-K^2*2^(-kappa-3)*(4*kappa^2+kappa+3));
temp2=a^2/g*(Omega+wref)*K*(kappa^2+2*kappa+3)/((kappa+1)*(kappa+2)*2^(kappa/2));
temp3=-1*a^2/g*K^2*(kappa+3)/2^(kappa+3);
Ampdim=href+temp1+temp2+temp3; % The dimensional amplitude of the R-H wave at latitude 45N
Amp=Ampdim/href; % The dimensionless amplitude of the R-H wave at latitude 45N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Used to set the amplitude without the Rossby-Haurwitx formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %% Get the height of the zonal flow at latitude pi/4 
%  Ampz=h0+w*Fr^2*(1/Ro+w)/4;
%  %% Modify this by a small amount
%  epsz=0.01;
%  Amp=Ampz*(1+epsz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use Amp to fix the size of epsilon so that the amplitudes for the linearised
% and R-H models are consistent at latitude 45N
epsilon=(Amp-h0-w*Fr^2*(1/Ro+w)/4)/baseheight;


% Write all the parameters to file for use in the nonlinear c++ code
% First clear the file contents by opening and closing without appending (ie writing)
fid=fopen('linparams.txt','wt');
fclose(fid);
% Now open again and append all the constants and variables
fid=fopen('linparams.txt','at');
fprintf(fid,'%d',N);
fprintf(fid,'\r');
fprintf(fid,'%d',kappa);
fprintf(fid,'\r');
fprintf(fid,'%.16e',h0);
fprintf(fid,'\r');
fprintf(fid,'%.16e',w);
fprintf(fid,'\r');
fprintf(fid,'%.16e',Amp);
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
fprintf(fid,'%.16e',pref);
fprintf(fid,'\r');
fprintf(fid,'%.16e',rhoref);
fprintf(fid,'\r');
fprintf(fid,'%.16e',Mb);
fprintf(fid,'\r');
fprintf(fid,'%.16e',g);
fprintf(fid,'\r');
fprintf(fid,'%.16e',Omega);
fprintf(fid,'\r');
fprintf(fid,'%.16e',aref);
fprintf(fid,'\r');
fprintf(fid,'%.16e',a);
fprintf(fid,'\r');
fprintf(fid,'%.16e',gamma);
fprintf(fid,'\r');
fclose(fid);


% Write the coefficients and the wavespeed c to file for use in the nonlinear c++ code
% First clear the file contents by opening and closing without appending (ie writing)
fid=fopen('lincoeffs.txt','wt');
fclose(fid);
% Now open again and append all the constants and variables
fid=fopen('lincoeffs.txt','at');
for i=1:N
    fprintf(fid,'%.16e',epsilon*P(i,1));
    fprintf(fid,'\r');
end
for i=1:N
    fprintf(fid,'%.16e',epsilon*Q(i,1));
    fprintf(fid,'\r');
end
for i=1:N
    fprintf(fid,'%.16e',epsilon*H(i,1));
    fprintf(fid,'\r');
end
fprintf(fid,'%.16e',c);
fprintf(fid,'\r');
fclose(fid);

% Return if we don't want to do any plotting
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start the plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The plotting grid points and meshes
eta = 0:2*pi/150:2*pi;
phi = [0:29*pi/(180*60):29*pi/180 pi/6:40*pi/(180*120):70*pi/180 71*pi/180:19*pi/(180*60):pi/2];
etamat=zeros(max(size(phi)),max(size(eta)));
phimat=zeros(max(size(phi)),max(size(eta)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the free surface height field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the c++ function to calculate the free surface height at all grid points
hfield=fsfield(eta,phi,H,kappa,w,Fr,Ro,epsilon,h0);

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
%% Plot the  sea level pressure field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the c++ function to calculate the sea level pressure at all grid points
pfield=pfield(eta,phi,H,kappa,w,Fr,Ro,epsilon,h0,Ma,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot with grid lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
[conts,hands]=contour(xgrid,ygrid,pfield,numcontours);
%COLORBAR('horiz')
axis off
clabel(conts,hands,'manual');

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
t1 = 'Polar stereographic pressure contours';
tit=strcat(t1,t2);
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
%% Plot the  sea level density field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the c++ function to calculate the sea level pressure at all grid points
rhofield=rhofield(eta,phi,H,kappa,w,Fr,Ro,epsilon,h0,Ma,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot with grid lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
[conts,hands]=contour(xgrid,ygrid,rhofield,numcontours);
%COLORBAR('horiz')
axis off
clabel(conts,hands,'manual');

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
t1 = 'Polar stereographic density contours';
tit=strcat(t1,t2);
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

% Call the c++ function to calculate the velocity components at all grid points
% Back up old xgrid and ygrid for later use
xgridold=xgrid;
ygridold=ygrid;

% Redefine the mesh resolution for better plots...
eta = 0:2*pi/100:2*pi;
phi = [0:29*pi/(180*30):29*pi/180 pi/6:40*pi/(180*50):70*pi/180 70*pi/180:19*pi/(180*30):pi/2];

ulamfield=ulam(eta,phi,P,kappa,w,epsilon);
uphifield=uphi(eta,phi,Q,kappa,epsilon);

etamat=zeros(max(size(phi)),max(size(eta)));
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
% Plot both free surface and velocity field
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

%%%%%%%%%%%%%%%%%%%%%
% Plot both pressure and velocity field
%%%%%%%%%%%%%%%%%%%%%
figure
axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
contour(xgridold,ygridold,pfield,numcontours)
%COLORBAR('horiz')
axis off
quiver(xgrid,ygrid,U,V);
t1 = 'Polar stereographic pressure and horizontal velocity vector field ';
tit=strcat(t1,t2);
title(tit)