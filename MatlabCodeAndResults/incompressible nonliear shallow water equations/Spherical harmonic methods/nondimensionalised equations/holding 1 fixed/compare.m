function [A,B,C] = compare

% Declare all constants as global variables

global N kappa g a Omega w tol MAT scale h0

% Declare and initialise all constants
%N=100; % Truncation level (user defined)
kappa=4; % Base wavenumber (user defined)
g=9.80616; % Gravitational acceleration
Omega=2*pi/(24*60*60); % Earth's angular velocity
a=6.37122*10^6; % Radius of the Earth
h0=8*10^3; % Polar free suface height
K=7.848*10^(-6); % Parameter used by Phillips to fix the amplitude
% Zonal flow parameters
w=7.848*10^(-6);%1/10*Omega; % Base zonal flow angular velocity (user defined)

% Set the number of contours for plotting
numcontours=20;

% Find the maximum value of the free surface at latitude pi/4
temp1=a^2/g*(w*(2*Omega+w)/4-K^2*2^(-kappa-3)*(4*kappa^2+kappa+3));
temp2=a^2/g*(Omega+w)*K*(kappa^2+2*kappa+3)/((kappa+1)*(kappa+2)*2^(kappa/2));
temp3=-1*a^2/g*K^2*(kappa+3)/2^(kappa+3);
Amp=h0+temp1+temp2+temp3; % The amplitude of the R-H wave at latitude 45N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Used to set the amplitude without the Rossby-Haurwitx formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the height of the zonal flow at latitude pi/4
%Ampz=h0+w*a^2*(2*Omega+w)/(4*g);
%% Modify this by a small amount
%epsz=0.000001;
%Amp=Ampz*(1+epsz);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Scale factor for nice determinant sizes
scale=a/1000;

% Quadrature error tolerance
tol=[];%10^(-10);

% Initialise the static elements of MAT by loading from file
load MAT;
% Find the truncation used to generate MAT
Nmat=length(MAT)/3;
% Using current value of w, assemble the matrix MAT,
% noting that B32 is stored in the place of C3
% and B33 is stored in the place of B2. Make sure we do the dependent
% entries first i.e. A3 and B3
%
% Do A3
MAT(2*Nmat+1:3*Nmat,1:Nmat) = -scale*kappa*(w*a*(2*Omega+w)/(2.0*g)*MAT(2*Nmat+1:3*Nmat,1:Nmat)+h0/a*MAT(1:Nmat,2*Nmat+1:3*Nmat));
% Do B3
MAT(2*Nmat+1:3*Nmat,Nmat+1:2*Nmat) = scale*((-3.0*w*a*(2*Omega+w)/(2.0*g))*MAT(2*Nmat+1:3*Nmat,Nmat+1:2*Nmat)+(w*a*(2*Omega+w)/(2.0*g))*MAT(2*Nmat+1:3*Nmat,2*Nmat+1:3*Nmat)+h0/a*(MAT(Nmat+1:2*Nmat,2*Nmat+1:3*Nmat)-MAT(Nmat+1:2*Nmat,Nmat+1:2*Nmat)));
% Do B1
MAT(1:Nmat,Nmat+1:2*Nmat) = -2*scale*(Omega+w)*MAT(1:Nmat,Nmat+1:2*Nmat);
% Do C1
MAT(1:Nmat,2*Nmat+1:3*Nmat) = -scale*kappa*g/a*MAT(1:Nmat,2*Nmat+1:3*Nmat);
% Do A2
MAT(Nmat+1:2*Nmat,1:Nmat) = -1*MAT(1:Nmat,Nmat+1:2*Nmat);
% Do C2
MAT(Nmat+1:2*Nmat,2*Nmat+1:3*Nmat) = scale*g/a*MAT(Nmat+1:2*Nmat,2*Nmat+1:3*Nmat);

% Set the new value for N (N<=Nmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE...ONLY USE EVEN VALUES FOR N SO THAT EACH SET OF COEFFS
% HAS THE SAME NUMBER OF COEFFS IN IT. Makes it easier for exporting them...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=10;
% check if MAT needs resizing
if N < Nmat
    MAT=resizeMAT(Nmat);
elseif N > Nmat
    Nmats = num2str(Nmat);
    strcat('N must be less than Nmat= ', Nmats,', setting N=Nmat')
    N=Nmat;
else
    % do nothing since N=Nmat...
end


%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

% Define the initial guess for the wavespeed
guess=1.8*10^(-6);

% Calculate the eigenvalue
eigval = eigsolve(guess)

%define the wavespeed c
c=eigval;

% Write all the parameters to file for use in the nonlinear c++ code
% First clear the file contents by opening and closing without appending (ie writing)
fid=fopen('linparams.txt','wt');
fclose(fid);
% Now open again and append all the constants and variables
fid=fopen('linparams.txt','at');
fprintf(fid,'%d',N/2);
fprintf(fid,'\r');
fprintf(fid,'%d',kappa);
fprintf(fid,'\r');
fprintf(fid,'%.16e',h0);
fprintf(fid,'\r');
fprintf(fid,'%.16e',w);
fprintf(fid,'\r');
fprintf(fid,'%.16e',Amp);
fprintf(fid,'\r');
fclose(fid);


%Calculate the matrix with detX(c)=0
X=calcmat(c);

%split off the right hand side (first column of lhs)
rhs=-1.0*MAT(2:end,1);

%resize the matrix
mat=MAT(2:end,2:end);

%solve for the remaining coeffs
y=mat\rhs;

scalfac=max(abs(y));

A=[1;y(1:N-1)];
A=A/scalfac;
B=y(N:2*N-1);
B=B/scalfac;
C=y(2*N:3*N-1);
C=C/scalfac;


% Find the maximum value of uphi at latitude 45N for the unscaled uphi velocity component
%maxuphival = abs(fminbnd('minuphi',0,2*pi,[],pi/4,B,kappa));%minuphi(eta,pi/4,B,kappa);
% Use Amp to fix the size of epsilon so that the amplitudes are consistent at latitude 45N
epsilon=(Amp-h0-w*a^2*(2*Omega+w)/(4*g))/h1(pi/4,C,kappa)


% Write the coefficients and the wavespeed c to file for use in the nonlinear c++ code
% First clear the file contents by opening and closing without appending (ie writing)
fid=fopen('lincoeffs.txt','wt');
fclose(fid);
% Now open again and append all the constants and variables
fid=fopen('lincoeffs.txt','at');
for i=1:N/2
    fprintf(fid,'%.16e',sqrt(2*pi)*epsilon*A(2*i-1));
    fprintf(fid,'\r');
end
for i=1:N/2
    fprintf(fid,'%.16e',sqrt(2*pi)*epsilon*B(2*i));
    fprintf(fid,'\r');
end
for i=1:N/2
    fprintf(fid,'%.16e',sqrt(2*pi)*epsilon*C(2*i-1));
    fprintf(fid,'\r');
end
fprintf(fid,'%.16e',c);
fprintf(fid,'\r');
fclose(fid);

eta = 0:2*pi/150:2*pi;
phi = [0:29*pi/(180*5):29*pi/180 pi/6:40*pi/(180*30):70*pi/180 71*pi/180:19*pi/(180*5):pi/2];




% base zonal flow max height
%maxh=w*a^2*(2*Omega+w)/(2*g)+h0;
% Declare the small amplitude purturbation parameter in terms of the max height
%epsilon=0.035*maxh;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start the plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the pressure field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the c++ function to calculate the free surface height at all grid points
hfield=fsfield(eta,phi,C,kappa,w,a,g,Omega,epsilon,h0);

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
ulamfield=ulam(eta,phi,A,kappa,w,a,epsilon);
uphifield=uphi(eta,phi,B,kappa,epsilon);

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