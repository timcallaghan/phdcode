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
P=nlincoeffs(1,1:M*N)';
Q=nlincoeffs(1,M*N+1:2*M*N)';
H=nlincoeffs(1,2*M*N+1:3*M*N+1)';
% Assign the wavespeed c
c=nlincoeffs(1,3*M*N+2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% We do all plotting here %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the number of contours for plotting
numcontours=20;

% The plotting grid points
eta = 0:2*pi/150:2*pi;
phi = [0:29*pi/(180*25):29*pi/180 pi/6:40*pi/(180*30):70*pi/180 71*pi/180:19*pi/(180*25):pi/2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the pressure field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the c++ function to calculate the free surface height at all grid points
hfield=fsfield(eta,phi,H,kappa,M,N);

% Make a mercator plot...
[Xmesh,Ymesh]=meshgrid(eta,phi);
figure
%axis([-1.1 1.1 -1.1 1.1]);
axis equal
hold on
contour(Xmesh,Ymesh,hfield,numcontours)
% Flip the ycoords and map the lower half of the picture
Ymesh=-Ymesh;
[Cont,Conthand]=contour(Xmesh,Ymesh,hfield,numcontours);
% Now we add in a box around the edge so we can easily crop it later on...
phimax=pi/2;
phimin=-pi/2;
etamin=0;
etamax=2*pi;
etavec=[etamin etamax];
phivec=[phimin phimax];
top=[phimax phimax];
bottom=[phimin phimin];
left=[etamin etamin];
right=[etamax etamax];
plot(left,phivec,'w');
plot(etavec,top,'w');
plot(right,phivec,'w');
plot(etavec,bottom,'w');
% Keep background color
%set(gcf, 'InvertHardCopy', 'off');
%axis([etamin etamax phimin phimax]);
axis off
axis square
axis tight
% Print the current figure to eps...we need to construct the file name from
% the parameters in the model really...ok for now as a test...
print -depsc picture4
% Close the figure
close
% Use the perl script fixbbox to remove the excess bounding box area
% Note that this script is specific to these graphics and wont work 
% for general exported graphics
'Fixing the bounding box'
!fixbbox.pl picture4.eps

% Scale the figure to the size require for our OpenGL textures
% This particular scaling makes 1024x1024 pixel images since original cropped eps
% is always 350x350...therefore scale=1024/350...etc
%!scale.pl 2.925714286 newpicture4.eps picture4.eps
%!scale.pl 1.462857143 newpicture4.eps picture4.eps
%!scale.pl 0.731428571 newpicture4.eps picture4.eps
% This makes 2048x2048 textures
'Scaling'
!scale.pl 5.851428571 newpicture4.eps picture4.eps

% Now use imagemagick to convert our picture to png format
'Converting to PNG'
!convert +antialias picture4.eps picture4.png

% For the sake of testing we make a small image first
%!scale.pl 1.462857143 newpicture4.eps smallpic4.eps
%!convert +antialias smallpic4.eps smallpic4.png

% Read in our new png image
'Reading in PNG'
PNG=imread('picture4.png');
% Create an alpha channel for the PNG
'Creating alpha channel'
Achan=makealpha(PNG,255);
% Write the image to file with the alpha channel included
'Writing Alpha PNG'
imwrite(PNG,'test.png','png','Alpha',Achan,'BitDepth',8);
% Convert the png to a TGA file
%'Converting to TGA'
%!convert test.png test.tga
