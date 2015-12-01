% Script to plot up the wavespeed v's amplitude
% relationships found from the nonlinear code

% This script uses a previous value of AmpsandC to plot so no real calculation here...

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

% Here we declare the number of Directories that we want to
% load data from.
StartDirect=1;
NumDirects=203;
% The number of gaps in the data set
numgaps=1;
% The current gap number that we are up to
currentgap=0;

% Since we want to store Ampe, Ampp, and Ampt and c we set aside room for this storage now.
% We will store in the order [c;Ampe;Ampp;Ampt] etc
% We need to subtract off the total number of gaps since we don't store data at the points
%AmpsandC=zeros(4,NumDirects-StartDirect+1-numgaps);
% Load AmpsandC
load AmpsandC


% Create vectors to plot the linearized solution.
% For kappa=5
lincval=(1.579134708860227e+000)*ones(1,NumDirects-StartDirect+1-numgaps);

% Now that we have all the data we do some plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Equatorial amplitude plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do plot with both line and 'X's at plot points
% Create a new held figure with grid
figure
hold on
grid on
% Plot each solution segment
plot(AmpsandC(2,1:101),AmpsandC(1,1:101),AmpsandC(2,1:101),AmpsandC(1,1:101),'rx');
plot(AmpsandC(2,102:end),AmpsandC(1,102:end),AmpsandC(2,102:end),AmpsandC(1,102:end),'rx');
% Put in linearised solution
plot(AmpsandC(2,:),lincval,'--');
% Label the figure
title('Wavespeed versus Equatorial Amplitude');
xlabel('Equatorial Amplitude in degress');
ylabel('Wavespeed');

% Now do plot with line only.
% Create a new held figure with grid
figure
hold on
grid on
% Plot each solution segment
plot(AmpsandC(2,1:101),AmpsandC(1,1:101));
plot(AmpsandC(2,102:end),AmpsandC(1,102:end));
% Put in linearised solution
plot(AmpsandC(2,:),lincval,'--');
% Label the figure
title('Wavespeed versus Equatorial Amplitude');
xlabel('Equatorial Amplitude in degress');
ylabel('Wavespeed');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Polar amplitude plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do plot with both line and 'X's at plot points
% Create a new held figure with grid
figure
hold on
grid on
% Plot each solution segment
plot(AmpsandC(3,1:101),AmpsandC(1,1:101),AmpsandC(3,1:101),AmpsandC(1,1:101),'rx');
plot(AmpsandC(3,102:end),AmpsandC(1,102:end),AmpsandC(3,102:end),AmpsandC(1,102:end),'rx');
% Put in linearised solution
plot(AmpsandC(3,:),lincval,'--');
% Label the figure
title('Wavespeed versus Polar Amplitude');
xlabel('Polar Amplitude in degrees');
ylabel('Wavespeed');

% Now do plot with line only.
% Create a new held figure with grid
figure
hold on
grid on
% Plot each solution segment
plot(AmpsandC(3,1:101),AmpsandC(1,1:101));
plot(AmpsandC(3,102:end),AmpsandC(1,102:end));
% Put in linearised solution
plot(AmpsandC(3,:),lincval,'--');
% Label the figure
title('Wavespeed versus Polar Amplitude');
xlabel('Polar Amplitude in degrees');
ylabel('Wavespeed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Averaged amplitude plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do plot with both line and 'X's at plot points
% Create a new held figure with grid
figure
hold on
grid on
% Plot each solution segment
plot(AmpsandC(4,1:101),AmpsandC(1,1:101),AmpsandC(4,1:101),AmpsandC(1,1:101),'rx');
plot(AmpsandC(4,102:end),AmpsandC(1,102:end),AmpsandC(4,102:end),AmpsandC(1,102:end),'rx');
% Put in linearised solution
plot(AmpsandC(4,:),lincval,'--');
% Label the figure
title('Wavespeed versus Combined Amplitude');
xlabel('Averaged Amplitude in degrees');
ylabel('Wavespeed');

% Now do plot with line only.
% Create a new held figure with grid
figure
hold on
grid on
% Plot each solution segment
plot(AmpsandC(4,1:101),AmpsandC(1,1:101));
plot(AmpsandC(4,102:end),AmpsandC(1,102:end));
% Put in linearised solution
plot(AmpsandC(4,:),lincval,'--');
% Label the figure
title('Wavespeed versus Combined Amplitude');
xlabel('Averaged Amplitude in degrees');
ylabel('Wavespeed');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot them all together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do plot with both line and 'X's at plot points
% Create a new held figure with grid
figure
hold on
grid on
% Plot each solution segment
plot(AmpsandC(2,1:101),AmpsandC(1,1:101),AmpsandC(2,1:101),AmpsandC(1,1:101),'rx');
plot(AmpsandC(2,102:end),AmpsandC(1,102:end),AmpsandC(2,102:end),AmpsandC(1,102:end),'rx');

plot(AmpsandC(3,1:101),AmpsandC(1,1:101),AmpsandC(3,1:101),AmpsandC(1,1:101),'rx');
plot(AmpsandC(3,102:end),AmpsandC(1,102:end),AmpsandC(3,102:end),AmpsandC(1,102:end),'rx');

plot(AmpsandC(4,1:101),AmpsandC(1,1:101),AmpsandC(4,1:101),AmpsandC(1,1:101),'rx');
plot(AmpsandC(4,102:end),AmpsandC(1,102:end),AmpsandC(4,102:end),AmpsandC(1,102:end),'rx');
% Put in linearised solution
plot(AmpsandC(2,:),lincval,'--');
% Label the figure
title('Wavespeed versus Amplitudes');
xlabel('Amplitude in degrees');
ylabel('Wavespeed');

% Now do plot with line only.
% Create a new held figure with grid
figure
hold on
grid on
% Plot each solution segment
plot(AmpsandC(2,1:101),AmpsandC(1,1:101));
plot(AmpsandC(2,102:end),AmpsandC(1,102:end));

plot(AmpsandC(3,1:101),AmpsandC(1,1:101));
plot(AmpsandC(3,102:end),AmpsandC(1,102:end));

plot(AmpsandC(4,1:101),AmpsandC(1,1:101));
plot(AmpsandC(4,102:end),AmpsandC(1,102:end));
% Put in linearised solution
plot(AmpsandC(2,:),lincval,'--');
% Label the figure
title('Wavespeed versus Amplitudes');
xlabel('Amplitude in degrees');
ylabel('Wavespeed');

