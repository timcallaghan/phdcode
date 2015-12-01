% Script to plot up the nonlinear solutions found
% using the c++ code for the compressible flow

% Clear all constants from memory
clear all

% Declare and initialise all constants
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
wbase=wref*aref/vref; % Dimensionless base zonal flow angular velocity
h0=1; % Dimensionless polar free suface height (h/href=1 is assumed at the poles)
a=aref/href; % Dimensionless value of the Earth's radius...relative to href...

% Now compute the base mass of the system and also the base pseudo mass for later calculation.
Mb=basemass(gamma,a,Ma,Fr);
Mbp=basepmass(gamma,a,Ma,Fr);


% Here we declare the number of Directories that we want to
% load data from.
StartDirect=1;
NumDirects=156;
% The number of gaps in the data set
numgaps=1;
% The current gap number that we are up to
currentgap=0;

% Load AmpsandC
load AmpsandC

% Create vectors to plot the linearized solution.
% For kappa=4
lincval=(9.840823862787393e-001)*ones(1,NumDirects-StartDirect+1-numgaps);

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
plot(AmpsandC(2,1:91),AmpsandC(1,1:91),AmpsandC(2,1:91),AmpsandC(1,1:91),'rx');
plot(AmpsandC(2,92:end),AmpsandC(1,92:end),AmpsandC(2,92:end),AmpsandC(1,92:end),'rx');
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
plot(AmpsandC(2,1:91),AmpsandC(1,1:91));
plot(AmpsandC(2,92:end),AmpsandC(1,92:end));
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
plot(AmpsandC(3,1:91),AmpsandC(1,1:91),AmpsandC(3,1:91),AmpsandC(1,1:91),'rx');
plot(AmpsandC(3,92:end),AmpsandC(1,92:end),AmpsandC(3,92:end),AmpsandC(1,92:end),'rx');
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
plot(AmpsandC(3,1:91),AmpsandC(1,1:91));
plot(AmpsandC(3,92:end),AmpsandC(1,92:end));
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
plot(AmpsandC(4,1:91),AmpsandC(1,1:91),AmpsandC(4,1:91),AmpsandC(1,1:91),'rx');
plot(AmpsandC(4,92:end),AmpsandC(1,92:end),AmpsandC(4,92:end),AmpsandC(1,92:end),'rx');
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
plot(AmpsandC(4,1:91),AmpsandC(1,1:91));
plot(AmpsandC(4,92:end),AmpsandC(1,92:end));
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
plot(AmpsandC(2,1:91),AmpsandC(1,1:91),AmpsandC(2,1:91),AmpsandC(1,1:91),'rx');
plot(AmpsandC(2,92:end),AmpsandC(1,92:end),AmpsandC(2,92:end),AmpsandC(1,92:end),'rx');

plot(AmpsandC(3,1:91),AmpsandC(1,1:91),AmpsandC(3,1:91),AmpsandC(1,1:91),'rx');
plot(AmpsandC(3,92:end),AmpsandC(1,92:end),AmpsandC(3,92:end),AmpsandC(1,92:end),'rx');

plot(AmpsandC(4,1:91),AmpsandC(1,1:91),AmpsandC(4,1:91),AmpsandC(1,1:91),'rx');
plot(AmpsandC(4,92:end),AmpsandC(1,92:end),AmpsandC(4,92:end),AmpsandC(1,92:end),'rx');
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
plot(AmpsandC(2,1:91),AmpsandC(1,1:91));
plot(AmpsandC(2,92:end),AmpsandC(1,92:end));

plot(AmpsandC(3,1:91),AmpsandC(1,1:91));
plot(AmpsandC(3,92:end),AmpsandC(1,92:end));

plot(AmpsandC(4,1:91),AmpsandC(1,1:91));
plot(AmpsandC(4,92:end),AmpsandC(1,92:end));
% Put in linearised solution
plot(AmpsandC(2,:),lincval,'--');
% Label the figure
title('Wavespeed versus Amplitudes');
xlabel('Amplitude in degrees');
ylabel('Wavespeed');

