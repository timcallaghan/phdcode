% Clear all variables in current workspace
%clear all

% Declare all constants as global variables
global N kappa w MATA MATB h0 Sr Fr Ro

% Declare and initialise all constants
N=5; % Truncation level (user defined)
kappa=4; % Base wavenumber (user defined)
g=9.80616; % Gravitational acceleration
Omega=2*pi/(24*60*60); % Earth's angular velocity
a=6.37122*10^6; % Radius of the Earth

% Declare characteristic scales
href=8*10^3; % Characteristic length scale (m)
cref=Omega/30; % Characteristic angular velocity scale of the Rossby wave (s^-1)
vref=40; % Characteristic velocity scale(m*s^-1)
wref=7.848*10^(-6); %Characteristic angular velocity scale of the zonal flow (s^-1)

% Define the nondimensional numbers
Sr=a*cref/vref; % Strouhal number
Fr=vref/sqrt(g*href); % Froude number
Ro=vref/(2*Omega*a); % Rossby number

% Dimensionless zonal flow parameters
w=wref*a/vref; % Dimensionless base zonal flow angular velocity
h0=1; % Dimensionless polar free suface height (h/href=1 is assumed at the poles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up a vector of w values to calculate c at.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wstart=1.3*w;
wend=0.0*w;
numwvals=50;
wvec=wstart:(wend-wstart)/numwvals:wend;
% Find the length of wvec
lenw=length(wvec);
% Initialise the cmat storage for 3 eigenvalue vectors
cmat=zeros(lenw,3*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Finished Initialisation     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise the two matrices
MATA=zeros(3*N);
MATB=zeros(3*N);

% Perform the loop and populate cmat
for ii=1:lenw
    % Set the zonal angular velocity to the current iterate value
    w = wvec(1,ii);
    
    % Now calculate the elements of MATA and MATB using our routines
    % and the current value of w.
    calcmatA;
    calcmatB;
    
    % Now solve the generalized eigenvalue problem. We only
    % need the eigenvalues so don't bother calculating the
    % eigenvectors.
    eigenvalues = eig(MATA,MATB);
    
    % Now we sort eigenvalues from smallest to largest
    eigenvalues=sortrows(eigenvalues,1);
    
    % Assign all the eigenvalues to cmat
    cmat(ii,:)=eigenvalues';
    % Go back and solve for more eigenvalues
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Finished eigenvalue/eigenvector calculation     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now that we have them all, plot them up...
figure
plot(wvec,cmat,'r')
grid on
hold on

% Put on Rossby-Haurwitz solution
rh=((kappa*(3+kappa)*wvec*vref/a-2*Omega)/((1+kappa)*(2+kappa)))/cref;
plot(wvec,rh,'--')

% Set up automatic labelling
kappas = num2str(kappa);
Ns = num2str(N);
wstarts = num2str(wstart);
wends = num2str(wend);

t3 = 'Plot of w v''s c''s for w=[';
t1 = strcat(t3,wends,',',wstarts,']');
t3 = ' with kappa=';
t2 = strcat(t1,t3,kappas);
t3 = ', and N=';
t2 = strcat(t2,t3,Ns);
title(t2)
ylabel('wavespeed c');
xlabel('zonal flow speed w');
