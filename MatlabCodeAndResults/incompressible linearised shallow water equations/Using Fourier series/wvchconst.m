% Clear all variables in current workspace
%clear all

% Declare all constants as global variables
global N kappa w MATA MATB h0 Sr Fr Ro

% Declare and initialise all constants
N=40; % Truncation level (user defined)
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
arad=a/href; % Dimensionless value of Earth's radius relative to href.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up a vector of w values to calculate c at.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%wstart=5.1*w;
wstart=5.186594338;
%wstart=6.7;
wend=0.0*w;
numwvals=100;
wvec=wstart:(wend-wstart)/numwvals:wend;
% Find the length of wvec
lenw=length(wvec);
% Initialise the cmat storage for all eigenvalue vectors
% and also room for different values of kappa.
cmat=zeros(lenw,3*N,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Finished Initialisation     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise the two matrices
MATA=zeros(3*N);
MATB=zeros(3*N);

for kappa=3:4
    % Perform the loop and populate cmat
    for ii=1:lenw
        % Set the zonal angular velocity to the current iterate value
        w = wvec(1,ii);
        % Now we need to get the value of h0 to make all volumes equal
        h0 = calch0(w,Fr,Ro,arad);
        
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
        cmat(ii,:,kappa-2)=eigenvalues';
        % Go back and solve for more eigenvalues
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Finished eigenvalue/eigenvector calculation     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now that we have them all, plot them up...
figure
hold on
% Plot only 1 row so we can set up legends
pkappa3=plot(wvec,cmat(:,1,1),'g');
% Plot the rest of the rows
plot(wvec,cmat(:,:,1),'g');
% Plot only 1 row so we can set up legends
pkappa4=plot(wvec,cmat(:,1,2),'r');
% Plot the rest of the rows
plot(wvec,cmat(:,:,2),'r');
%legend([pkappa3,pkappa4],'kappa=3','kappa=4',2);
grid on
% Put on Rossby-Haurwitz solution
kappa=3;
rh=((kappa*(3+kappa)*wvec*vref/a-2*Omega)/((1+kappa)*(2+kappa)))/cref;
rhkappa3=plot(wvec,rh,'--');
kappa=4;
rh=((kappa*(3+kappa)*wvec*vref/a-2*Omega)/((1+kappa)*(2+kappa)))/cref;
rhkappa4=plot(wvec,rh,'--');
legend([pkappa3,pkappa4,rhkappa3],'kappa=3','kappa=4','R-H solutions',2);
% Set the axis for better viewing
%axis([wend,wstart,-20,20]);
axis([wend,wstart,-5,15]);
% Set up automatic labelling
kappas = num2str(kappa);
Ns = num2str(N);
wstarts = num2str(wstart);
wends = num2str(wend);

t3 = 'Plot of w v''s c''s for w=[';
t1 = strcat(t3,wends,',',wstarts,']');
%t3 = ' with kappa=';
%t2 = strcat(t1,t3,kappas);
t3 = ' and N=';
t2 = strcat(t1,t3,Ns);

%t2 = strcat(t2,t3,Ns);
title(t2)
ylabel('wavespeed c');
xlabel('zonal flow speed w');

