% Clear all variables in current workspace
clear all
tic
% Declare all constants as global variables
global N kappa w MATA MATB h0 Sr Fr Ro Ma gamma

% Declare and initialise all constants
N=50; % Truncation level (user defined)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up a vector of w values to calculate c at.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%wstart=5.1*w;
wstart=5.186594338;
%wstart=6.7;
wend=0.0*w;
numwvals=200;
wvec=wstart:(wend-wstart)/numwvals:wend;
% Find the length of wvec
lenw=length(wvec);
% Initialise the cmat storage for all eigenvalue vectors
% and also room for different values of kappa.
cmat=zeros(lenw,3*N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Finished Initialisation     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise the two matrices
MATA=zeros(3*N);
MATB=zeros(3*N);

for kappa=4:4
    % Perform the loop and populate cmat
    for ii=1:lenw
        % Set the zonal angular velocity to the current iterate value
        w = wvec(1,ii);
        % Now find the new value of h0 that will make the new mass equal in size
        % to the mass from the base compressible calculations. We use the in-built
        % root finding procedure to do this.
        [h0,errorval,exitflag,output] = fzero(@masserror,h0,[],w,a,gamma,Fr,Ma,Ro,Mbp);
        
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
        cmat(ii,:,kappa-3)=eigenvalues';
        % Go back and solve for more eigenvalues
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Finished eigenvalue/eigenvector calculation     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We need to find the indices of each of the dominant eigenvalues
% for each of kappa=3,4 and 5. Do this with a search.
% k4mineigindex=1;
% mineigenvalue=abs(cmat(1,1,1));
% for i=2:3*N
%     if (abs(cmat(1,i,1)) < mineigenvalue)
%         mineigenvalue=abs(cmat(1,i,1));
%         k4mineigindex=i;
%     end
% end

% This test is more robust and will always find the correct
% minimum eigenvalue.
% First find a value that has a positive c val at the right
% end of the graph.
k4mineigindex=1;
while cmat(1,k4mineigindex,1) < 0,
    k4mineigindex = k4mineigindex + 1;
end
% Now k4mineigindex is the index of the first encountered
% positive eigenvalue so store this val.
mineigenvalue=cmat(1,k4mineigindex,1);
% Now loop through the rest of the values and discard any that are negative
for i=k4mineigindex+1:3*N
    % first check for values that are negative
    if cmat(1,i,1) < 0
        % do nothing
    % Now check if we have a new minimum value    
    elseif cmat(1,i,1) < mineigenvalue
        mineigenvalue=cmat(1,i,1);
        k4mineigindex=i;
    else
        % Do nothinf\g
    end
end


% Now that we have them all, plot them up...
figure
hold on
grid on
% Plot kappa=4 line so we can set up legends
%pkappa4=plot(wvec,cmat(:,k4mineigindex,1),'b');
plot(wvec,cmat(:,1:end,1),'b');


% Set the axis for better viewing
%axis([wend,wstart,-20,20]);
axis([wend,wstart,-2,13.5]);
%axis([wend,wstart,-80,80]);
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
%title(t2)
ylabel('Wave speed, c');
xlabel('Zonal flow angular speed, w');

