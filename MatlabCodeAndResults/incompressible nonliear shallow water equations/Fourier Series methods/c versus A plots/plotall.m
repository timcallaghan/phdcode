% Script to plot up the wavespeed v's amplitude
% relationships found from the nonlinear code

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
NumDirects=134;
% The number of gaps in the data set
numgaps=4;
% The current gap number that we are up to
currentgap=0;

% Since we want to store Ampe, Ampp, and Ampt and c we set aside room for this storage now.
% We will store in the order [c;Ampe;Ampp;Ampt] etc
% We need to subtract off the total number of gaps since we don't store data at the points
AmpsandC=zeros(4,NumDirects-StartDirect+1-numgaps);

% The plotting grid points and meshes
eta = 0:2*pi/1440:2*pi;
%eta = 0:2*pi/40:2*pi;
phi = [0:29*pi/(180*60):29*pi/180 pi/6:40*pi/(180*120):70*pi/180 71*pi/180:19*pi/(180*60):pi/2];
%phi = [0:29*pi/(180*10):29*pi/180 pi/6:40*pi/(180*20):70*pi/180 71*pi/180:19*pi/(180*10):pi/2];
eatmat=zeros(max(size(phi)),max(size(eta)));
phimat=zeros(max(size(phi)),max(size(eta)));

tic

% Now we loop through each directory and find the Amplitudes and wavespeed
% values for each data set.
for direct=StartDirect:NumDirects
    direct
    % We need to test if we are on a gap and if so, increment the number of gaps passed
    % so far by 1
    % Note that these values are HARD-CODED....they are relative to each run!
    if (direct == 64 | direct == 82 | direct == 102)
        currentgap = currentgap+1;
    else  
        % First we need to find the base name of our current directory.
        % This depends on the value of direct
        if direct < 10
            basedirect=strcat('step0',num2str(direct));
        else
            basedirect=strcat('step',num2str(direct));
        end
        % Now the we have our base directory we need to make our file
        % command strings to open the files
        oldparams=strcat(basedirect,'/oldparams.txt');
        linparams=strcat(basedirect,'/linparams.txt');
        current=strcat(basedirect,'/current.txt');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Read in all the parameters %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get the truncation levels used
        MN=dlmread(oldparams,'\r');
        % Assign the parameters
        M=MN(1,1); % The nonlinear longitudinal truncation.
        N=MN(2,1); % The nonlinear latitudinal truncation.
        
        % Get the other model parameters (only read first 5 rows)
        params=dlmread(linparams,'\r',[0 0 4 0]);
        % Assign all the parameters
        kappa=params(1,1); % Number of wavelengths around a latitude circle
        wpercent=params(3,1); % Nondimensional Base zonal flow angular velocity (as a percentage of w defined above)   
        % Modify the value of w based on wpercent
        w=wpercent*wbase;
        % Get the new value of h0
        h0=params(4,1);
        
        % Read in the coefficients calculated from the nonlinear model
        nlincoeffs=dlmread(current,' ');
        % Assign the coeffs to their respective vectors
        P=nlincoeffs(1,1:M*N)';
        Q=nlincoeffs(1,M*N+1:2*M*N)';
        H=nlincoeffs(1,2*M*N+1:3*M*N+1)';
        % Assign the wavespeed c
        c=nlincoeffs(1,3*M*N+2);  
        
        % Call the c++ function to calculate the free surface height at all grid points
        hfield=fsfield(eta,phi,H,kappa,M,N);
        
        [rows,cols]=size(hfield);
        % Set up the mesh grids
        for i=1:cols;
            phimat(:,i)=phi';
        end   
        for i=1:rows
            etamat(i,:)=eta;
        end
        % Use the polar stereographic projection
        xgrid=etamat;
        ygrid=phimat;
        for i=1:rows
            for j=1:cols
                xgrid(i,j)=cos(phimat(i,j))*cos(etamat(i,j))/(1+sin(phimat(i,j)));
                ygrid(i,j)=cos(phimat(i,j))*sin(etamat(i,j))/(1+sin(phimat(i,j)));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Find the contour that has been deformed from the
        %%% pi/4 contour of the base zonal flow
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%% Zonal flow level
        hz=h0+w*Fr^2/4*(1/Ro+w);
        % Find the contour matrix for the given contour level
        figure
        axis([-1.1 1.1 -1.1 1.1]);
        axis equal
        hold on
        [Conts,Hands]=contour(xgrid,ygrid,hfield,[hz hz]);
        % Close the figure
        close
        % Strip off the starting info on what the matix looks like...
        % This is ok for now but WILL need altering if the contours split up and form
        % more than 1 continuous curve
        Contnew=Conts(:,2:end);
        % Find the length of Contnew and make a radius vector
        lencontnew=length(Contnew);
        radiusvec=sqrt(Contnew(1,:).^2+Contnew(2,:).^2);
        % Find the maxium and minimum radii
        [Y,I]=max(radiusvec);
        [YY,II]=min(radiusvec);
        % Now I and II hold the index of the maximum/minimum entries in Contnew so we 
        % get the xmax,ymax,xmin and ymin points as follows
        xmax=Contnew(1,I);
        ymax=Contnew(2,I);
        xmin=Contnew(1,II);
        ymin=Contnew(2,II);
        % We now need to convert these to degrees using the reverse of the
        % polar stereographic transform...
        lammax=atan(ymax/xmax);
        lammin=atan(ymin/xmin);
        phimax=pi/2.0-2.0*atan(Y);
        phimin=pi/2.0-2.0*atan(YY);
        % Convert to degrees and store in a vector
        %Amppointvec=[lammax*180/pi phimax*180/pi;lammin*180/pi phimin*180/pi];
        % Also find the amplitudes for later use
        Ampe=(pi/4-phimax)*180/pi;
        Ampp=(phimin-pi/4)*180/pi;
        Ampt=(phimin-phimax)/2.0*180/pi;
        %Ampvec=[Ampe;Ampp;Ampt];
        % Store the results in AmpsandC
        AmpsandC(1,direct-StartDirect+1-currentgap)=c;
        AmpsandC(2,direct-StartDirect+1-currentgap)=Ampe;
        AmpsandC(3,direct-StartDirect+1-currentgap)=Ampp;
        AmpsandC(4,direct-StartDirect+1-currentgap)=Ampt;
        % Go back and complete the run
    end
end

% Save the results...
save AmpsandC AmpsandC

return

% Create vectors to plot the linearized solution with
% For kappa=4
%lincval=(9.5521611943036622e-001)*ones(1,NumDirects-StartDirect+1);
lincval=(3.9519291252918765e-001)*ones(1,NumDirects-StartDirect+1-numgaps);
% For kappa=5
%lincval=(1.579134708860227e+000)*ones(1,NumDirects-StartDirect+1);

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
plot(AmpsandC(2,1:43),AmpsandC(1,1:43),AmpsandC(2,1:43),AmpsandC(1,1:43),'rx');
plot(AmpsandC(2,44:71),AmpsandC(1,44:71),AmpsandC(2,44:71),AmpsandC(1,44:71),'rx');
plot(AmpsandC(2,72:88),AmpsandC(1,72:88),AmpsandC(2,72:88),AmpsandC(1,72:88),'rx');
plot(AmpsandC(2,89:107),AmpsandC(1,89:107),AmpsandC(2,89:107),AmpsandC(1,89:107),'rx');
plot(AmpsandC(2,108:130),AmpsandC(1,108:130),AmpsandC(2,108:130),AmpsandC(1,108:130),'rx');
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
plot(AmpsandC(2,1:43),AmpsandC(1,1:43));
plot(AmpsandC(2,44:71),AmpsandC(1,44:71));
plot(AmpsandC(2,72:88),AmpsandC(1,72:88));
plot(AmpsandC(2,89:107),AmpsandC(1,89:107));
plot(AmpsandC(2,108:130),AmpsandC(1,108:130));
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
plot(AmpsandC(3,1:43),AmpsandC(1,1:43),AmpsandC(3,1:43),AmpsandC(1,1:43),'rx');
plot(AmpsandC(3,44:71),AmpsandC(1,44:71),AmpsandC(3,44:71),AmpsandC(1,44:71),'rx');
plot(AmpsandC(3,72:88),AmpsandC(1,72:88),AmpsandC(3,72:88),AmpsandC(1,72:88),'rx');
plot(AmpsandC(3,89:107),AmpsandC(1,89:107),AmpsandC(3,89:107),AmpsandC(1,89:107),'rx');
plot(AmpsandC(3,108:130),AmpsandC(1,108:130),AmpsandC(3,108:130),AmpsandC(1,108:130),'rx');
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
plot(AmpsandC(3,1:43),AmpsandC(1,1:43));
plot(AmpsandC(3,44:71),AmpsandC(1,44:71));
plot(AmpsandC(3,72:88),AmpsandC(1,72:88));
plot(AmpsandC(3,89:107),AmpsandC(1,89:107));
plot(AmpsandC(3,108:130),AmpsandC(1,108:130));
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
plot(AmpsandC(4,1:43),AmpsandC(1,1:43),AmpsandC(4,1:43),AmpsandC(1,1:43),'rx');
plot(AmpsandC(4,44:71),AmpsandC(1,44:71),AmpsandC(4,44:71),AmpsandC(1,44:71),'rx');
plot(AmpsandC(4,72:88),AmpsandC(1,72:88),AmpsandC(4,72:88),AmpsandC(1,72:88),'rx');
plot(AmpsandC(4,89:107),AmpsandC(1,89:107),AmpsandC(4,89:107),AmpsandC(1,89:107),'rx');
plot(AmpsandC(4,108:130),AmpsandC(1,108:130),AmpsandC(4,108:130),AmpsandC(1,108:130),'rx');
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
plot(AmpsandC(4,1:43),AmpsandC(1,1:43));
plot(AmpsandC(4,44:71),AmpsandC(1,44:71));
plot(AmpsandC(4,72:88),AmpsandC(1,72:88));
plot(AmpsandC(4,89:107),AmpsandC(1,89:107));
plot(AmpsandC(4,108:130),AmpsandC(1,108:130));
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
plot(AmpsandC(2,1:43),AmpsandC(1,1:43),AmpsandC(2,1:43),AmpsandC(1,1:43),'rx');
plot(AmpsandC(2,44:71),AmpsandC(1,44:71),AmpsandC(2,44:71),AmpsandC(1,44:71),'rx');
plot(AmpsandC(2,72:88),AmpsandC(1,72:88),AmpsandC(2,72:88),AmpsandC(1,72:88),'rx');
plot(AmpsandC(2,89:107),AmpsandC(1,89:107),AmpsandC(2,89:107),AmpsandC(1,89:107),'rx');
plot(AmpsandC(2,108:130),AmpsandC(1,108:130),AmpsandC(2,108:130),AmpsandC(1,108:130),'rx');

plot(AmpsandC(3,1:43),AmpsandC(1,1:43),AmpsandC(3,1:43),AmpsandC(1,1:43),'rx');
plot(AmpsandC(3,44:71),AmpsandC(1,44:71),AmpsandC(3,44:71),AmpsandC(1,44:71),'rx');
plot(AmpsandC(3,72:88),AmpsandC(1,72:88),AmpsandC(3,72:88),AmpsandC(1,72:88),'rx');
plot(AmpsandC(3,89:107),AmpsandC(1,89:107),AmpsandC(3,89:107),AmpsandC(1,89:107),'rx');
plot(AmpsandC(3,108:130),AmpsandC(1,108:130),AmpsandC(3,108:130),AmpsandC(1,108:130),'rx');

plot(AmpsandC(4,1:43),AmpsandC(1,1:43),AmpsandC(4,1:43),AmpsandC(1,1:43),'rx');
plot(AmpsandC(4,44:71),AmpsandC(1,44:71),AmpsandC(4,44:71),AmpsandC(1,44:71),'rx');
plot(AmpsandC(4,72:88),AmpsandC(1,72:88),AmpsandC(4,72:88),AmpsandC(1,72:88),'rx');
plot(AmpsandC(4,89:107),AmpsandC(1,89:107),AmpsandC(4,89:107),AmpsandC(1,89:107),'rx');
plot(AmpsandC(4,108:130),AmpsandC(1,108:130),AmpsandC(4,108:130),AmpsandC(1,108:130),'rx');
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
plot(AmpsandC(2,1:43),AmpsandC(1,1:43));
plot(AmpsandC(2,44:71),AmpsandC(1,44:71));
plot(AmpsandC(2,72:88),AmpsandC(1,72:88));
plot(AmpsandC(2,89:107),AmpsandC(1,89:107));
plot(AmpsandC(2,108:130),AmpsandC(1,108:130));

plot(AmpsandC(3,1:43),AmpsandC(1,1:43));
plot(AmpsandC(3,44:71),AmpsandC(1,44:71));
plot(AmpsandC(3,72:88),AmpsandC(1,72:88));
plot(AmpsandC(3,89:107),AmpsandC(1,89:107));
plot(AmpsandC(3,108:130),AmpsandC(1,108:130));

plot(AmpsandC(4,1:43),AmpsandC(1,1:43));
plot(AmpsandC(4,44:71),AmpsandC(1,44:71));
plot(AmpsandC(4,72:88),AmpsandC(1,72:88));
plot(AmpsandC(4,89:107),AmpsandC(1,89:107));
plot(AmpsandC(4,108:130),AmpsandC(1,108:130));
% Put in linearised solution
plot(AmpsandC(2,:),lincval,'--');
% Label the figure
title('Wavespeed versus Amplitudes');
xlabel('Amplitude in degrees');
ylabel('Wavespeed');

% % Polar amplitude plot
% figure
% hold on
% grid on
% plot(AmpsandC(3,:),AmpsandC(1,:),AmpsandC(3,:),AmpsandC(1,:),'rx',AmpsandC(3,:),lincval,'--');
% title('Wavespeed versus Polar Amplitude');
% xlabel('Polar Amplitude in degrees');
% ylabel('Wavespeed');
% 
% % Averaged amplitude plot
% figure
% hold on
% grid on
% plot(AmpsandC(4,:),AmpsandC(1,:),AmpsandC(4,:),AmpsandC(1,:),'rx',AmpsandC(4,:),lincval,'--');
% title('Wavespeed versus Combined Amplitude');
% xlabel('Combined Amplitude in degrees');
% ylabel('Wavespeed');
% 
%     
% % Plot them all together
% figure
% hold on
% grid on
% plot(AmpsandC(2,:),AmpsandC(1,:),AmpsandC(2,:),AmpsandC(1,:),'rx',AmpsandC(3,:),AmpsandC(1,:),AmpsandC(3,:),AmpsandC(1,:),'rx',AmpsandC(4,:),AmpsandC(1,:),AmpsandC(4,:),AmpsandC(1,:),'rx',AmpsandC(2,:),lincval,'--');
% title('Wavespeed versus Amplitudes');
% xlabel('Amplitude in degrees');
% ylabel('Wavespeed');

toc