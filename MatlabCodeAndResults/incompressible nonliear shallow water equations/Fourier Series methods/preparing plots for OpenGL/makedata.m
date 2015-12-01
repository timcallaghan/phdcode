% Script to plot up the nonlinear solutions found
% using the c++ code...it uses an automation approach
% to generate sequential *.png files for use in the
% OpenGL viewer program...

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
NumDirects=9;
% The number of gaps in the data set
numgaps=4;
% The current gap number that we are up to
currentgap=0;

% The plotting grid points
eta = 0:2*pi/1440:2*pi;
phi = [0:29*pi/(180*60):29*pi/180 pi/6:40*pi/(180*120):70*pi/180 71*pi/180:19*pi/(180*60):pi/2];
%eta = 0:2*pi/150:2*pi;
%phi = [0:29*pi/(180*25):29*pi/180 pi/6:40*pi/(180*30):70*pi/180 71*pi/180:19*pi/(180*25):pi/2];

% Set the number of contours for plotting
numcontours=20;

% The latitude and longitude mesh grids
% We are making a mercator plot here...
[Xmesh,Ymesh]=meshgrid(eta,phi);



% Now we loop through each directory and generate a plot if it's
% not a gap point...
for direct=StartDirect:NumDirects
    direct
    % We need to test if we are on a gap and if so, increment the number of gaps passed
    % so far by 1
    % Note that these values are HARD-CODED....they are relative to each run!
    if (direct == 44 | direct == 73 | direct == 91 | direct == 111)
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% We do all plotting here %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Call the c++ function to calculate the free surface height at all grid points
        hfield=fsfield(eta,phi,H,kappa,M,N);
        
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
        
        
        % Print the current figure to eps. Since we are not storing these pictures
        % we can overwrite them each time...
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

        % Read in our new png image
        'Reading in PNG'
        PNG=imread('picture4.png');
        % Create an alpha channel for the PNG
        'Creating alpha channel'
        Achan=makealpha(PNG,255);
        
        % Write the image to file with the alpha channel included. We 
        % need to construct the filename from
        % the current directory location...this is given by
        % direct-StartDirect+1-currentgap
        pngname = strcat(num2str(direct-StartDirect+1-currentgap),'.png');
        'Writing Alpha PNG'
        imwrite(PNG,pngname,'png','Alpha',Achan,'BitDepth',8);
        % Go back and complete the rest of the pictures...
    end
end

% Now change directories and do the other plotting
%cd('C:\matlabR12\work\compressible nonlinear shallow water equations\c versus A plots\kappa=4w=1.25\M=N=15')
% call the plotting program to generate c v's A plots
%plotall