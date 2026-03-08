function fvec = residual(x)
%/////////////////////
%// The residual function that returns the residual vector
%// for the equations of motion.

global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed Vzon intol a

% Initialise the residual vector
fvec = zeros(numunknowns,1);
phival=0;
vars=zeros(9,1);
wavespeed=x(numunknowns+1,1);

for i=1:M
    %// This loop corresponds to eta
	for j=1:N
	    %// This loop corresponds to phi
		%// Get the phi grid point
		phival = phi(j,1);
		%// Set up the field variables vector
		vars(1,1)=ulameval(x,i,j,phival);
		vars(2,1)=uphieval(x,i,j);
		vars(3,1)=heval(x,i,j,phival);
		%// eta derivatives
		vars(4,1)=ulamevaleta(x,i,j);
		vars(5,1)=uphievaleta(x,i,j);
		vars(6,1)=hevaleta(x,i,j);
		%// phi derivatives
		vars(7,1)=ulamevalphi(x,i,j,phival);
		vars(8,1)=uphievalphi(x,i,j);
		vars(9,1)=hevalphi(x,i,j,phival);
        
%         %%%%%%%%%% Testing %%%%%%%%
%         % Make the file name from the current loop parameters
%         is = num2str(i);
%         js = num2str(j);
%         filetype='.txt';
%         t1=strcat(is,js,filetype);
%         % Write the field vars to file for comparison
%         % in Mathematica
%         dlmwrite(t1,vars,'\t')
        
		%// Calculate the function values
		%// at the current colloction point.
        fvec((i-1)*3*M+(j-1)*3+1,1) = mass(vars,phival,wavespeed);
        fvec((i-1)*3*M+(j-1)*3+2,1) = lammomen(vars,phival,wavespeed);
		fvec((i-1)*3*M+(j-1)*3+3,1) = phimomen(vars,phival,wavespeed);
    end
end
% Calculate the current volume for the free surface
Vnl=2.0*kappa*a^2*dblquad(@intgrand,0,pi/2,0,pi/kappa,intol,@adaptlob,x(2*M*N+1:numunknowns,1),M,N,kappa,a);
% Now put in the mass specification condition
fvec(numunknowns,1)=1.0-Vnl/Vzon;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions to calculate all the field vars and their derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = ulameval(x,i,j,phival)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise to the zonal flow
    value = w*cos(phival);
    % Add in the rest of the series
    for m=1:M
	    for n=1:N
		    value = value + x((m-1)*N+n,1)*Ceta(m,i)*Cphi1(n,j);
        end
    end
    
function value = uphieval(x,i,j)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise to the zonal flow
    value = 0.0;
    % Add in the rest of the series
    for m=1:M
	    for n=1:N
            value = value + x(M*N+(m-1)*N+n,1)*Seta(m,i)*Sphi2(n,j);
        end
    end
    
function value = heval(x,i,j,phival)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise to zero.
    value = 0.0;
    % Add in eta-independent series components
    % First add in constant term
    value = value + x(2*M*N+1,1);
    % Now add in the rest of the terms
    for n=1:N
        value = value + x(2*M*N+n+1,1)*Cphi2(n,j);
    end
    % Add in the case when n=1
    for m=1:M-1
        value = value - x(2*M*N+m*N+2,1)*Ceta(m,i)*(Cphi2(1,j)+1.0);
    end
    % Add in the rest of the series
    for m=1:M-1
	    for n=2:N
		    value = value + x(2*M*N+m*N+n+1,1)*Ceta(m,i)*((-1)^n)*(Cphi2(n,j)+Cphi2(n-1,j));
        end
    end

%// Functions for evaluating the eta
%// derivaties of the field variables at a given 
%// point in the flow grid.

function value = ulamevaleta(x,i,j)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise the sum
    value = 0.0;
    % Add in the rest of the series
    for m=1:M
	    for n=1:N
		    value = value - kappa*m*x((m-1)*N+n,1)*Seta(m,i)*Cphi1(n,j);
        end
    end

function value = uphievaleta(x,i,j)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise the sum
    value = 0.0;
    % Add in the rest of the series
    for m=1:M
	    for n=1:N
            value = value + kappa*m*x(M*N+(m-1)*N+n,1)*Ceta(m,i)*Sphi2(n,j);
        end
    end

function value = hevaleta(x,i,j)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise the sum
    value = 0.0;
    % Add in the case when n=1
    for m=1:M-1
        value = value + kappa*m*x(2*M*N+m*N+2,1)*Seta(m,i)*(Cphi2(1,j)+1.0);
    end
    % Add in the rest of the series
    for m=1:M-1
	    for n=2:N
		    value = value - kappa*m*x(2*M*N+m*N+n+1,1)*Seta(m,i)*(-1)^(n)*(Cphi2(n,j)+Cphi2(n-1,j));
        end
    end

%// Functions for evaluating the phi
%// derivaties of the field variables at a given 
%// point in the flow grid.


function value = ulamevalphi(x,i,j,phival)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise to the zonal flow
    value = -w*sin(phival);
    % Add in the rest of the series
    for m=1:M
	    for n=1:N
		    value = value - (2*n-1)*x((m-1)*N+n,1)*Ceta(m,i)*Sphi1(n,j);
        end
    end

function value = uphievalphi(x,i,j)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise to the zonal flow
    value = 0.0;
    % Add in the rest of the series
    for m=1:M
	    for n=1:N
            value = value + 2*n*x(M*N+(m-1)*N+n,1)*Seta(m,i)*Cphi2(n,j);
        end
    end


function value = hevalphi(x,i,j,phival)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
    % Initialise to the zonal flow
    value = 0.0;
    % Add in eta-independent series components
    % Now add in the rest of the terms
    for n=1:N
        value = value - 2*n*x(2*M*N+n+1,1)*Sphi2(n,j);
    end
    % Add in the case when n=1
    for m=1:M-1
        value = value + x(2*M*N+m*N+2,1)*Ceta(m,i)*(2*Sphi2(1,j)+0.0);
    end
    % Add in the rest of the series
    for m=1:M-1
	    for n=2:N
		    value = value - x(2*M*N+m*N+n+1,1)*Ceta(m,i)*((-1)^n)*(2*n*Sphi2(n,j)+2*(n-1)*Sphi2(n-1,j));
        end
    end    

%% The functions to compute the dynamical physics
function value = mass(vars,phival,wavespeed)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
	value = 0.0;
	value = (vars(1,1)-Sr*wavespeed*cos(phival))*vars(6,1)+cos(phival)*vars(2,1)*vars(9,1)+vars(3,1)*(vars(4,1)+cos(phival)*vars(8,1)-sin(phival)*vars(2,1));

function value = lammomen(vars,phival,wavespeed)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
	value = 0.0;
	value = (vars(1,1)-Sr*wavespeed*cos(phival))*vars(4,1)+cos(phival)*vars(2,1)*vars(7,1)-(cos(phival)/Ro+vars(1,1))*sin(phival)*vars(2,1)+vars(6,1)/Fr^2;

function value = phimomen(vars,phival,wavespeed)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed
	value = 0.0;
	value = (vars(1,1)-Sr*wavespeed*cos(phival))*vars(5,1)+cos(phival)*vars(2,1)*vars(8,1)+(cos(phival)/Ro+vars(1,1))*sin(phival)*vars(1,1)+cos(phival)*vars(9,1)/Fr^2;

