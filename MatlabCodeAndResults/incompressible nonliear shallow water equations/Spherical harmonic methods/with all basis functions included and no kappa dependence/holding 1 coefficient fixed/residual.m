function fvec = residual(x)
%/////////////////////////////////////////////////////////
%// The residual function that returns the residual vector
%// for the equations of motion.
%/////////////////////////////////////////////////////////

global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

% Initialise the residual vector
fvec = zeros(numunknowns,1);
% Declare a variable to hold the value of phi
phival=0;
% Allocate space for our field variables
vars=zeros(9,1);
% Get the current value of the wavespeed
wavespeed=x(numunknowns+1,1);

for i=1:2*M
    % This loop corresponds to eta. We have 2*M eta grid points.
	for j=1:M
	    % This loop corresponds to phi. We have M phi grid points.
		% Get the value of phi at the 
        % current collocation grid point
		phival = phi(j,1);
		% Set up the field variables vector
		vars(1,1)=ulameval(x,i,j,phival);
		vars(2,1)=uphieval(x,i,j);
		vars(3,1)=heval(x,i,j,phival);
		% eta derivatives
		vars(4,1)=ulamevaleta(x,i,j);
		vars(5,1)=uphievaleta(x,i,j);
		vars(6,1)=hevaleta(x,i,j);
		% phi derivatives
		vars(7,1)=ulamevalphi(x,i,j,phival);
		vars(8,1)=uphievalphi(x,i,j);
		vars(9,1)=hevalphi(x,i,j,phival);
        
        %%%%%%%%%% Testing %%%%%%%%
        % Make the file name from the current loop parameters
        %is = num2str(i);
        %js = num2str(j);
        %filetype='.txt';
        %t1=strcat(is,js,filetype);
        % Write the field vars to file for comparison
        % in Mathematica
        %dlmwrite(t1,vars,'\t')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
		% Calculate the function values
		% at the current colloction point.
		fvec((i-1)*3*M+(j-1)*3+1,1) = lammomen(vars,phival,wavespeed);
		fvec((i-1)*3*M+(j-1)*3+2,1) = phimomen(vars,phival,wavespeed);
		fvec((i-1)*3*M+(j-1)*3+3,1) = mass(vars,phival,wavespeed);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions to calculate all the field vars and their derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = ulameval(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

    value = 0.0;
    m=0;
    l=0;
    % First add in those basis functions that are using
    % basis recombination
    for l=2:2*M+1
	   if mod(l,2)==0
		   	value = value + x(l/2,1)*P(M+2,l-1,j)*sqrt(1.0/(2.0*pi));
        end
    end
    % Now add on the basis functions that are not
    % using recombination
    % First add in cosine based basis functions
    for m=1:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(1+m*M+(l-m)/2,1)*P(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end
    % Now add in sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(M*M+(m-1)*M+1+(l-m)/2,1)*P(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end 
    % Finally, add in the base zonal flow structure
    value = value + w*cos(phival);


function value = uphieval(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

    value = 0.0;
    m=0;
    l=0;
    % First add in those basis functions that are using
    % basis recombination
    for l=2:2*M+1
	   if mod(l,2)~=0
		   	value = value + x(2*M*M+(l-1)/2,1)*P(M+2,l-1,j)*sqrt(1.0/(2.0*pi));
        end
    end
    % Now add on the basis functions that are not
    % using recombination
    % First add in cosine based basis functions
    for m=1:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)~=0
		    	value = value + x(2*M*M+1+m*M+(l-m-1)/2,1)*P(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end
    % Now add in sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)~=0
		    	value = value + x(3*M*M+(m-1)*M+1+(l-m-1)/2,1)*P(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end     
    
    
function value = heval(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

    value = 0.0;
    m=0;
    l=0;
    % We don't have any basis functions that are using
    % basis recombination so ignore this step...

    % Now add on the basis functions that are not
    % using recombination
    % First add in cosine based basis functions
    for m=0:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(4*M*M+1+m*M+(l-m)/2,1)*P(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end
    % Now add in sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(5*M*M+(m-1)*M+1+(l-m)/2,1)*P(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end 
    % Finally, add in the base zonal flow structure
    value = value + h0+w*Fr^2*(1/Ro+w)/2.0*cos(phival)^2;
    
      
% Functions for evaluating the eta
% derivaties of the field variables at a given 
% point in the flow grid.

function value = ulamevaleta(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
    
    value = 0.0;
    m=0;
    l=0;
    % No contribution from the m=0 case...
    
    % Now add on the basis functions that are not
    % using recombination
    % First add in derivatives of the cosine based basis functions
    for m=1:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(1+m*M+(l-m)/2,1)*P(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end
    % Adjust for the negative sign from the derivative
    value = -1.0*value;
    % Now add in derivatives of the sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(M*M+(m-1)*M+1+(l-m)/2,1)*P(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end 
    % No eta derivative of the base zonal flow structure  
    
function value = uphievaleta(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
    
    value = 0.0;
    m=0;
    l=0;
    % No contribution from the m=0 case...
    
    % Now add on the basis functions that are not
    % using recombination
    % First add in derivatives of the cosine based basis functions
    for m=1:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)~=0
		    	value = value + x(2*M*M+1+m*M+(l-m-1)/2,1)*P(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end
    % Adjust for the negative sign from the derivative
    value = -1.0*value;
    % Now add in derivatives of the sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)~=0
		    	value = value + x(3*M*M+(m-1)*M+1+(l-m-1)/2,1)*P(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end     
      
function value = hevaleta(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

    value = 0.0;
    m=0;
    l=0;
    % We don't have any basis functions that are using
    % basis recombination so ignore this step...

    % Now add on the basis functions that are not
    % using recombination
    % First add in derivatives of the cosine based basis functions
    for m=0:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(4*M*M+1+m*M+(l-m)/2,1)*P(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end
    % Adjust for the negative sign from the derivative
    value = -1.0*value;
    % Now add in derivatives of the sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(5*M*M+(m-1)*M+1+(l-m)/2,1)*P(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end 
    % No eta derivative of the base zonal flow structure 
    
% Functions for evaluating the phi
% derivaties of the field variables at a given 
% point in the flow grid.


function value = ulamevalphi(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
    
    value = 0.0;
    m=0;
    l=0;
    % First add in those basis functions that are using
    % basis recombination
    for l=2:2*M+1
	   if mod(l,2)==0
		   	value = value + x(l/2,1)*Pp(M+2,l-1,j)*sqrt(1.0/(2.0*pi));
        end
    end
    % Now add on the basis functions that are not
    % using recombination
    % First add in phi derivatives of the cosine based basis functions
    for m=1:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(1+m*M+(l-m)/2,1)*Pp(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end
    % Now add in phi derivatives of the sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(M*M+(m-1)*M+1+(l-m)/2,1)*Pp(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end 
    % Finally, add in the phi derivaive of the base zonal flow structure
    value = value-1.0*w*sin(phival);
    
function value = uphievalphi(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
    
    value = 0.0;
    m=0;
    l=0;
    % First add in those basis functions that are using
    % basis recombination
    for l=2:2*M+1
	   if mod(l,2)~=0
		   	value = value + x(2*M*M+(l-1)/2,1)*Pp(M+2,l-1,j)*sqrt(1.0/(2.0*pi));
        end
    end
    % Now add on the basis functions that are not
    % using recombination
    % First add in phi derivatives of the cosine based basis functions
    for m=1:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)~=0
		    	value = value + x(2*M*M+1+m*M+(l-m-1)/2,1)*Pp(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end
    % Now add in phi derivatives of the sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)~=0
		    	value = value + x(3*M*M+(m-1)*M+1+(l-m-1)/2,1)*Pp(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end     
    
function value = hevalphi(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
    
    value = 0.0;
    m=0;
    l=0;
    % We don't have any basis functions that are using
    % basis recombination so ignore this step...

    % Now add on the basis functions that are not
    % using recombination
    % First add in the phi derivative of the cosine based basis functions
    for m=0:M-1
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(4*M*M+1+m*M+(l-m)/2,1)*Pp(m+1,l-m+1,j)*C(m+1,i);
          end
        end
    end
    % Now add in the phi derivative of the sine based basis functions
    for m=1:M
	    for l=m:2*M+m-1
		    if mod(l-m,2)==0
		    	value = value + x(5*M*M+(m-1)*M+1+(l-m)/2,1)*Pp(m+1,l-m+1,j)*S(m+1,i);
          end
        end
    end 
    % Finally, add in the phi derivative of the base zonal flow structure
    value = value-1.0*w*Fr^2*(1/Ro+w)*cos(phival)*sin(phival);
    
%% The functions to compute the dynamical physics
function value = lammomen(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	value = (vars(1,1)-Sr*wavespeed*cos(phival))*vars(4,1)+cos(phival)*vars(2,1)*vars(7,1)-(cos(phival)/Ro+vars(1,1))*sin(phival)*vars(2,1)+vars(6,1)/Fr^2;

function value = phimomen(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	value = (vars(1,1)-Sr*wavespeed*cos(phival))*vars(5,1)+cos(phival)*vars(2,1)*vars(8,1)+(cos(phival)/Ro+vars(1,1))*sin(phival)*vars(1,1)+cos(phival)*vars(9,1)/Fr^2;

function value = mass(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	value = (vars(1,1)-Sr*wavespeed*cos(phival))*vars(6,1)+cos(phival)*vars(2,1)*vars(9,1)+vars(3,1)*(vars(4,1)+cos(phival)*vars(8,1)-sin(phival)*vars(2,1));
