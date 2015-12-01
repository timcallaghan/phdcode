function [fvec,jac]=residjac(x)
global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed fixedval

%/////////////////////
%// The residual function that returns the residual vector
%// for the equations of motion.

% Recombine the fixed coeff into the coeffs vector
tempx=zeros(numunknowns+1,1);
for z=1:numunknowns+1
    if (z<fixed)
        tempx(z,1)=x(z,1);
    elseif (z>fixed)
        tempx(z,1)=x(z-1,1);
    else
        % Insert fixedval at the position of fixed
        tempx(z,1)=fixedval;
    end
end
   

% Initialise the residual vector
fvec = zeros(numunknowns,1);
phival=0;
vars=zeros(9,1);
i=0;
j=0;
wavespeed=tempx(numunknowns+1,1);

for i=1:M
    %// This loop corresponds to eta
	for j=1:M
	    %// This loop corresponds to phi
		%// Get the phi grid point
		phival = phi(j,1);
		%// Set up the field variables vector
		vars(1,1)=ulameval(tempx,i,j,phival);
		vars(2,1)=uphieval(tempx,i,j);
		vars(3,1)=heval(tempx,i,j,phival);
		%// eta derivatives
		vars(4,1)=ulamevaleta(tempx,i,j);
		vars(5,1)=uphievaleta(tempx,i,j);
		vars(6,1)=hevaleta(tempx,i,j);
		%// phi derivatives
		vars(7,1)=ulamevalphi(tempx,i,j,phival);
		vars(8,1)=uphievalphi(tempx,i,j);
		vars(9,1)=hevalphi(tempx,i,j,phival);
        
        %%%%%%%%%% Testing %%%%%%%%
        % Make the file name from the current loop parameters
        %is = num2str(i);
        %js = num2str(j);
        %filetype='.txt';
        %t1=strcat(is,js,filetype);
        % Write the field vars to file for comparison
        % in Mathematica
        %dlmwrite(t1,vars,'\t')
        
		%// Calculate the function values
		%// at the current colloction point.
		fvec((i-1)*3*M+(j-1)*3+1,1) = lammomen(vars,phival,wavespeed);
		fvec((i-1)*3*M+(j-1)*3+2,1) = phimomen(vars,phival,wavespeed);
		fvec((i-1)*3*M+(j-1)*3+3,1) = mass(vars,phival,wavespeed);
    end
end
% Calculate the jacobian matrix
jac=jacobian(tempx);
%// Calculate the amplitude condition
%fvec(numunknowns,1)=ampcond(coeffs); Not needed when holding one fixed



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions to calculate all the field vars and their derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = ulameval(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

    value = 0.0;
    m=0;
    l=0;
    for m=1:M
	    for l=kappa*m:2*M+kappa*m-1
		    if mod(l-m*kappa,2)==0
		    	value = value + x(1+(m-1)*M+(l-kappa*m)/2,1)*P(m,l-kappa*m+1,j)*C(m,i);
          end
        end
    end
    value = value + w*cos(phival);


function value = uphieval(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

    value = 0.0;
    m=0;
    l=0;
    for m=1:M
	    for l=kappa*m:2*M+kappa*m-1
            if mod(l-m*kappa,2) ~= 0
			    value = value + x(1+M*M+(m-1)*M+(l-kappa*m-1)/2,1)*P(m,l-kappa*m+1,j)*S(m,i);
          end
        end
    end


function value = heval(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

    value = 0.0;
	m=0;
    l=0;
	for m=1:M
		for l=kappa*m:2*M+kappa*m-1
			if mod(l-m*kappa,2)==0
				value = value + x(1+2*M*M+(m-1)*M+(l-kappa*m)/2,1)*P(m,l-kappa*m+1,j)*C(m,i);
            end
        end
    end

	value = value + h0+w*Fr^2*(1/Ro+w)/2.0*cos(phival)^2;


%// Functions for evaluating the eta
%// derivaties of the field variables at a given 
%// point in the flow grid.

function value = ulamevaleta(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	m=0;
    l=0;
	for m=1:M
		for l=kappa*m:2*M+kappa*m-1
			if mod(l-m*kappa,2)==0
				value = value + m*x(1+(m-1)*M+(l-kappa*m)/2,1)*P(m,l-kappa*m+1,j)*S(m,i);
            end
        end
    end

	value = -1.0*value;

function value = uphievaleta(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	m=0;
    l=0;
	for m=1:M
		for l=kappa*m:2*M+kappa*m-1
			if mod(l-m*kappa,2)~=0
				value = value + m*x(1+M*M+(m-1)*M+(l-kappa*m-1)/2,1)*P(m,l-kappa*m+1,j)*C(m,i);
            end
        end
    end

function value = hevaleta(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	m=0;
    l=0;
	for m=1:M
		for l=kappa*m:2*M+kappa*m-1
			if mod(l-m*kappa,2)==0
				value = value + m*x(1+2*M*M+(m-1)*M+(l-kappa*m)/2,1)*P(m,l-kappa*m+1,j)*S(m,i);
            end
        end
    end

	value = -1.0*value;

%// Functions for evaluating the phi
%// derivaties of the field variables at a given 
%// point in the flow grid.


function value = ulamevalphi(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	m=0;
    l=0;
	for m=1:M
		for l=kappa*m:2*M+kappa*m-1
			if mod(l-m*kappa,2)==0
				value = value + x(1+(m-1)*M+(l-kappa*m)/2,1)*Pp(m,l-kappa*m+1,j)*C(m,i);
            end
        end
    end

	value = value-1.0*w*sin(phival);

function value = uphievalphi(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	m=0;
    l=0;
	for m=1:M
		for l=kappa*m:2*M+kappa*m-1
			if mod(l-m*kappa,2)~=0
				value = value + x(1+M*M+(m-1)*M+(l-kappa*m-1)/2,1)*Pp(m,l-kappa*m+1,j)*S(m,i);
            end
        end
    end


function value = hevalphi(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	m=0;
    l=0;
	for m=1:M
		for l=kappa*m:2*M+kappa*m-1
			if mod(l-m*kappa,2)==0
				value = value + x(1+2*M*M+(m-1)*M+(l-kappa*m)/2,1)*Pp(m,l-kappa*m+1,j)*C(m,i);
            end
        end
    end
    
	value = value-1.0*w*Fr^2*(1/Ro+w)*cos(phival)*sin(phival);

    
% The amplitude condition
function value = ampcond(x)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	%// Evalute the free surface value for the
	%// current coeff values
	value = 0.0;
	factor = sqrt(1.0/(2.0*pi));
	m=0;
    l=0;
	for m=1:M
		for l=kappa*m:2*M+kappa*m-1
			if mod(l-m*kappa,2)==0
				value = value + x(1+2*M*M+(m-1)*M+(l-kappa*m)/2,1)*Pamp(m,l-kappa*m+1);
            end
        end
    end

	value = value*factor;
	value = value + h0+w*Fr^2*(1/Ro+w)/4.0;
	value = value-Amp;

    
%% The functions to compute the dynamical physics
function value = lammomen(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	value = kappa*(vars(1,1)-Sr*wavespeed*cos(phival))*vars(4,1)+cos(phival)*vars(2,1)*vars(7,1)-(cos(phival)/Ro+vars(1,1))*sin(phival)*vars(2,1)+kappa*vars(6,1)/Fr^2;

function value = phimomen(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	value = kappa*(vars(1,1)-Sr*wavespeed*cos(phival))*vars(5,1)+cos(phival)*vars(2,1)*vars(8,1)+(cos(phival)/Ro+vars(1,1))*sin(phival)*vars(1,1)+cos(phival)*vars(9,1)/Fr^2;

function value = mass(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
	value = 0.0;
	value = kappa*(vars(1,1)-Sr*wavespeed*cos(phival))*vars(6,1)+cos(phival)*vars(2,1)*vars(9,1)+vars(3,1)*(kappa*vars(4,1)+cos(phival)*vars(8,1)-sin(phival)*vars(2,1));
