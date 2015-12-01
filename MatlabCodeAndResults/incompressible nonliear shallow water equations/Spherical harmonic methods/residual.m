function fvec = residual(coeffs)

global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp

%/////////////////////
%// The residual function that returns the residual vector
%// for the equations of motion.

% Initialise the residual vector
fvec = zeros(numunknowns,1);
phival=0;
vars=zeros(9,1);
i=0;
j=0;
wavespeed=coeffs(numunknowns,1);

for i=1:M
    %// This loop corresponds to eta
	for j=1:M
	    %// This loop corresponds to phi
		%// Get the phi grid point
		phival = phi(j,1);
		%// Set up the field variables vector
		vars(1,1)=ulameval(coeffs,i,j,phival);
		vars(2,1)=uphieval(coeffs,i,j);
		vars(3,1)=heval(coeffs,i,j,phival);
		%// eta derivatives
		vars(4,1)=ulamevaleta(coeffs,i,j);
		vars(5,1)=uphievaleta(coeffs,i,j);
		vars(6,1)=hevaleta(coeffs,i,j);
		%// phi derivatives
		vars(7,1)=ulamevalphi(coeffs,i,j,phival);
		vars(8,1)=uphievalphi(coeffs,i,j);
		vars(9,1)=hevalphi(coeffs,i,j,phival);
		%// Calculate the function values
		%// at the current colloction point.
		fvec((i-1)*3*M+(j-1)*3+1,1) = lammomen(vars,phival,wavespeed);
		fvec((i-1)*3*M+(j-1)*3+2,1) = phimomen(vars,phival,wavespeed);
		fvec((i-1)*3*M+(j-1)*3+3,1) = mass(vars,phival,wavespeed);
    end
end
%// Calculate the amplitude condition
fvec(numunknowns,1)=ampcond(coeffs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions to calculate all the field vars and their derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = ulameval(x,i,j,phival)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp

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
    value = value + w*a*cos(phival);


function value = uphieval(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp

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
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp

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

	value = value + h0+w*a*a*(2.0*Omega+w)/(2.0*g)*cos(phival)^2;


%// Functions for evaluating the eta
%// derivaties of the field variables at a given 
%// point in the flow grid.

function value = ulamevaleta(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
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
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
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
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
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
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
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

	value = value-1.0*w*a*sin(phival);

function value = uphievalphi(x,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
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
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
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
    
	value = value-1.0*w*a*a*(2.0*Omega+w)/g*cos(phival)*sin(phival);

    
% The amplitude condition
function value = ampcond(x)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
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
	value = value + h0+w*a*a*(2.0*Omega+w)/(4.0*g);
	value = value-Amp;

    
%% The functions to compute the dynamical physics
function value = lammomen(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
	value = 0.0;
	value = kappa*(vars(1,1)/a-wavespeed*cos(phival))*vars(4,1)+cos(phival)*vars(2,1)*vars(7,1)/a-(2.0*Omega*sin(phival)*cos(phival)+vars(1,1)*sin(phival)/a)*vars(2,1)+kappa*g*vars(6,1)/a;

function value = phimomen(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
	value = 0.0;
	value = kappa*(vars(1,1)/a-wavespeed*cos(phival))*vars(5,1)+cos(phival)*vars(2,1)*vars(8,1)/a+(2.0*Omega*sin(phival)*cos(phival)+vars(1,1)*sin(phival)/a)*vars(1,1)+g*cos(phival)*vars(9,1)/a;

function value = mass(vars,phival,wavespeed)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
	value = 0.0;
	value = kappa*(vars(1,1)/a-wavespeed*cos(phival))*vars(6,1)+cos(phival)*vars(2,1)*vars(9,1)/a+vars(3,1)/a*(kappa*vars(4,1)+cos(phival)*vars(8,1)-sin(phival)*vars(2,1));
