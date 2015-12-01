function df = jacobian(x)			   

global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
%/* Computes the jacobian matrix using
%analytical expressions for the entries. This should
%be correct to machine precision.
%*/
phival=0;
wavespeed=x(numunknowns+1,1);
vars=zeros(9,1);
i=0;
j=0;

for i=1:M
    %// This loop corresponds to eta
	for j=1:M
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
        
        %// Calculate the function derivatives
		%// at the current colloction point in 
		%// a row-wise manner.
		df((i-1)*3*M+(j-1)*3+1,:) = dlammomen(wavespeed,vars,phival,i,j);
		df((i-1)*3*M+(j-1)*3+2,:) = dphimomen(wavespeed,vars,phival,i,j);
		df((i-1)*3*M+(j-1)*3+3,:) = dmass(wavespeed,vars,phival,i,j);
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

%// Functions that evaluate the derivatives
%// of the equations of motion with respect to a coefficient
%// at a given point.

function y = dlammomen(wavespeed,vars,phival,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
    
    % Initialise the row vector that will store all the derivatives
    yall=zeros(1,numunknowns+1);
    % The derivatives without the fixed coeff
    y=zeros(1,numunknowns);
    
    % Terms from using the chain rule on the lambda momentum equation
	A1=kappa*vars(4,1)-sin(phival)*vars(2,1);
	A2=kappa*(vars(1,1)-Sr*wavespeed*cos(phival));
	A3=cos(phival)*vars(2,1);
	B1=cos(phival)*vars(7,1)-(cos(phival)/Ro+vars(1,1))*sin(phival);
	C2=kappa/Fr^2;
	D = -1.0*kappa*Sr*cos(phival)*vars(4,1);
       
	for row=1:M
		for col=1:M
			yall(1,M*(row-1)+col)=A1*P(row,2*col-1,j)*C(row,i)-A2*row*P(row,2*col-1,j)*S(row,i)+A3*Pp(row,2*col-1,j)*C(row,i);
			yall(1,M*M+M*(row-1)+col)=B1*P(row,2*col,j)*S(row,i);
			yall(1,2*M*M+M*(row-1)+col)= -1.0*C2*row*P(row,2*col-1,j)*S(row,i);
        end
    end
	%// Fill out the derivative with respect to the wavespeed c
	yall(1,numunknowns+1)=D;
    
    % Now fill out y without translating the fixed coeff derivative
    for z=1:numunknowns+1
        if (z<fixed)
            y(1,z)=yall(1,z);
        elseif (z>fixed)
            y(1,z-1)=yall(1,z);
        else
             % Do nothing since we don't want to translate the fixed coeff
        end
    end

function y = dphimomen(wavespeed,vars,phival,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
    
    % Initialise the row vector that will store all the derivatives
    yall=zeros(1,numunknowns+1);
    % The derivatives without the fixed coeff
    y=zeros(1,numunknowns);
    
    % Terms from using the chain rule on the phi momentum equation
    A1=kappa*vars(5,1)+(cos(phival)/Ro+2.0*vars(1,1))*sin(phival);
	B1=cos(phival)*vars(8,1);
	B2=kappa*(vars(1,1)-Sr*wavespeed*cos(phival));
	B3=cos(phival)*vars(2,1);
	C3=cos(phival)/Fr^2;
	D = -1.0*kappa*Sr*cos(phival)*vars(5,1);
       
	for row=1:M
		for col=1:M
			yall(1,M*(row-1)+col)=A1*P(row,2*col-1,j)*C(row,i);
			yall(1,M*M+M*(row-1)+col)=B1*P(row,2*col,j)*S(row,i)+B2*row*P(row,2*col,j)*C(row,i)+B3*Pp(row,2*col,j)*S(row,i);
			yall(1,2*M*M+M*(row-1)+col)= C3*Pp(row,2*col-1,j)*C(row,i);
        end
    end
	%// Fill out the derivative with respect to the wavespeed c
	yall(1,numunknowns+1)=D;
    
    % Now fill out y without translating the fixed coeff derivative
    for z=1:numunknowns+1
        if (z<fixed)
            y(1,z)=yall(1,z);
        elseif (z>fixed)
            y(1,z-1)=yall(1,z);
        else
             % Do nothing since we don't want to translate the fixed coeff
        end
    end

function y = dmass(wavespeed,vars,phival,i,j)
    global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed
    
    % Initialise the row vector that will store all the derivatives
    yall=zeros(1,numunknowns+1);
    % The derivatives without the fixed coeff
    y=zeros(1,numunknowns);
    
    % Terms from using the chain rule on the mass equation
    A1=kappa*vars(6,1);
	A2=vars(3,1)*kappa;
	B1=cos(phival)*vars(9,1)-vars(3,1)*sin(phival);
	B3=vars(3,1)*cos(phival);
	C1=kappa*vars(4,1)+cos(phival)*vars(8,1)-sin(phival)*vars(2,1);
	C2=kappa*(vars(1,1)-Sr*wavespeed*cos(phival));
	C3= cos(phival)*vars(2,1);
	D = -1.0*kappa*Sr*cos(phival)*vars(6,1);
       
	for row=1:M
		for col=1:M
			yall(1,M*(row-1)+col)=A1*P(row,2*col-1,j)*C(row,i)-A2*row*P(row,2*col-1,j)*S(row,i);
			yall(1,M*M+M*(row-1)+col)=B1*P(row,2*col,j)*S(row,i)+B3*Pp(row,2*col,j)*S(row,i);
			yall(1,2*M*M+M*(row-1)+col)= C1*P(row,2*col-1,j)*C(row,i)-C2*row*P(row,2*col-1,j)*S(row,i)+C3*Pp(row,2*col-1,j)*C(row,i);
        end
    end
	%// Fill out the derivative with respect to the wavespeed c
	yall(1,numunknowns+1)=D;
    
    % Now fill out y without translating the fixed coeff derivative
    for z=1:numunknowns+1
        if (z<fixed)
            y(1,z)=yall(1,z);
        elseif (z>fixed)
            y(1,z-1)=yall(1,z);
        else
             % Do nothing since we don't want to translate the fixed coeff
        end
    end

