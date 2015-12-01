function df = jacobian(x)			   

global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed Vzon intol a

%/* Computes the jacobian matrix using
%analytical expressions for the entries. This should
%be correct to machine precision.
%*/
df=zeros(numunknowns);
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
        
        %// Calculate the function derivatives
		%// at the current colloction point in 
		%// a row-wise manner.
		df((i-1)*3*M+(j-1)*3+1,:) = dmass(wavespeed,vars,phival,i,j);
        df((i-1)*3*M+(j-1)*3+2,:) = dlammomen(wavespeed,vars,phival,i,j);
		df((i-1)*3*M+(j-1)*3+3,:) = dphimomen(wavespeed,vars,phival,i,j);	
    end
end
% Calculate the Jacobian elements for the volume specification equation
% Set up some memory
yall=zeros(1,numunknowns+1);
% Here we define the jacobian integral tolerance
jacintol=10^(-1);
% First do the terms with m=0
mval=0;
for nval=0:N
    T1=I10(nval,a,Vzon);
    T2=I20(nval,a,Vzon,x(2*M*N+1:numunknowns,1),N);
    T3=-2.0*kappa/Vzon*dblquad(@I3,0,pi/2,0,pi/kappa,jacintol,@adaptlob,x(2*M*N+1:numunknowns,1),M,N,kappa,a,nval,mval);
    yall(1,2*M*N+nval+1)=T1+T2+T3;
    %df(numunknowns,2*M*N+nval+1)=T1+T2+T3;
end
for mval=1:M-1
    for nval=1:N
        T2=I2i(mval,nval,a,Vzon,x(2*M*N+1:numunknowns,1),N);
        T3=-2.0*kappa/Vzon*dblquad(@I3,0,pi/2,0,pi/kappa,jacintol,@adaptlob,x(2*M*N+1:numunknowns,1),M,N,kappa,a,nval,mval);
        yall(1,2*M*N+mval*N+nval+1)=T2+T3;
        %df(numunknowns,2*M*N+mval*N+nval+1)=T2+T3;
    end
end
% Now fill out df without translating the fixed coeff derivative
for z=1:numunknowns+1
    if (z<fixed)
        df(numunknowns,z)=yall(1,z);
    elseif (z>fixed)
        df(numunknowns,z-1)=yall(1,z);
    else
        % Do nothing since we don't want to translate the fixed coeff
    end
end

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

%// Functions that evaluate the derivatives
%// of the equations of motion with respect to a coefficient
%// at a given point.

function y = dmass(wavespeed,vars,phival,i,j)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed Vzon intol a
    
    % Initialise the row vector that will store all the derivatives
    yall=zeros(1,numunknowns+1);
    % The derivatives without the fixed coeff...we remove this later
    y=zeros(1,numunknowns);
    
    % Terms from using the chain rule on the mass equation
    A1=vars(6,1);
	A2=vars(3,1);
	B1=cos(phival)*vars(9,1)-vars(3,1)*sin(phival);
	B3=vars(3,1)*cos(phival);
	C1=vars(4,1)+cos(phival)*vars(8,1)-sin(phival)*vars(2,1);
	C2=vars(1,1)-Sr*wavespeed*cos(phival);
	C3=cos(phival)*vars(2,1);
	D = -1.0*Sr*cos(phival)*vars(6,1);
    
    % First we need to calculate the dervative for the uncommon indices in each series.
    % That is, the indices that are not part of a double series expression.
    % Terms related to m=0
    % Set up derivative with respect to the constant term H00
    yall(1,2*M*N+1)=C1;
    % Now do the rest of the eta-independent terms
    for n=1:N
        yall(1,2*M*N+n+1)=C1*Cphi2(n,j)-C3*2*n*Sphi2(n,j);
    end
    % Terms related to m=M
    for n=1:N
        yall(1,M*N-N+n)=A1*Ceta(M,i)*Cphi1(n,j)-A2*kappa*M*Seta(M,i)*Cphi1(n,j);
        yall(1,2*M*N-N+n)=B1*Seta(M,i)*Sphi2(n,j)+B3*2*n*Seta(M,i)*Cphi2(n,j);
    end
    % Now do the common series components
    for m=1:M-1
        for n=1:N
            if n==1
                % We need to use an alternative method for the H(i,1) derivatives
                yall(1,(m-1)*N+n)=A1*Ceta(m,i)*Cphi1(n,j)-A2*kappa*m*Seta(m,i)*Cphi1(n,j);
                yall(1,M*N+(m-1)*N+n)=B1*Seta(m,i)*Sphi2(n,j)+B3*2*n*Seta(m,i)*Cphi2(n,j);
                yall(1,2*M*N+m*N+n+1)=-C1*Ceta(m,i)*(Cphi2(1,j)+1.0)+C2*kappa*m*Seta(m,i)*(Cphi2(1,j)+1.0)+C3*Ceta(m,i)*2*Sphi2(1,j);
            else
                yall(1,(m-1)*N+n)=A1*Ceta(m,i)*Cphi1(n,j)-A2*kappa*m*Seta(m,i)*Cphi1(n,j);
                yall(1,M*N+(m-1)*N+n)=B1*Seta(m,i)*Sphi2(n,j)+B3*2*n*Seta(m,i)*Cphi2(n,j);
                yall(1,2*M*N+m*N+n+1)=C1*Ceta(m,i)*((-1)^n)*(Cphi2(n,j)+Cphi2(n-1,j))-C2*kappa*m*Seta(m,i)*(-1)^(n)*(Cphi2(n,j)+Cphi2(n-1,j))-C3*Ceta(m,i)*((-1)^n)*(2*n*Sphi2(n,j)+2*(n-1)*Sphi2(n-1,j));
            end
        end
    end
    % Fill out the derivative with respect to the wavespeed c
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
    
    

function y = dlammomen(wavespeed,vars,phival,i,j)
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed Vzon intol a
    
    % Initialise the row vector that will store all the derivatives
    yall=zeros(1,numunknowns+1);
    % The derivatives without the fixed coeff...we remove this later
    y=zeros(1,numunknowns);
    
    % Terms from using the chain rule on the lambda momentum equation
	A1=vars(4,1)-sin(phival)*vars(2,1);
	A2=vars(1,1)-Sr*wavespeed*cos(phival);
	A3=cos(phival)*vars(2,1);
	B1=cos(phival)*vars(7,1)-(cos(phival)/Ro+vars(1,1))*sin(phival);
	C2=1/Fr^2;
	D = -1.0*Sr*cos(phival)*vars(4,1);
    
    % First we need to calculate the dervative for the uncommon indices in each series.
    % That is, the indices that are not part of a double series expression.
    % Terms related to m=0
    % Set up derivative with respect to the constant term H00
    yall(1,2*M*N+1)=0.0;
    % Now do the rest of the eta-independent terms
    for n=1:N
        yall(1,2*M*N+n+1)=0.0;
    end
    % Terms related to m=M
    for n=1:N
        yall(1,M*N-N+n)=A1*Ceta(M,i)*Cphi1(n,j)-A2*kappa*M*Seta(M,i)*Cphi1(n,j)-A3*(2*n-1)*Ceta(M,i)*Sphi1(n,j);
        yall(1,2*M*N-N+n)=B1*Seta(M,i)*Sphi2(n,j);
    end
    % Now do the common series components
    for m=1:M-1
        for n=1:N
            if n==1
                % We need to use an alternative method for the H(i,1) derivatives
                yall(1,(m-1)*N+n)=A1*Ceta(m,i)*Cphi1(n,j)-A2*kappa*m*Seta(m,i)*Cphi1(n,j)-A3*(2*n-1)*Ceta(m,i)*Sphi1(n,j);
                yall(1,M*N+(m-1)*N+n)=B1*Seta(m,i)*Sphi2(n,j);
                yall(1,2*M*N+m*N+n+1)=C2*kappa*m*Seta(m,i)*(Cphi2(1,j)+1.0);
            else
                yall(1,(m-1)*N+n)=A1*Ceta(m,i)*Cphi1(n,j)-A2*kappa*m*Seta(m,i)*Cphi1(n,j)-A3*(2*n-1)*Ceta(m,i)*Sphi1(n,j);
                yall(1,M*N+(m-1)*N+n)=B1*Seta(m,i)*Sphi2(n,j);
                yall(1,2*M*N+m*N+n+1)=-C2*kappa*m*Seta(m,i)*(-1)^(n)*(Cphi2(n,j)+Cphi2(n-1,j));
            end
        end
    end
    % Fill out the derivative with respect to the wavespeed c
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
    global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed Vzon intol a
    
    % Initialise the row vector that will store all the derivatives
    yall=zeros(1,numunknowns+1);
    % The derivatives without the fixed coeff...we remove this later
    y=zeros(1,numunknowns);
    
    % Terms from using the chain rule on the phi momentum equation
    A1=vars(5,1)+(cos(phival)/Ro+2.0*vars(1,1))*sin(phival);
	B1=cos(phival)*vars(8,1);
	B2=vars(1,1)-Sr*wavespeed*cos(phival);
	B3=cos(phival)*vars(2,1);
	C3=cos(phival)/Fr^2;
	D = -1.0*Sr*cos(phival)*vars(5,1);
    
    % First we need to calculate the dervative for the uncommon indices in each series.
    % That is, the indices that are not part of a double series expression.
    % Terms related to m=0
    % Set up derivative with respect to the constant term H00
    yall(1,2*M*N+1)=0.0;
    % Now do the rest of the eta-independent terms
    for n=1:N
        yall(1,2*M*N+n+1)=-C3*2*n*Sphi2(n,j);
    end
    % Terms related to m=M
    for n=1:N
        yall(1,M*N-N+n)=A1*Ceta(M,i)*Cphi1(n,j);
        yall(1,2*M*N-N+n)=B1*Seta(M,i)*Sphi2(n,j)+B2*kappa*M*Ceta(M,i)*Sphi2(n,j)+B3*2*n*Seta(M,i)*Cphi2(n,j);
    end
    % Now do the common series components
    for m=1:M-1
        for n=1:N
            if n==1
                % We need to use an alternative method for the H(i,1) derivatives
                yall(1,(m-1)*N+n)=A1*Ceta(m,i)*Cphi1(n,j);
                yall(1,M*N+(m-1)*N+n)=B1*Seta(m,i)*Sphi2(n,j)+B2*kappa*m*Ceta(m,i)*Sphi2(n,j)+B3*2*n*Seta(m,i)*Cphi2(n,j);
                yall(1,2*M*N+m*N+n+1)=C3*Ceta(m,i)*2*Sphi2(1,j);
            else
                yall(1,(m-1)*N+n)=A1*Ceta(m,i)*Cphi1(n,j);
                yall(1,M*N+(m-1)*N+n)=B1*Seta(m,i)*Sphi2(n,j)+B2*kappa*m*Ceta(m,i)*Sphi2(n,j)+B3*2*n*Seta(m,i)*Cphi2(n,j);
                yall(1,2*M*N+m*N+n+1)=-C3*Ceta(m,i)*((-1)^n)*(2*n*Sphi2(n,j)+2*(n-1)*Sphi2(n-1,j));
            end
        end
    end
    % Fill out the derivative with respect to the wavespeed c
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
