function df = fdjac(x,fvec)
% Computes forward-difference approximation to Jacobian. On
% input, x[1..n] is the point at which the Jacobian is to be
% evaluated, fvec[1..n] is the vector of function values at
% the point, and residual(x) is a user-supplied routine that
% returns the vector of functions at x. On output, df[1..n][1..n]
% is the Jaconian array.

% Global variables
global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro fixed

% Approximate square root of the machine precision.
EPS=1.0e-5; 

h=0;
temp=0;
n=numunknowns;
f=zeros(n,1);
fd=zeros(n,n);
for j=1:n+1
    if j < fixed
        %'made it in here'
        temp=x(j,1);
	    h=EPS*abs(temp);
	    if (h == 0.0) 
            h=EPS;
        end
	    % Trick to reduce finite precision error.
	    x(j,1)=temp+h;
	    h=x(j,1)-temp;
        f=residual(x);
	    x(j,1)=temp;
        % Calculate the FD jacobian column-wise
        df(:,j)=(f-fvec)/h;
    elseif j > fixed
        %'made it to here'
        temp=x(j,1);
	    h=EPS*abs(temp);
	    if (h == 0.0) 
            h=EPS;
        end
	    % Trick to reduce finite precision error.
	    x(j,1)=temp+h;
	    h=x(j,1)-temp;
        f=residual(x);
	    x(j,1)=temp;
        % Calculate the FD jacobian column-wise
        df(:,j-1)=(f-fvec)/h;
    else
        %'got here too'
        %do nothing since we are currently on the coeff that is fixed
    end
end

