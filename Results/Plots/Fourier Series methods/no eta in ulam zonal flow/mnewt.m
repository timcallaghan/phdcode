function [y,errf,errx] = mnewt(ntrial,x,tolx,tolf)
% Given and initial guess x[1..n] for a root in n dimensions, take
% ntrial Newton-Raphson steps to improve the root. Stop if the root
% converges in either summed absolute variabe increments tolx or summed
% absolute function values tolf.

% Global variables
global M N kappa w h0 phi Ceta Seta Cphi1 Sphi1 Cphi2 Sphi2 numunknowns Sr Fr Ro fixed Vzon intol a maxhalvings

% Initialise the coeff vector
y=x;
errx=0;
errf=0;
n=numunknowns;
p=zeros(n,1);
delx=zeros(n,1);
delxnew=zeros(n+1,1);
fvec=zeros(n,1);
fjac=zeros(n,n);
for k=1:ntrial
	% Get the function values at y for
    % fvec and the Jacobian matrix.
    fvec=residual(y);
	%// Check function convergence
	errf=sum(abs(fvec))
	if errf <= tolf
        'Tolf reached'
        return;
    end
    % Use analytical jacobian
    fjac=jacobian(y);
    % Use Finite difference jacobian
    %fjac=fdjac(y,fvec);
    %'Condition number='
    %cond(fjac)

	%// Right hand side of linear equations.
	p = -fvec;
	% Solve the linear equations.
    %delx=fjac\p;
    % Can use a QR factorisation here instead.
    [Q,R]=qr(fjac);
    ytemp = Q'*p;
    delx = R\ytemp;
    errx=norm(delx,1)
    % Translate this to a vector with all coeffs
    for i=1:n+1
        if (i<fixed)
            delxnew(i,1)=delx(i,1);
        elseif (i>fixed)
            delxnew(i,1)=delx(i-1,1);
        else
             % Insert a zero at the position of fixed
            delxnew(i,1)=0;
        end
    end
    % Check for infinite coeffs
    if (sum(abs(delxnew)) == Inf)
        'Infinte values have been calculated...returning'
        return;
    end
    % Implement the halving technique
    % First calulate ytemp
    ytemp= y + delxnew;
    % Now calculate the residuals for this updated y value
    fvec=residual(ytemp);
    % and the one-norm error
    errytemp=sum(abs(fvec));
    % Initialise counter
    halvits=0;
    % Check to see if we need to use halving and if so do it
    while ((errytemp > errf) & halvits <=maxhalvings)
        delxnew=delxnew/2;
        ytemp= y + delxnew;
        fvec=residual(ytemp);
        errytemp=sum(abs(fvec));
        halvits = halvits +1;
    end
    % If maximum halving iterations was reached then we need to exit with an error
    if (halvits==maxhalvings+1)
        'maximum halvings reached...start from new guess'
        return;
    end
    % Print out halvings
    halvits
    
    % Update the vector y
    y = y + delxnew;
    % Check for convergence
	errx=sum(abs(delxnew))
    '*************'
	if errx <= tolx
        'Tolx reached'
        return;
    end
end