function y = dlegendre(n,phival)
%This program is a MODIFICATION of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%       ===============================================
%       Purpose: Compute Legendre polynomial derivatives Pn'(x)
%       Input :  x --- Argument of Pn(x)
%                n --- Degree of Pn(x) ( n = 0,1,...)
%       Output:  PD(n) --- Pd(x)
%       ===============================================

%// Convert from phival to x
x = sin(phival);

% Define the value of Pn(x) for n=0
p0=1.0d0;
% Define the value for n=1
p1=x;

% declare pf
pf =[];

% If n is 0 or 1 we can return immediately.
if n==0
    y = 0.0d0;
elseif n==1
    y = 1.0*cos(phival);
else
    for  k=2:n;
        if k==n
            % Do nothing...calculate this step once the loop has terminated...see below
        else
            pf=(2.0d0.*k-1.0d0)./k.*x.*p1-(k-1.0d0)./k.*p0;
            p0=p1;
            p1=pf;
        end
    end;
    pf=(2.0d0.*n-1.0d0)./n.*x.*p1-(n-1.0d0)./n.*p0;
    if (abs(x) == 1.0d0);
        y=cos(phival)*0.5d0.*x.^(n+1).*n.*(n+1.0d0);
    else;
        y=cos(phival)*n.*(p1-x.*pf)./(1.0d0-x.*x);
    end;
end

