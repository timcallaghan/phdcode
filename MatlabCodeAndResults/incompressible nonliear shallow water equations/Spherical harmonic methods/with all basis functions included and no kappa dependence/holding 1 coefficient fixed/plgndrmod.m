function y = plgndrmod(l,m,x)
% Computes the associated Legrendre polynomials
% P(m,l){x}. Here m and l are integers satisfying
% 0<=m<=l, while x lies in the range -1<=x<=1
% Note that this is in the special form for
% computing the derivative w.r.t. phi...
% i.e. it has been divided through by sqrt(1-x^2)
i=0;
ll=0;

fact=0;
pll=0;
pmm=0;
pmmp1=0;
somx2=0;

%// Compute P(m,m)
pmm=-1.0;
if (m>1)
	somx2=sqrt((1.0-x)*(1.0+x));
	fact=3.0;
	for i=2:m
		pmm = -fact*somx2*pmm;
		fact = fact + 2.0;
    end
end
if (l==m)
	y = pmm;
    return;
else
	%// Compute P(m,m+1)
	pmmp1=x*(2*m+1)*pmm;
	if l==(m+1)
		y = pmmp1;
        return;
	else
		%// Compute P(m,l), l>m+1
		for ll=m+2:l
			pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
			pmm=pmmp1;
			pmmp1=pll;
        end
        y = pll;
        return;
    end
end
