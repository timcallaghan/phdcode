function y = plgndr(l,m,x)
% Computes the associated Legrendre polynomials
% P(m,l){x}. Here m and l are integers satisfying
% 0<=m<=l, while x lies in the range -1<=x<=1
% We include the normalising term as well.
i=0;
ll=0;	
fact=0;
pll=0;
pmm=0;
pmmp1=0;
somx2=0;

normalize = sqrt((2.0*l+1.0)/2.0*factrl(l-m)/factrl(l+m));

%// Compute P(m,m)
pmm=1.0;
if (m>0)
	somx2=sqrt((1.0-x)*(1.0+x));
	fact=1.0;
	for i=1:m
		pmm = -fact*somx2*pmm;
		fact = fact+2.0;
    end
end
if (l==m)
	y = normalize*pmm;
    return;
else
	%// Compute P(m,m+1)
	pmmp1=x*(2*m+1)*pmm;
	if l==(m+1)
		y = normalize*pmmp1;
        return;
    else
		%// Compute P(m,l), l>m+1
		for ll=m+2:l
			pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
			pmm=pmmp1;
			pmmp1=pll;
        end
        y = normalize*pll;
        return;
    end
end