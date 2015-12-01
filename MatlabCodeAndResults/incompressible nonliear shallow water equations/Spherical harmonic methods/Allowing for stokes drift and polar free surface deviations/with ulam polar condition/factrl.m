function y = factrl(n)

%// Returns the value of n! as a floating point number.
ntop=4;
%// Fill in table only as required...
a=zeros(33,1);
a(1:5,1)=[1.0;1.0;2.0;6.0;24.0];
j=0;
if (n > 32)
    %// Larger value than size of table is required so use
	%// the gamma function to evaluate it. 
	%// Note it may overflow
	y = exp(gammln(n+1.0));
    return;
end

while (ntop<n)
    %// Fill in table up to desired value.
	ntop = ntop+1;
    j=ntop;
	a(ntop+1,1)=a(j,1)*ntop;
end

y = a(n+1,1);