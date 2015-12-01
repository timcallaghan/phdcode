function yy = gammln(xx)

	
%// Returns the value of ln[gamma(xx)] for xx > 0;
j=0;
x=0;
y=0;
tmp=0;
ser=0;
	
cof = [76.18009172947146; -86.50532032941677 ;24.01409824083091;-1.231739572450155;0.128650973866179e-2;-0.5395239384953e-5];

y=xx;
x=xx;
tmp=x+5.5;
tmp = tmp -(x+0.5)*log(tmp);
ser=1.000000000190015;
for j=0:5
    y = y+1;
	ser = ser + cof(j+1,1)/y;
end

yy = -tmp+log(2.5066282746310005*ser/x);