function y = f1(h0,w,Fr,Ro,a)
%F1 Summary of this function goes here
%  Detailed explanation goes here

y = 2*pi/3*(2*h0^3+2*h0^2*w*Fr^2*(1/Ro+w)+4/5*h0*w^2*Fr^4*(1/Ro+w)^2+4/35*w^3*Fr^6*(1/Ro+w)^3);

y = y + 2*pi*a^2*(2*h0+2/3*w*Fr^2*(1/Ro+w));

y = y + 2*pi*a*(2*h0^2+4/3*h0*w*Fr^2*(1/Ro+w)+4/15*w^2*Fr^4*(1/Ro+w)^2);