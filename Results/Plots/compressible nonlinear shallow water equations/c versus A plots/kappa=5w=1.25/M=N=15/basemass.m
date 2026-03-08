function y = basemass(gamma,a,Ma,Fr)
% Computes the base mass for the compressible
% shallow water equations with a constant free
% surface height.

B=(gamma-1)/gamma*(Ma/Fr)^2;
f1=B*(a+1);

temp1=4*pi*(gamma-1)*B^((3-2*gamma)/(gamma-1))/(gamma*(2*gamma-1)*(3*gamma-2));

temp2=2*(f1*(gamma-1))^2+2*B*f1*(gamma-1)*gamma*a+gamma*(2*gamma-1)*(B*a)^2;

y = temp1*temp2;