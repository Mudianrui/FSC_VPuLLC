function [ds] = trns01_pd(x)
if x<=0.002 || x>=1
    ds=0;
else
    ds = -1/(exp((1-2*x)/(x*(1-x)))+1)^2 * exp((1-2*x)/(x*(1-x))) * (-2*x*(1-x)-(1-2*x)*(1-2*x))/(x^2*(1-x)^2);
end
