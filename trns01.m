function [s] = trns01(x)
if x<=0
    s=0;
elseif x>=1
    s=1;
else
    s = 1/(exp((1-2*x)/(x*(1-x)))+1);
end
