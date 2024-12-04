function T1=smoothTfun1_pd(z)
if z>0
    T1=exp((z-1)/z)/(z^2);
else
    T1=0;
end
