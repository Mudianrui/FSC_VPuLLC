function [sys,x0,str,ts] = S1_plant_vw(t,x,u,flag)
switch flag
case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1
    sys=mdlDerivatives(t,x,u);
case 3
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 2;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [0,0];
str = [];
ts  = [0 0];

function sys=mdlDerivatives(t,x,u)
v = x(1);
w = x(2);
ur = u(1);
ul = u(2);
url = [ur;ul];
D = [1 1;1 -1];
us = D*url;
P = 1/180*pi*sin(t);
a1 = -4.5;
b1 = 3.3;
c1 = 0.1;
d1 = 0.1*sin(t);
a2 = -2.8;
b2 = 12.5;
d2 = 0.1*cos(t);
uv = us(1);
uw = us(2);
dv = a1*v + b1*uv + c1*sin(P) + d1;
dw = a2*w + b2*uw + d2;
sys(1) = dv;
sys(2) = dw;

function sys=mdlOutputs(~,x,~)
sys(1) = x(1);
sys(2) = x(2);
