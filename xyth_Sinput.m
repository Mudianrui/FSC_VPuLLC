function [sys,x0,str,ts] = xyth_Sinput(t,x,u,flag)
switch flag
case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
case 3
    sys=mdlOutputs(t,x,u);
case {1,2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end


function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 46;
sizes.NumInputs      = 16;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [];

global Nn
Nn = 5;

global isClockwise Lr Lx Ly Ll Lt SR Ssum
isClockwise = 1;%1 clockwise, -1 counterclockwise
% xe0 = 0;
% ye0 = 0;
Lr = 0.4;
Lx = 3;
Ly = 2.5;
Ll = 0.4;
Lt = 0.6;
SR = [Lx,pi*Lr/2,Ly,pi*Lr/2,Ll,pi*Lr/2,Lt,pi*Lr/2,Lx-2*Ll-4*Lr];%,pi*r/2,Lt,pi*r/2,Ll,pi*r/2,Ly,pi*r/2
Ssum = sum(SR)+sum(SR(2:8));

global SN SCycles
SN = zeros(Nn,1);
SCycles = zeros(Nn,1);


function sys=mdlOutputs(~,~,u)
t = u(1);
global Nn
x = zeros(Nn,1);
y = zeros(Nn,1);
% th = zeros(nn,1);
ni = 1; x(ni) = u(ni*3-1); y(ni) = u(ni*3); %th(ni) = u(ni*3+1);
ni = 2; x(ni) = u(ni*3-1); y(ni) = u(ni*3); %th(ni) = u(ni*3+1);
ni = 3; x(ni) = u(ni*3-1); y(ni) = u(ni*3); %th(ni) = u(ni*3+1);
ni = 4; x(ni) = u(ni*3-1); y(ni) = u(ni*3); %th(ni) = u(ni*3+1);
ni = 5; x(ni) = u(ni*3-1); y(ni) = u(ni*3); %th(ni) = u(ni*3+1);

global vd dvd
vd = 0.8+0.2*sin(t);
dvd = 0.2*cos(t);

s = zeros(Nn,1);
global SN SCycles
% sys(1)
sys(1)=t;
% sys(2)——sys(10)
ni = 1; sSNSC = LaneCaculate(x(ni),y(ni),SN(ni),SCycles(ni));%sSNSC=[xd,yd,thd,vd,wd,isClockwise*ep,RegionNum,S,SN,SCycles];
srt = sSNSC(1:end-2); SN(ni) = sSNSC(end-1); SCycles(ni) = sSNSC(end);
s(ni) = sSNSC(end-2);
sys=[sys,srt,dvd];
% sys(11)——sys(19)
ni = 2; sSNSC = LaneCaculate(x(ni),y(ni),SN(ni),SCycles(ni));
srt = sSNSC(1:end-2); SN(ni) = sSNSC(end-1); SCycles(ni) = sSNSC(end);
s(ni) = sSNSC(end-2);
sys=[sys,srt,s(ni-1)-s(ni)];
% sys(20)——sys(28)
ni = 3; sSNSC = LaneCaculate(x(ni),y(ni),SN(ni),SCycles(ni));
srt = sSNSC(1:end-2); SN(ni) = sSNSC(end-1); SCycles(ni) = sSNSC(end);
s(ni) = sSNSC(end-2);
sys=[sys,srt,s(ni-1)-s(ni)];
% sys(29)——sys(37)
ni = 4; sSNSC = LaneCaculate(x(ni),y(ni),SN(ni),SCycles(ni));
srt = sSNSC(1:end-2); SN(ni) = sSNSC(end-1); SCycles(ni) = sSNSC(end);
s(ni) = sSNSC(end-2);
sys=[sys,srt,s(ni-1)-s(ni)];
% sys(38)——sys(46)
ni = 5; sSNSC = LaneCaculate(x(ni),y(ni),SN(ni),SCycles(ni));
srt = sSNSC(1:end-2); SN(ni) = sSNSC(end-1); SCycles(ni) = sSNSC(end);
s(ni) = sSNSC(end-2);
sys=[sys,srt,s(ni-1)-s(ni)];
