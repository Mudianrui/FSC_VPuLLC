function [sys,x0,str,ts] = S0_controller(t,x,u,flag)
%% Scontroller
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
%% mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 35;
sizes.NumInputs      = 19;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [];
global IFe0_0 e0_0
IFe0_0 = 0;
e0_0 = 0;
global ve0
ve0 = 0;

sn = 5;
global IFes_0 es_0
IFes_0 = zeros(sn,1);
es_0 = zeros(sn,1);
global IFdes_0 des_0
IFdes_0 = zeros(sn,1);
des_0 = zeros(sn,1);
global se0
se0 = zeros(sn,1);
global IFeth_0 eth_0
IFeth_0 = zeros(sn,1);
eth_0 = zeros(sn,1);
global rhoth0
rhoth0 = zeros(sn,1);

function sys=mdlOutputs(~,~,u)
si = 1+0;
%% input
t = u(1);
xr = u(2);
yr = u(3);
thr = u(4);
vr = u(5);
wr = u(6);
ep = u(7);
Rn = u(8);
S = u(9);
dvr = u(10);
x = u(11);
y = u(12);
th = u(13);
v = u(14);
w = u(15);
Izetap0 = u(16);
hatvartheta0 = u(17);
Izetapth = u(18);
hatvarthetath = u(19);
the = thNormalization(th-thr);
rn = wr/vr;% = 1/Lr;
ds = v*cos(the)/(1-ep*rn);

%% controller
T = 5;
Lw = 0.05;%Lane constraint: half of lane width
lambda = 1;
alpha = 0.6;
beta = 3;
Gamma = trns01((T-t)/T);
dGamma = -trns01_pd((T-t)/T)/T;

%% ACC
epsilon0 = 0.1;
k0alpha = 1;
k0beta = 1;
sigma0 = 1;
gamma0 = 1;
kappa0alpha = 1;
kappa0beta = 1;
e0 = v-vr;
global IFe0_0 e0_0
if IFe0_0==0
    IFe0_0 = 1;
    e0_0 = e0;
end
e0W = e0-e0_0*Gamma;
global ve0
if ve0==0
    ve0 = abs(e0)+epsilon0;
end
rho0_0 = 2;%ve0*2;
rho0 = (rho0_0-epsilon0)*smoothTfun1((T-t)/T)+epsilon0;
drho0 = -(rho0_0-epsilon0)*smoothTfun1_pd((T-t)/T)/T;
zeta0 = constraintTF(e0W,rho0);
zetap0 = xppow(zeta0);
[mu,upsilon] = constraintTF_pd(e0W,rho0);
mu0 = mu;
upsilon0 = -mu*dGamma*e0_0+upsilon*drho0;
S0 = zeta0+lambda*Izetap0;
phi0 = sum([abs(v),abs(dvr),abs((upsilon0+lambda*zetap0)/mu0)])+2;
uv = -(k0alpha*sgn(S0,alpha)+k0beta*sgn(S0,beta))/mu0-sigma0*S0*mu0*(hatvartheta0*phi0)^2;
dhatvartheta0 = gamma0*abs(S0*mu0)*phi0-kappa0alpha*hatvartheta0^alpha-kappa0beta*hatvartheta0^beta;
% uv = -10*zeta0;

%% LK
epsilonth = 0.2*pi/4;
kthalpha = 0.4;
kthbeta = kthalpha;
sigmath = 1;
gammath = 1;
kappathalpha = kthalpha;
kappathbeta = kthalpha;
kpi = 1;%kpi<=1
Dw = 0.3;
thd = thr+atan(tan(kpi*pi/4)*ep/Lw);
ethB = th-thd;
eth = thNormalization(ethB);
global IFeth_0 eth_0
if IFeth_0(si)==0
    IFeth_0(si) = 1;
    eth_0(si) = eth;
end
ethW = eth;
global rhoth0
if rhoth0(si)==0
    rhoth0(si) = abs(eth)+epsilonth;
end
rhoth_0 = pi/4;
rhoth = (rhoth_0-epsilonth)*smoothTfun1((T-t)/T)+epsilonth;
drhoth = -(rhoth_0-epsilonth)*smoothTfun1_pd((T-t)/T)/T;
zetath = constraintTF(ethW,rhoth)+Dw*w;
zetapth = xppow(zetath);
[mu,upsilon] = constraintTF_pd(ethW,rhoth);
muth = mu;
upsilonth = upsilon*drhoth;
Ath = w-ds*rn+tan(kpi*pi/4)/Lw/(1+(tan(kpi*pi/4)*ep/Lw)^2)*v*sin(the);
Bth = Dw;
Sth = zetath+lambda*Izetapth;
phith = sum([abs(muth*Ath),abs(upsilonth+lambda*zetapth),Bth*abs(w),Bth]);
uw = -(kthalpha*sgn(Sth,alpha)+kthbeta*sgn(Sth,beta))/Bth-sigmath*Sth*(hatvarthetath*phith)^2/Bth;
dhatvarthetath = gammath*abs(Sth)*phith-kappathalpha*hatvarthetath^alpha-kappathbeta*hatvarthetath^beta;

%% output
D = [1 1;1 -1];
url = D\[uv;uw];
% sys(1)——sys(3)
sys = [url(1),url(2),t];
% sys(4)——sys(11)
sys = [sys,xr,yr,thr,vr,wr,ep,Rn,S];
if Rn==1
    th = thNormalization(th+pi)-pi;
end
% sys(12)——sys(16)
sys = [sys,x,y,th,v,w];
% sys(17)——sys(28)
sys = [sys,-Lw,Lw,e0,-rho0,rho0,thd,eth,-rhoth,rhoth,0,0,0];
% sys(29)——sys(30)
sys = [sys,uv,uw];
% sys(31)——sys(34)
sys = [sys,zetap0,dhatvartheta0,zetapth,dhatvarthetath];
% sys(35)
sys = [sys,ds];
