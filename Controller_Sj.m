function sys = Controller_Sj(si,u)
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
Sd = u(10);
x = u(11);
y = u(12);
th = u(13);
v = u(14);
w = u(15);
Izetaps = u(16);
hatvarthetas = u(17);
Izetapth = u(18);
hatvarthetath = u(19);
des_0i = u(20);
vj_1 = u(21);
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
epsilons = 0.1;
ksalpha = 10;
ksbeta = 10;
sigmas = 1;
gammas = 1;
kappasalpha = 1;
kappasbeta = 1;
d = 0.25;
h = 0.25;
wb = 1;
esB = Sd;
sd = d+h*ds;
es = esB-sd;
global IFes_0 es_0
if IFes_0(si)==0
    IFes_0(si) = 1;
    es_0(si) = es;
end
global IFdes_0 des_0
if IFdes_0(si)<3
    IFdes_0(si) = IFdes_0(si)+1;
    des_0(si) = des_0i;
end
esu = (es_0(si)+(wb*es_0(si)+des_0(si))*t)*exp(-wb*t);
esW = es-Gamma*esu;
global se0
if se0(si)==0
    se0(si) = abs(es)+epsilons;
end
rhos_0 = 2;%se0(si)*2;
rhos = (rhos_0-epsilons)*smoothTfun1((T-t)/T)+epsilons;
drhos = -(rhos_0-epsilons)*smoothTfun1_pd((T-t)/T)/T;
zetas = constraintTF(esW,rhos);
zetaps = xppow(zetas);
[mu,upsilon] = constraintTF_pd(esW,rhos);
mus = mu;
upsilons = -mu*(Gamma*(des_0(si)-(wb*es_0(si)+des_0(si))*wb*t)*exp(-wb*t)+dGamma*esu)+upsilon*drhos;
As = vj_1-ds+h*v*sin(the)*((w-ds*rn)*(1-ep*rn)+rn*v*cos(the))/(1-ep*rn)^2;
Bs = h*cos(the)/(1-ep*rn);
Ss = zetas+lambda*Izetaps;
phis = sum([abs(As),Bs*abs(v),Bs,abs((upsilons+lambda*zetaps)/mus)]);
uv = (ksalpha*sgn(Ss,alpha)+ksbeta*sgn(Ss,beta))/Bs/mus+sigmas*Ss*mus*(hatvarthetas*phis)^2/Bs;
dhatvarthetas = gammas*abs(Ss*mus)*phis-kappasalpha*hatvarthetas^alpha-kappasbeta*hatvarthetas^beta;

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
% sys(4)——sys(12)
sys = [sys,xr,yr,thr,vr,wr,ep,Rn,S,Sd];
if Rn==1
    th = thNormalization(th+pi)-pi;
end
% sys(13)——sys(17)
sys = [sys,x,y,th,v,w];
% sys(18)——sys(29)
sys = [sys,-Lw,Lw,es,-rhos,rhos,the,eth,-rhoth,+rhoth,0,0,0];
% sys(30)——sys(31)
sys = [sys,uv,uw];
% sys(32)——sys(35)
sys = [sys,zetaps,dhatvarthetas,zetapth,dhatvarthetath];
% sys(36)
sys = [sys,ds];
