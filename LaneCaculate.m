function [s] = LaneCaculate(x,y,SN,SCycles)

global vd isClockwise Lr Lx Ly Ll Lt SR Ssum
% Dx = Lx/2-Ly-4*r+2^0.5*((Ll+4*r)*(Ly-Ll+2*r))^0.5;%交点
% Dy = Dx-Lx/2+Ly/2;
if (sqrt((x-(Lx/2-2*Lr-Ll))^2+(y-(Ly/2-Lt))^2)-Lr>=y-(-Ly/2-Lr) && x-Lx/2<=-y-Ly/2 && Lx/2-2*Lr-Ll<=x && x<=Lx/2) ||...
    (y<=(-Lt-2*Lr)/2 && -(Lx/2-2*Lr-Ll)<=x && x<=Lx/2-2*Lr-Ll) ||...
    (sqrt((x+(Lx/2-2*Lr-Ll))^2+(y-(Ly/2-Lt))^2)-Lr>=y-(-Ly/2-Lr) && x+Lx/2>=y+Ly/2 && -Lx/2<=x && x<=-(Lx/2-2*Lr-Ll))
    RegionNum = 1;
    xd = x;
    yd = -Ly/2-Lr;
    thd = -pi;
    wd = 0;
    ep = y-yd;
    sn = Lx/2-x;
elseif x<=-Lx/2 && y<=-Ly/2
    RegionNum = 2;
    rx = -Lx/2;
    ry = -Ly/2;
elseif x>=Lx/2 && y<=-Ly/2
    RegionNum = -2;
    rx = Lx/2;
    ry = -Ly/2;
elseif (sqrt((x+(Lx/2-2*Lr-Ll))^2+(y-(Ly/2-Lt))^2)-Lr>=x+Lx/2+Lr && y+Ly/2>=x+Lx/2 && -Ly/2<=y && y<=Ly/2-Lt) ||...
        (x+Lx/2<=-y+Ly/2 && Ly/2-Lt<=y && y<=Ly/2)
    RegionNum = 3;
    xd = -Lx/2-Lr;
    yd = y;
    thd = pi/2;
    wd = 0;
    ep = x-xd;
    sn = y+Ly/2;
elseif (sqrt((x-(Lx/2-2*Lr-Ll))^2+(y-(Ly/2-Lt))^2)-Lr>=-x+Lx/2+Lr && y+Ly/2>=-x+Lx/2 && -Ly/2<=y && y<=Ly/2-Lt) ||...
        (x-Lx/2>=y-Ly/2 && Ly/2-Lt<=y && y<=Ly/2)
    RegionNum = -3;
    xd = Lx/2+Lr;
    yd = y;
    thd = -pi/2;
    wd = 0;
    ep = xd-x;
    sn = Ly/2-y;
elseif x<=-Lx/2 && y>=Ly/2
    RegionNum = 4;
    rx = -Lx/2;
    ry = Ly/2;
elseif x>=Lx/2 && y>=Ly/2
    RegionNum = -4;
    rx = Lx/2;
    ry = Ly/2;
%5 or -5
elseif (y-Ly/2>=x+Lx/2-Ll && x+Lx/2>=-y+Ly/2 && -Lx/2<=x && x<=-Lx/2+Ll) ||...
        (y-Ly/2>=-x-Lx/2+Ll && x-Lx/2<=y-Ly/2 && Lx/2-Ll<=x && x<=Lx/2)
    if x<0
        RegionNum = 5;
        sn = x+Lx/2;
    else
        RegionNum = -5;
        sn = x-(Lx/2-Ll);
    end
    xd = x;
    yd = Ly/2+Lr;
    thd = 0;
    wd = 0;
    ep = -y+yd;
elseif -Lx/2+Ll<=x && x<=0 && y>=Ly/2
    RegionNum = 6;
    rx = -Lx/2+Ll;
    ry = Ly/2;
elseif 0<=x && x<=Lx/2-Ll && y>=Ly/2
    RegionNum = -6;
    rx = Lx/2-Ll;
    ry = Ly/2;
elseif y-(Ly/2-Lt)<=x-(-Lx/2+Ll+2*Lr) && y-(Ly/2-Lt)<=-x-(-Lx/2+Ll+2*Lr) &&...
        -Lx/2+2*Lr+Ll<=x && x<=Lx/2-2*Lr-Ll && (-Lt-2*Lr)/2<y
    RegionNum = 9;
    xd = x;
    yd = Ly/2-Lt-Lr;
    thd = 0;
    wd = 0;
    ep = -y+yd;
    sn = x-(-Lx/2+Ll+2*Lr);
else%7 8 and -7 -8
    if x<0
        if y>Ly/2-Lt
            RegionNum = 7;
            xd = -Lx/2+Ll+Lr;
            yd = y;
            thd = -pi/2;
            wd = 0;
            ep = -x+xd;
            sn = Ly/2-y;
        else
            RegionNum = 8;
            rx = -Lx/2+Ll+2*Lr;
            ry = Ly/2-Lt;
        end
    else
        if y>Ly/2-Lt
            RegionNum = -7;
            xd = Lx/2-Ll-Lr;
            yd = y;
            thd = pi/2;
            wd = 0;
            ep = x-xd;
            sn = y-(Ly/2-Lt);
        else
            RegionNum = -8;
            rx = Lx/2-Ll-2*Lr;
            ry = Ly/2-Lt;
        end
    end
end
switch RegionNum
    case {2,4,6,-6,-4,-2}
        th = atan2(y-ry,x-rx);
        xd = rx+Lr*cos(th);
        yd = ry+Lr*sin(th);
        thd = th-pi/2;
        wd = -vd/Lr;
        ep = Lr-sqrt((x-rx)^2+(y-ry)^2);
        sn = (pi/2-thNormalization4(th))*Lr;
    case {8,-8}
        th = atan2(y-ry,x-rx);
        xd = rx+Lr*cos(th);
        yd = ry+Lr*sin(th);
        thd = th+pi/2;
        wd = vd/Lr;
        ep = sqrt((x-rx)^2+(y-ry)^2)-Lr;
        sn = thNormalization4(th)*Lr;
end
if RegionNum==1
    Sn = sn;
elseif RegionNum>1
    Sn = sum(SR(1:RegionNum-1))+sn;
elseif RegionNum==-8
    Sn = sum(SR)+sn;
else
    Sn = sum(SR)+sum(SR(-(RegionNum-1):8))+sn;
end

% Cycle to continuous
if Sn-SN>Ssum/2
    SCycles = SCycles-1;
elseif Sn-SN<-Ssum/2
    SCycles = SCycles+1;
end
SN = Sn;
S = SCycles*Ssum+Sn;

if isClockwise==1
    thd = thNormalization(thd);
else
    thd = thNormalization(thd+pi);
end

s=[xd,yd,thd,vd,wd,isClockwise*ep,RegionNum,S,SN,SCycles];
