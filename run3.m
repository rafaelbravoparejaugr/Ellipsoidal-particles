global a b req alpha betarho betarho2 betaom betaom2


options=optimset('MaxFunEvals',1000000,'MaxIter',100000,'TolFun',5E-10);   % Option to display output

x0(1:2)=-2.02;

xn1 = fsolve(@geo3,x0,options);  

L = 2*req;

figure(1)

%%%%%CIRCULO 1
xc1=0;
yc1=0;
r11=req;
r12=req;
thetam1=0;

n = 50;  k=0:n;  fi=2*pi*k/n; 
x=xc1+r11*cos(fi)*cos(thetam1) - r12*sin(fi)*sin(thetam1); 
y=yc1+r11*cos(fi)*sin(thetam1) + r12*sin(fi)*cos(thetam1);
axis([-1 2 -1 1]);
plot(x,y,'r');
axis([-1 2 -1 1]);    
hold('on')    
axis([-1 2 -1 1])

hold('on')

%%%%%CIRCULO 2
xc2=L;
yc2=0;
r21=req;
r22=req;
thetam2=0;

n = 50;  k=0:n;  fi=2*pi*k/n; 
x=xc2+r21*cos(fi)*cos(thetam2) - r22*sin(fi)*sin(thetam2); 
y=yc2+r21*cos(fi)*sin(thetam2) + r22*sin(fi)*cos(thetam2);
axis([-1 2 -1 1]);
hold('on')
plot(x,y,'g');
axis([-1 2 -1 1]);    
hold('on')    
axis([-1 2 -1 1])


%%%%%ELIPSE

xp = a*cos(xn1(1)); 
yp = b*sin(xn1(1));
d1 = sqrt(xp^2+yp^2);

xce=L-(d1*cos(betaom)+req*cos(betarho));
yce=-(d1*sin(betaom)+req*sin(betarho));


r1=a;
r2=b;
thetame=alpha;


n = 50;  k=0:n;  fi=2*pi*k/n; 
x=xce+r1*cos(fi)*cos(thetame) - r2*sin(fi)*sin(thetame); 
y=yce+r1*cos(fi)*sin(thetame) + r2*sin(fi)*cos(thetame);
hold('on')
axis([-1 2 -1 1]);
hold('on')
plot(x,y);
axis([-1 2 -1 1]);    
hold('on')    
axis([-1 2 -1 1])


 






