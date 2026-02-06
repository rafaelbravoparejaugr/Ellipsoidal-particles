

function F=geo3(x)

global a b req alpha betarho betarho2 betaom betaom2


%x(1) -> betapsi1
%x(2) -> betapsi2

betarho  = atan2(a*sin(x(1)),b*cos(x(1))) + alpha;
betarho2 = atan2(a*sin(x(2)),b*cos(x(2))) + alpha;


xp = a*cos(x(1)); 
yp = b*sin(x(1));

d1 = sqrt(xp^2+yp^2);

xp = a*cos(x(2)); 
yp = b*sin(x(2));

d2 = sqrt(xp^2+yp^2);


betaom   = atan2(b*sin(x(1)),a*cos(x(1))) + alpha;
betaom2  = atan2(b*sin(x(2)),a*cos(x(2))) + alpha;


F(1) = -req*cos(betarho)-d1*cos(betaom)+d2*cos(betaom2)+req*cos(betarho2)+2*req;

F(2) = -req*sin(betarho)-d1*sin(betaom)+d2*sin(betaom2)+req*sin(betarho2);


d1
d2
betarho
betarho2




    
    
    
    
    
    
    
    
    
    
    
    