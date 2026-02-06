clear all;

%datos fijos del problema
dt     = 0.001;
d0     = 0.005;
g      = -9.81;
ji     = 0.25;
k      = 1;
ndisks = 3;
nbeta  = 0.5;
ngamma = 1.0;
nu     = 10^(-3);
phi    = 15.0; %angulo de rozamiento en grados!!!!
rhos   = 2500;
rhof   = 10;
rug    = 0.1;


%datos del programa

nstart = 1;
nfiles = 1000;

%geometria
f      = 0.45;
a      = 0.01;   %semieje mayor
b      = f*a;   %semieje menor
req    = a;   %radio de las esferas
alpha  = 45/180*pi; %inclinacion del elipsoide
%betaom = 350/180*pi-alpha; %pivoting angle
%beta   = betaom;
aincl  = 0;%angulo de inclinacion del lecho
aa     = a;
bb     = b;
rreq   = req;

%abrir ficheros de salida resultplot para gnuplot
[fichresults,texto]     = fopen('fichresults','w');
[fichresultsplot,texto] = fopen('fichresultsplot','w');

%I es el factor de escala
escinic = 0.7;
escincr = 0.02;
escfin  = 0.7;

% velocidad us inicial e incremento
usinic  = 0.8;
usincr  = 0.005;

run3;

N=1; %contador inicial para las variables results
for I=escinic:escincr:escfin
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fichero .tem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fichtem,texto] = fopen('cont.tem','w');

a   = aa*I;
b   = bb*I;
c   = b;
req = rreq*I;
tol = 1.0*req;
r11 = req;
r12 = req;
r21 = req;
r22 = req;


%llamada a la rutina que calcula los centros de esferas y elipse


p1 = [xc2*I;yc2*I];
p2 = [xc1*I;yc1*I];
p3 = [xce*I;yce*I];

mrot = [cos(aincl) -sin(aincl); sin(aincl) cos(aincl)];

p1rot = mrot*p1;
p2rot = mrot*p2;
p3rot = mrot*p3;

fprintf(fichtem,'%i 0 0 0 0 0 0 0 0 0 %i \n',ndisks,ndisks);

fprintf(fichtem,'%f %f %f %f %f  \n',p1rot(1),p1rot(2),r11,r12,thetam1);
fprintf(fichtem,'%i %i %i %i \n',1,1,1,1);
fprintf(fichtem,'%f %f %f %f %f  \n',p2rot(1),p2rot(2),r21,r22,thetam2);
fprintf(fichtem,'%i %i %i %i \n',1,1,2,1);
fprintf(fichtem,'%f %f %f %f %f  \n',p3rot(1),p3rot(2),a,b,thetame);
fprintf(fichtem,'%i %i %i %i \n',1,1,3,1);

fprintf(fichtem,'%i %i \n',1,5);
fprintf(fichtem,'%i    \n',1);
fprintf(fichtem,'%i %i %i \n',4,0,0);

p1aux=[(xc1*I-r11);0.25*r12];
p1auxrot=mrot*p1aux;
fprintf(fichtem,'%f %f %f \n',p1auxrot(1),p1auxrot(2),2);
p1aux=[(xc1*I-r11);-r12];
p1auxrot=mrot*p1aux;
fprintf(fichtem,'%f %f %f \n',p1auxrot(1),p1auxrot(2),2);
p1aux=[(xc2*I+r21);-r22];
p1auxrot=mrot*p1aux;
fprintf(fichtem,'%f %f %f \n',p1auxrot(1),p1auxrot(2),2);
p1aux=[(xc2*I+r21);0.25*r22];
p1auxrot=mrot*p1aux;
fprintf(fichtem,'%f %f %f \n',p1auxrot(1),p1auxrot(2),2);
fprintf(fichtem,'%f    \n',0);
fprintf(fichtem,'%f      ',0);

st  = fclose(fichtem); 

dist = 0;

us = usinic;

while dist<tol
    
us  = us + usincr ;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fichero .par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fichpar,texto] = fopen('cont.par','w');

fprintf(fichpar,'1 \n');
fprintf(fichpar,'%f  \n',nbeta);
fprintf(fichpar,'%f  \n',ngamma);
fprintf(fichpar,'%f  \n',ji);
fprintf(fichpar,'%f  \n',rhof);
fprintf(fichpar,'%f  \n',rug);
fprintf(fichpar,'%f  \n',nu);
fprintf(fichpar,'%f  \n',us);
fprintf(fichpar,'100000000  \n');
fprintf(fichpar,'%f  \n',4/3*pi*req^2*50000000*abs(g));
fprintf(fichpar,'%f \n',dt);
fprintf(fichpar,'%f \n',d0);
fprintf(fichpar,'0 0    \n');
fprintf(fichpar,'%f %f  \n',0,g);
fprintf(fichpar,'0 \n');
fprintf(fichpar,'%f \n', rhos);
fprintf(fichpar,'2 \n');
fprintf(fichpar,'1    %f 0 \n',phi);
fprintf(fichpar,'2    %f 0 \n',phi);
fprintf(fichpar,'-1 \n');
fprintf(fichpar,'-1 \n');
fprintf(fichpar,'-1');

st  = fclose(fichpar);

%llamada a dda

%pause

[stat1,stat2] = unix('./dda');

%lectura de resultados

for J = nstart : nfiles

KK  = strcat('Onor',int2str(J));
KK1 = strcat(KK,'.rs');
KK1
  
[fichgeo,texto] = fopen(KK1,'r');
 
            nd(1:11) = fscanf(fichgeo,'%lg',11);
                           
            for JJ = 1 : nd(1)
            aux(1:9) = fscanf(fichgeo,'%lg',9);
            Coord(J,JJ,1) = aux(1);
            Coord(J,JJ,2) = aux(2);
            Coord(J,JJ,3) = aux(3);
            Coord(J,JJ,4) = aux(4);
            Coord(J,JJ,5) = aux(5);
            end % I                      
                        
            final = J;

            aux(1:2) = fscanf(fichgeo,'%lg',2);
            
            mm=2*aux(1)+2+3*aux(2)-2;
            
            aux1(1:mm) = fscanf(fichgeo,'%lg',mm);
                       
            time(J) = aux1(mm);
                        
            st  = fclose(fichgeo);             

            dist=sqrt((Coord(J,3,2)-Coord(1,3,2))^2+(Coord(J,3,1)-Coord(1,3,1))^2);
            
            if dist > tol
            break;
            end %if dist > tol
                
end % J



us
tol
dist
pause

if dist<tol
unix('rm *.rs');
unix('rm *.df');
unix('rm *.par');
end


end %while dist<tol

usinic  = us-5*usincr;

% if (vcr-7*vincr)<0
% vcrinic=0;
% else
%usinic=us-10*usincr;
% vcrinic=vcr-2*vincr;
% end %if
% vcrinic=0.05;

unix('rm *.tem');

vol = 4/3*pi*a*b*b;

%diametro equivalente

deq  = (6*vol/pi)^(1/3);

if us*deq*rhof/nu<4
cad   = 'laminar';
else
cad   = 'turbulent';
end


result(N,1) = I;
result(N,2) = us;
result(N,3) = us;
result(N,4) = rhof*us*deq/nu;
result(N,5) = us^2*rhof/((rhos - rhof)*abs(g)*deq);

fprintf(fichresults,'%e %e %e %e %e %s \n',result(N,1),result(N,2),result(N,3),result(N,4),result(N,5), cad);

fprintf(fichresultsplot,'%e %e %e %e %e \n',result(N,1),result(N,2),result(N,3),result(N,4),result(N,5));

N = N + 1;


end %I

fclose(fichresults);
fclose(fichresultsplot);


              
