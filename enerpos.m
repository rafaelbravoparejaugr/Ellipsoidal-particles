%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Macro que obtiene las posiciones y nerg?as del problema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear();


time(1)=0;

nstart=1;
nfiles=2100;
nevery=50;

Ener(1:nfiles)=0;

figure(1)

nn=0;

for JJ = nstart:nevery: nfiles

 nn=nn+1;
   
 KK  = Strcat('Onor',int2str(JJ));
 KK1 = Strcat(KK,'.rs');
 KK1
 
 LL  = Strcat('Onor',int2str(JJ));
 LL1 = Strcat(LL,'.df');
 LL1  
 
 [fich1,texto] = fopen(KK1,'r');
 [fich2,texto] = fopen(LL1,'r');
 
            nd(1:11) = fscanf(fich1,'%lg',11);
                     
      
            for II = 1 : nd(1)
            aux(1:9) = fscanf(fich1,'%lg',9);
            Coord(nn,II,1) = aux(1);
            Coord(nn,II,2) = aux(2);
            Coord(nn,II,3) = aux(3);
            Coord(nn,II,4) = aux(4); 
            Coord(nn,II,5) = aux(5);             
            end % II                      
            
            
            final = nn;

            aux(1:2) = fscanf(fich1,'%lg',2);
            
            mm=2*aux(1)+2+3*aux(2)-2;
            
            aux1(1:mm) = fscanf(fich1,'%lg',mm);
                       
            time(nn) = aux1(mm);
            
            st  = fclose(fich1);  
            
%Lectura de fichero de velocidades          
            
            aux(1:12) = fscanf(fich2,'%lg',12);
            
            rho=aux(11);
            
            for II=1:aux(12)
            aux(1:3) = fscanf(fich2,'%lg',3);             
            end %II  
            
            aux(1) = fscanf(fich2,'%lg',1); 
            
            Lx(nn)=0;
            Ly(nn)=0;
            JJang(nn)=0;
            
            for II=1:aux(1)
            Velo(nn,II,1:4) = fscanf(fich2,'%lg',4);
            Lx(nn) = Lx(nn) + Velo(nn,II,2);
            Ly(nn) = Ly(nn) + Velo(nn,II,3); 
            JJang(nn) = JJang(nn) + Velo(nn,II,4);             
            Velo(nn,II,2)=Velo(nn,II,2)/(rho*pi*Coord(nn,II,3)*Coord(nn,II,4));
            Velo(nn,II,3)=Velo(nn,II,3)/(rho*pi*Coord(nn,II,3)*Coord(nn,II,4));   
            Velo(nn,II,4)=Velo(nn,II,4)/(0.25*rho*pi*Coord(nn,II,3)*Coord(nn,II,4)*(Coord(nn,II,3)^2+Coord(nn,II,4)^2));               
            end %II               
          
            
            Ener(nn)=0;
            
            for II=1:aux(1)
            Ener(nn)=Ener(nn)+0.5*(rho*pi*Coord(nn,II,3)*Coord(nn,II,4)*Velo(nn,II,2)^2+rho*pi*Coord(nn,II,3)*Coord(nn,II,4)*Velo(nn,II,3)^2+0.25*rho*pi*Coord(nn,II,3)*Coord(nn,II,4)*(Coord(nn,II,3)^2+Coord(nn,II,4)^2)*Velo(nn,II,4)^2)-aux(8)*Coord(nn,II,1)-aux(9)*Coord(nn,II,2);
            end %II
 
for II = 1 : nd(1)
    xc=Coord(nn,II,1);
    yc=Coord(nn,II,2);
    r1=Coord(nn,II,3);
    r2=Coord(nn,II,4);
    thetam=Coord(nn,II,5);
    n = 50;  k=0:n;  fi=2*pi*k/n; 
    x=xc+r1*cos(fi)*cos(thetam) - r2*sin(fi)*sin(thetam); 
    y=yc+r1*cos(fi)*sin(thetam) + r2*sin(fi)*cos(thetam);
%    axis([-0.5 2 -0.5 2]);
    plot(x,y);
%    axis([-0.5 2 -0.5 2]);    
hold('on')   
%    axis([-0.5 2 -0.5 2]);
end % II  
M(JJ)=getframe;
%    axis([-0.5 2 -0.5 2]);
hold('off')

end % JJ








              


