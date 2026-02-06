[fich1,texto1] = fopen('dts','w');
        
dt=5.0e-4;

nfich=1000;

for I=1:nfich
    
fprintf(fich1,'%e \n',dt*(I-1));
  
I
end %I

fclose(fich1)
