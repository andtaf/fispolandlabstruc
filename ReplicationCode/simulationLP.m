%function [DIFF,SIMUL] = simulationLP(BETAMAT,elast,gamma,values,neq,lag,x,restr,shock2,numatt,futsim)
% Writes vectors for the simulations and their differences
% AT 2021

function [DIFF,SIMUL] = simulationLP(BETAMAT,elast,gamma,values,neq,lag,x,restr,shock2,numatt,futsim)

numatt2 = numatt/size(shock2,2);

BETAEV = [];
sizediff = size(values,2)/2;

if size(values,2) ~= 0

for val=1:size(values,2)
   for h=1:gamma
       BETAEV(:,:,val,h) = evalplp(BETAMAT(h).BETAMAT(:,:),values(:,val),lag,x,restr,numatt,futsim,size(shock2,2));
       [SIMULA LAGEFFA] = simulateLP(BETAEV(:,:,val,h),neq,elast,shock2,restr,numatt,futsim);
       SIMULB(:,:,val,h) = SIMULA(:,:);
       LAGEFFB(:,:,val,:,h) = LAGEFFA(:,:,:);
   end  
end

SIMUL(:,:,:,:)=SIMULB(:,:,:,:);

 
if isempty(futsim) == 0
for ii = 1:numatt2
    SIMUL(:,:,:,ii+1)=SIMULB(:,:,:,ii+1)+LAGEFFB(:,:,:,ii,1);
end
end
% Storing Results
for val=1:sizediff
for h=1:gamma
    DIFF(:,:,val,h) = SIMUL(:,:,1+2*(val-1),h) - SIMUL(:,:,2+2*(val-1),h);
end
end

else
    
val=1;
   for h=1:gamma
       BETAEV(:,:,val,h) = evalplp(BETAMAT(h).BETAMAT(:,:),0,lag,x,restr,numatt,futsim,size(shock2,2));
       [SIMULA LAGEFFA] = simulateLP(BETAEV(:,:,val,h),neq,elast,shock2,restr,numatt,futsim);
       SIMULB(:,:,val,h) = SIMULA(:,:);
       LAGEFFB(:,:,val,:,h) = LAGEFFA(:,:,:);
   end  


SIMUL(:,:,:,:)=SIMULB(:,:,:,:);

 
if isempty(futsim) == 0
for ii = 1:numatt2
    SIMUL(:,:,:,ii+1)=SIMULB(:,:,:,ii+1)+LAGEFFB(:,:,:,ii,1);
end
end
% Storing Results


val=1;
for h=1:gamma
    DIFF(:,:,val,h) = 0;
end


end

end
