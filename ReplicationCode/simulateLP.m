% [IRF LAGEFF]=simulateLP(beta,neq,elast,shock2,restr,numatt,futsim)
% Computes the effects of plans from the betas on variables and Theta
% AT 2021

function [IRF LAGEFF]=simulateLP(beta,neq,elast,shock2,restr,numatt,futsim)

natt = sum(restr)/length(shock2)-1; % number of expected movements for type of plan
nev = length(shock2);
betasim = [];
betaforsim = [];
betaforsim2 = [];

for i = 1:nev
    elastsim = [];
    betasel = [];
    elastsim = [1 elast(:,1+natt*(i-1):natt*i)];
    betasel = [beta(i,:) ; beta(1+nev+natt*(i-1):nev+natt*i,:)]; 
    betasim = [betasim; elastsim * betasel];    
end

if isempty(futsim) == 0
    
betaforsim = [];
betaforsim2 = [];
count = 2;
for i=1:numatt
    betaforsel = [];
    betforsel = beta(i+nev*(1+natt),:);
    betaforsim = [betaforsim; elastsim(count) * betforsel];  
    if mod(i,nev) == 1
        count = count;
    else
        count = count +1;
    end
end  

for h = 1:natt
    betaforsim2(:,:,h) = betaforsim(1+(h-1)*nev:h*nev,:)';
end
LAGEFF(:,:,:)=betaforsim2(:,:,:);
else
    LAGEFF=[];
end
% OUTPUT
%-------------------------------------------------------------
% Storing IRFs
      
IRF(1:neq,1:nev)=betasim(1:nev,:)'; % in std dev.



%=============================================================