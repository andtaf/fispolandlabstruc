function [MEDIAN, STD] = distr(SIMULI,values,neq,shock2,centered,pct,gamma)

STDU = [];
STDL = [];
MEDIAN = [];

  % Storing Results
for val=1:values
    for i=1:neq
    for j=1:length(shock2)
        for h = 1 : gamma
        if centered==0
            [STDU(i,j,val,h), STDL(i,j,val,h)]=minprctile(SIMULI(i,j,val,h,:),pct);
        elseif centered==1
            STDU(i,j,val,h)=prctile(SIMULI(i,j,val,h,:),100-pct/2);
            STDL(i,j,val,h)=prctile(SIMULI(i,j,val,h,:),pct/2);
        end
    MEDIAN(i,j,val,h)=median(SIMULI(i,j,val,h,:));
    end
    end
    end
    
end
STD(:,:,1,:,:)=STDU;
STD(:,:,2,:,:)=STDL;
end