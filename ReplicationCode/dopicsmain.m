function dopicsmain(DIFF,DIFFSTD,MA,CUMSTD,shock,values,neq, Names)

% Creates (cumulative) IRFs ordered by values at which the interactions are
% supposed to be evaluated.
ende=size(values,2)/2;  
count2 = 0;
for start = 1:ende

for i = 1:length(shock)
position = shock(:,i);

fig_no=i+count2;
figure(fig_no);
                                                                
var=0;
count=0;

 

for numpic=1:neq
var=var+1;

for r=1:2
IR(:,r)=MA(var,position,r+count2,:);
CR(:,1,r)=CUMSTD(var,position,1,r+count2,:);
CR(:,2,r)=CUMSTD(var,position,2,r+count2,:);
end

IR(:,3)=DIFF(var,position,start,:);
CR(:,1,3)=DIFFSTD(var,position,1,start,:);
CR(:,2,3)=DIFFSTD(var,position,2,start,:);



impcompareoutout(IR,CR,var,neq,numpic,position,Names); 
end
end
count=count+2;
count2=count2+2;
end



set(figure(1), 'PaperUnits', 'centimeters');
set(figure(1), 'PaperPosition', [0 0 20 15])

set(figure(2), 'PaperUnits', 'centimeters');
set(figure(2), 'PaperPosition', [0 0 20 15])

set(figure(3), 'PaperUnits', 'centimeters');
set(figure(3), 'PaperPosition', [0 0 20 15])

set(figure(4), 'PaperUnits', 'centimeters');
set(figure(4), 'PaperPosition', [0 0 20 15])



end