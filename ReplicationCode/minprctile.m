%
%function [U L]=minprctile(X,pctile)
%
% =============================================================================
% Returns the upper and lower value of the values in X which define the 
% (1-pctile) perecent of the values. Rather than taking the 2.5% and 97.5 
% limits in the case of pctile=5, the values L and U are chosen to minimze 
% the distance of U and L mainting the area between U and L equal to 1-pctile.  
% ------------------------------------------------------------------------------
% Calls on the function prctile
% ============================================================================
% See also prctile
%
% Written by P. Towbin and S. Weber 21.01.2009
%

function [U L]=minprctile(X,pctile)


% "exact" is the decimeter unit at which the minimization is to be done
% the smaller it is chosen the finer the grid to fidn the minimizing area.
exact  = 0.1;

length = pctile/exact; 

lower=zeros(length+1,1);
upper=zeros(length+1,1);

for i=1:length+1
    pctile2 = exact*(i-1);
    lower(i,1)=prctile(X,pctile2);
    upper(i,1)=prctile(X,100-pctile+pctile2);
end
dist  = upper - lower;
[Y I] = min(dist);

L=lower(I,1);
U=upper(I,1);

% ============================================================================
