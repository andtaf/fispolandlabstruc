

function evalbeta=evalplp(beta,values,lag,x,restr,numatt,futsim,nshock)

nint=size(values,1);	 % number of interaction terms
neq=size(beta,2);		 % number of equations
nnlin = sum(restr(1,:));
nvar = size(x,2);
if values == 0
    values = [];
    nint= 0;
end
values=[1
	values];
betaeval = [];
% Evaluates the betas at the wished values 


for k=1:neq
for i=1:nnlin
betaeval(i,k)=beta(1+(1+nint)*(i-1):(1+nint)*i,k)' * values;
end
end

if isempty(futsim) == 0
    
for i=1:numatt
	betaeval(i+nnlin,k)=beta(futsim(:,1+(i-1)*length(values):i*length(values)),k)' * values;
end
end

if nnlin == 0
    betaeval = beta;
end

evalbeta=betaeval;


% for i=1:neq
%     constanteval(1,i)=beta(size(beta,1)-nint:size(beta,1),i)'*values;
% end
% evalbeta=[evalbeta
%           constanteval];
      %============================================================================


