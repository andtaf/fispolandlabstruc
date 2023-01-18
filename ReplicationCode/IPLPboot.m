% [BETAMAT ERRORMAT STERRMAT]=IPLPboot(parametric,Idata,ydata,xdata,errors,beta,number,I,lag,b)
% bootstrap AT 2021

function [BETAMAT ERRORMAT STERRMAT]=IPLPboot(parametric,Idata,ydata,xdata,errors,beta,number,I,lag,b)
                                                                    

% DETERMING ESTIMATION APPROACH
%-------------------------------------------------------------------------------------------------------------------
method1 = 'nearest';
method2 = 'nearest';
xsel=xdata;
yest=ydata - errors;
n = size(errors,1);
k = fix(n/b);
J = round(1+(n-b)*rand(1,k));
time = size(ydata,1)/I;
nint = size(Idata,2);

if time ~= floor(time);
%    display('Error: The time dimension is not identical across individuals!!')
end

%allocate space
BETAMAT(number).BETAMAT=zeros(size(beta,1),size(beta,2)); 
ERRORMAT(number).ERRORMAT= zeros(size(errors,1),size(errors,2));
STERRMAT(number).STERRMAT=zeros(size(beta,1),size(beta,2));


% CONSTRUCTING ARTIFICIAL DATA
%--------------------------------------------------------------------------
[Ty Ky]=size(ydata);
[Tx Kx]=size(xdata);
[Ti  Ki]=size(Idata);

% START BOOTSTRAP
%-------------------------------------------------------------------------
for nobs = 1:number
clear yart errors_new
    
% DRAW NEW ERRORS from OLD ONES if parametric =0 else from std normal
%-------------------------------------------------------------------------------------------
for ercol = 1:Ky
    if parametric==0
        c_errors = errors(:,ercol);
    else           
        c_errors = normrnd(0,kron(ones(size(errors(:,ercol),1),1),std(errors(:,ercol))));
    end
           
% DRAWING ERRORS
    erros_new = [];
    errors_indicator = [];
    psi = 1;
        while psi < k+1
        ind=randsample(n-b,1,true);
        block = [];
        for bl = 1:b
            block = [block ; ind+bl-1];
        end
        errors_indicator = [errors_indicator ; block];
        psi = psi+1;
        end
    if size(yest,1) == size(errors_indicator,1)
    else
       errors_adj = randsample(size(c_errors,1),(size(yest,1) - size(errors_indicator,1)),true);
       errors_indicator = [errors_indicator ; errors_adj];
    end
    for j=1:size(errors_indicator,1);errors_new(j,:)=c_errors(errors_indicator(j,1),:);end;
   
   c_indicator(:,1)=randsample(size(c_errors,1),size(c_errors,1),true);                                              
   for j=1:size(c_indicator,1);c_new_error(j,:)=c_errors(c_indicator(j,1),:);end;

   yart(:,ercol) = yest(:,ercol) + errors_new;
end 

% REPEAT VAR ESTIMATION
%-----------------------------------------------------------------------------
[betaX,sterrX,errorsX]=estimatePLPboot(yart,xsel,lag);

% STORE RESULTS
%----------------------------------------------------------------------------
BETAMAT(nobs).BETAMAT=betaX; 
ERRORMAT(nobs).ERRORMAT =errorsX; 
STERRMAT(nobs).STERRMAT=sterrX; 

clear sterrX betaX errorsX
end		% end of bootstrapping

%===========================================================================

