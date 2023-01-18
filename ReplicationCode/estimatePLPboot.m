%[beta,sterr,error2]=estimatePLPboot(ydata,xdata,lag)
% estimate regressions in the bootsrap (without cleaning variables for
% missing obs. AT 2021
function [beta,sterr,error2]=estimatePLPboot(ydata,xdata,lag)


%----------------------------------------------------------------------------------------------------
% MANIPULATE DATA
%---------------------------------------------------------------------------------------------------
X1=xdata; 
X=X1;
Y=ydata;
[~, Ky]=size(ydata);



for k =1:Ky
        res(k).beta          = (X1'*X1)^(-1)*X1'*Y(:,k);   
        error(k).ols         =Y(:,k)-X1*res(k).beta;
        errorunres(k).ols    = error(k).ols;
        variance(k).ols = (error(k).ols'*error(k).ols)/(size(Y(:,k),1)-Ky*lag-1);   
        std(k).error    = (diag((X'*X)^(-1))*variance(k).ols).^0.5;
end

beta=zeros(size(res(Ky).beta,1),Ky);
sterr=zeros(size(std(Ky).error,1),Ky);
errors2=zeros(size(error(Ky).ols,1),Ky);

for k=1:Ky
    error2(:,k)  = error(k).ols;
    beta(:,k)=res(k).beta;
    sterr(:,k)=std(k).error;
end


%=====================================================================================================