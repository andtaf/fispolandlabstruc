% [beta,sterr,error2,Ysel,Xsel]=estimatePLP(ydata,xdata,I,lag,Idata,restr,demean,hor,cont,addlag)
% Estimate the interacted regression. AT 2021


function [beta,sterr,error2,Ysel,Xsel]=estimatePLP(ydata,xdata,I,lag,Idata,restr,demean,hor,cont,addlag)


%----------------------------------------------------------------------------------------------------
% MANIPULATE DATA
%---------------------------------------------------------------------------------------------------
xdata=xdata; 
ydata=ydata;
hor=hor;
[Ty Ky]=size(ydata);
[Tx Kx]=size(xdata) ;
[Ti  Ki]=size(Idata);
time = size(ydata,1)/I;

if time ~= floor(time)
    display('Error: The time dimension is not identical across individuals!!')
end

%---------------------------------------------------------------------------------------------------------
% CREATE INTERACTION TERMS
% ORDERED ACCORDING TO [X1 X1*I1  X1*I2 ... X1*IKi  X2 X2*I1 ...XKx*IKi]
%-----------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------
% CREATE REGRESSORS
%----------------------------------------------------------------------------------------------------

Y = [];
error2 = [];
Y = paneldiff(ydata,I,hor);
Y = [Y(hor:end,:); NaN(hor-1,Ky)];
Yn = [];
adjsize = size(Y,1)/I;
tempX=[];
tempY=[];
tempI=[];
%-----------------------------------------------------------------------------------------------------
% Demean
%----------------------------------------------------------------------------------------------------
Xh=[];


for nx=1:Kx %cycles over endogenous variables

   if restr(1,nx) == 1
      Xh=[Xh xdata(:,nx)];
   for nint=1:Ki %cycles over exogenous drivers of interaction
            
      Xh=[Xh xdata(:,nx).*Idata(:,nint)];  % this selects trimmed values for each country/lag/variable 
                                                                                    % triple (ydata selects country-specific data while trimmering
                                                                                    % them according to the lag, nx the variable) and multiplies it
                                                                                    % by the respective trimmed values for each country/lag/exg. 
                                                                                    % variable -> works stairwise I obtain a 
                                                                                    % (time-lag-1)*(1+num of int) matrix
   end
   else
      Xh=[Xh xdata(:,nx)];

   end     
end


X = Xh;

for jj = 1:size(X,2)
meanX(:,jj) = nanmean(X(:,jj));
end
for jj = 1:Ky
meanY(:,jj) = nanmean(Y(:,jj));
end
for jj = 1:Ki
meanI(:,jj) = nanmean(Idata(:,jj));
end


if isempty(Idata) == 1
    
    if demean ==1
    for kk=1:I
        newX = X(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryX = [];
        newY = Y(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryY = [];
        for jj = 1:size(X,2)
         % if estimations is fixed effect
            cntryX(:,jj)=newX(:,jj)-nanmean(newX(:,jj)) + meanX(:,jj);
        end
        for jj = 1:size(Y,2)
         % if estimations is fixed effect
            cntryY(:,jj)=newY(:,jj)-nanmean(newY(:,jj)) + meanY(:,jj);
        end
        
        tempX = [tempX
                 cntryX];           
        tempY = [tempY
                 cntryY]; 
    end
else  % if estimations is not fixed effect
    for kk=1:I
        newX = X(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryX = [];
        newY = Y(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryY = [];
        
        for jj = 1:size(X,2)
         % if estimations is fixed effect
            cntryX(:,jj)=newX(:,jj);
        end
        

        for jj = 1:size(Y,2)
            cntryY(:,jj)=newY(:,jj)
        end
        
        tempX = [tempX
                 cntryX];      
   
        tempY = [tempY
                 cntryY];             
    end
end

X = tempX;
Idt = [];
Y = tempY;


else
if demean ==1
    for kk=1:I
        newX = X(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryX = [];
        newIdt = Idata(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryI = [];
        newY = Y(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryY = [];
        for jj = 1:size(X,2)
         % if estimations is fixed effect
            cntryX(:,jj)=newX(:,jj)-nanmean(newX(:,jj)) + meanX(:,jj);
        end
        
        for jj = 1:size(Idata,2)
             cntryI(:,jj)=newIdt(:,jj)-nanmean(newIdt(:,jj)) + meanI(:,jj);
        end
        
        for jj = 1:size(Y,2)
         % if estimations is fixed effect
            cntryY(:,jj)=newY(:,jj)-nanmean(newY(:,jj)) + meanY(:,jj);
        end
        
        tempX = [tempX
                 cntryX];      
        tempI = [tempI
                 cntryI];        
        tempY = [tempY
                 cntryY]; 
    end
else  % if estimations is not fixed effect
    for kk=1:I
        newX = X(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryX = [];
        newIdt = Idata(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryI = [];
        newY = Y(1+(kk-1)*adjsize:kk*adjsize,:);
        cntryY = [];
        
        for jj = 1:size(X,2)
         % if estimations is fixed effect
            cntryX(:,jj)=newX(:,jj);
        end
        
        for jj = 1:size(Idata,2)
             cntryI(:,jj)=newIdt(:,jj);
        end

        for jj = 1:size(Y,2)
            cntryY(:,jj)=newY(:,jj)
        end
        
        tempX = [tempX
                 cntryX];      
        tempI = [tempI
                 cntryI];      
        tempY = [tempY
                 cntryY];             
    end
end

X = tempX;
Idt = tempI;
Y = tempY;

end

[Tx2 Kx2] = size(X);
for i=1:lag+1 %foreach lag creates a structrue Xn composed of the matrix of interacted regressors X
    Xn(i).X = [];
end


for i=1:lag+1 % cycles over lags
    if i == lag+1
    Xlag = [];
   
    Xlag = panellag(X(:,1:addlag*(1+Ki)), I, i-1);
     
    if i == lag+1 && cont ~= 0

            for actual_country=1:I
                  start=1+(actual_country-1)*time;
                  ende=actual_country*time;
                  Xa(start:ende-1,:)=X(1+start:ende,Kx2-(cont-1):Kx2);
                  Xa(ende,:)=nan;
            end
            Xlag=[Xlag Xa];
    end        
    else
    Xlag = [];
   
    Xlag = panellag(X, I, i-1);
     
    end
    XX(i).X = Xlag; %this saves a different matrix for each lag
end

%  added to include interaction -> trimmers the interaction variables data
Idata2=[];
Idata2 = Idt;

% Create intercept
Xd = ones(size(XX(1).X,1),1);
X=[];
% Create intercept interactions
Xinter=[];
for nint=1:Ki
   Xinter=[Xinter Xd.*Idata2(:,nint)];  %  added to include interaction
end

 X=[X Xinter]; %  added to include interaction
%don't include constant (for fixed effects)
%X=[Xinter];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS COMMAND CREATES THE MATRIX OF REGRESSORS                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Create Regressors (NOTE: THE ORDER IS REVERSED, i.e. XX(2).X is the first
% lag in a PVAR(2) and XX(1).X is the second lag!!)
for i=1:lag+1
    X = [XX(lag+2-i).X X];   
end


%---------------------------------------------------------------------------------------
% ADJUSTING FOR NaN
%---------------------------------------------------------------------------------------
tempX=[];
tempY=[];
tic;
for kk=1:I
    newX = X(1+(kk-1)*adjsize:kk*adjsize,:);
    helpY = Y(1+(kk-1)*adjsize:kk*adjsize,:);
    cntryX = [];
      
    idX = any(isnan(newX), 2);
    idY = any(isnan(helpY), 2);
    
    id = idX + idY;

    [poscan,~] = find(any(id ~= 0,2));
    
    newX(poscan,:) = [];
    helpY(poscan,:) = [];
    
    tempX = [tempX
             newX];
         
    tempY = [tempY
             helpY];
         
end

X = tempX;
Y = tempY;
toc;

clear tempX tempY

% Taking away the lag zero part

% if ord_first == 1
%  Idata_adj=X(:,1:Kx*(Ki+1)); %this is the matrix of contemporaneous effects!!!!!!!
%  X=X(:,1+Kx*(Ki+1):end);
% else
% end

%-----------------------------------------------------------------------------------------------------
% DEFINE RESTRICTIONS
%---------------------------------------------------------------------------------------------------
int=Ki; % Ki is the number of interacting variables
restr_new=restr;

        

%----------------------------------------------------------------------------------------------------
% OLS Eq-by-Eq ESTIMATION
%----------------------------------------------------------------------------------------------------
%saveX=X;
%saveY=Y;
% UNRESTRICTED ESTIMATION

X1=[X ones(size(X,1),1) ];

for k =1:Ky
        res(k).beta          = (X1'*X1)^(-1)*X1'*Y(:,k);   
        error(k).ols         =Y(:,k)-X1*res(k).beta;
        errorunres(k).ols    = error(k).ols;
        variance(k).ols = (error(k).ols'*error(k).ols)/(size(Y(:,k),1)-Ky*lag-1);   
        std(k).error    = (diag((X'*X)^(-1))*variance(k).ols).^0.5;
end

% for k=1:Ky
%         rres(k).beta     = res(k).beta-(X'*X)^(-1)*R(k).R'*((R(k).R*(X'*X)^(-1)*R(k).R')^(-1))*(R(k).R*res(k).beta-r(k).r);   
%         rerror(k).ols    = Y(:,k)-X*rres(k).beta;
%         rvariance(k).ols = (rerror(k).ols'*rerror(k).ols)/(size(Y(:,k),1)-Ky*lag-1);   
%         rstd(k).error    = (diag((X'*X)^(-1))*rvariance(k).ols).^0.5;
% 
% %     Yh = Idata_adj(:,1+(k-1)*(nint+1):k*(nint+1)); 
% %     X = [Yh  X];
% end

%--------------------------------------------------------------------------
% Testing joint (interaction and variable) Exogeneity (if imposed)
%--------------------------------------------------------------------------
% numberofres=(1+nint)*(Ky-1)*lag+nint+nint;
% for exn=1:Ky
%     Fvalue(1,exn,hor)=((rerror(exn).ols'*rerror(exn).ols-errorunres(exn).ols'*errorunres(exn).ols)/numberofres)*((size(Y(exn).Y,1)-Ky*lag-1)/(errorunres(exn).ols'*errorunres(exn).ols));
% end

%-----------------------------------------------------------------------------------------------------
% STORING RESULTS
%----------------------------------------------------------------------------------------------------

beta=zeros(size(res(Ky).beta,1),Ky);
sterr=zeros(size(std(Ky).error,1),Ky);
errors2=zeros(size(error(Ky).ols,1),Ky);
Ysel = Y;
Xsel = X1;

for k=1:Ky
    error2(:,k)  = error(k).ols;
    beta(:,k)=res(k).beta;
    sterr(:,k)=std(k).error;
end


%=====================================================================================================