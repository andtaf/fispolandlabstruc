
% This code estimates the non linear effects of fiscal plans according to
% labour market institutions. This code builds on Tobin and Weber 2009
% code of IPVAR, to which I'm deeply indebted.
% 
% The code requires
% 
% y: a vector of dependent variables
% 
% x: a vector of Xs
% 
% elast: the vector of elasticities (Theta in the paper)
% 
% Idata: the vector of interactions (Gamma in the paper)
% 
% I: the number of countries
% 
% hominterc: if = 1  inclusion of heterogeneous intercept
% 
% demean: if = 1 fixed effects
% 
% hor: lenght of IRF
% 
% lag: number of lags of X
% 
% parametric : = 0 to bootstraped CI, =1 for parametrics
% 
% restr: Declare with Xs need to be interacted [=1 interacted, =0 not interacted]
% 
% centered: = 0 median, = 1 mean
% 
% pct: = 10 90 CI
% 
% block: size of the block bootstrap
% 
% shock: Variables that compose the plan
% 
% shock2: how many shocks
% 
% rep: number of bootstrap replications
% 
% nint: number of interaction terms
% 
% neq: number of ebndogenous variables
% 
% nvar: number of explanatori variables
% 
% nnlin: number of interactions
% 
% cont: declare how many variables need to be enter contemporaneously (order them last in the vector of X
% 
% futsim: an indicator to exclude specific colums so that only plans announced at t-i and implemented at t are included in lag t-1

% The estimation is in two steps: first, a regression with boostrapped CI
% is estimated, as Eq. 5 in the paper. Second, the estimates are used to
% simulate the effects of plans with the help of Theta, the vector of
% elasticities.



%% Start

tic;
clear; 
close all;
cd 'C:\Users\at.eco\OneDrive - CBS - Copenhagen Business School\Non Linearities\ReplicationCode'
%% LOAD DATA
load data
%%%%%%%%%%%%%%%%
%% Define variables
gdp = table2array(data(:,'gdp'));
cpi = table2array(data(:,'cpi'));
lt_r = table2array(data(:,'lt_r'));
st_r = table2array(data(:,'st_r'));
gpb_pot = table2array(data(:,'gpb_pot'));
glob_cycle = table2array(data(:,'glob_cycle'));
cci = table2array(data(:,'cci'));
tbu = table2array(data(:,'tbu'));
ebu = table2array(data(:,'ebu'));
td_att1 = table2array(data(:,'td_att1'));
ed_att1 = table2array(data(:,'ed_att1'));
td_att2 = table2array(data(:,'td_att2'));
ed_att2 = table2array(data(:,'ed_att2'));
td_att3 = table2array(data(:,'td_att3'));
ed_att3 = table2array(data(:,'ed_att3'));
td_att4 = table2array(data(:,'td_att4'));
ed_att4 = table2array(data(:,'ed_att4'));
td_att5 = table2array(data(:,'td_att5'));
ed_att5 = table2array(data(:,'ed_att5'));
tfp = table2array(data(:,'tfp'));
wage = table2array(data(:,'wage'));
in = table2array(data(:,'in'));
cons = table2array(data(:,'cons'));
epl = table2array(data(:,'epl'));
brr = table2array(data(:,'brr'));
un_r = table2array(data(:,'un_r'));
lf_r = table2array(data(:,'lf_r'));
year = table2array(data(:,'year'));
country = table2array(data(:,'country'));
pmr = table2array(data(:,'pmr'));
emp = table2array(data(:,'emp'));
hours = table2array(data(:,'hours'));

%%  Preliminary steps

cont = 2; %declare how many variables need to be enter contemporaneously (order them last in the vector of X

rng('default') %set random number for replication

I = 17; % number of countries

fullx = [year tbu ebu td_att1 td_att2 td_att3 ed_att1 ed_att2 ed_att3 ...
    panellag(paneldiff(gdp,I,1),I,1) panellag(paneldiff(cpi,I,1),I,1) ...
    panellag(paneldiff(un_r,I,1),I,1)...
    panellag(paneldiff(gpb_pot,I,1),I,1) panellag(paneldiff(glob_cycle,I,1),I,1) ...
    panellag(paneldiff(lt_r,I,1),I,1)]; %vecor of Xs

% trimmer the dataset
idx = any(fullx==2019,2);
[del,~] = find(idx);
x=fullx;
x(del,:) = [];
x(:,1) = [];
c_id = country;
c_id(del,:) = [];

xtime = size(x,1)/17; %declare time observations per country

y = [gdp in cons cpi un_r./100 wage st_r./100 tfp ];
ynms = array2table(y,'VariableNames',{'Output', 'Investment', 'Consumption', 'Prices', 'Un. Rate', 'Wages','Interest Rate',  'TFP'});
y(del,:) = [];


%elasticities of theta

elast11 = regress(td_att1,tbu);

elast12 = regress(td_att2,tbu);

elast13 = regress(td_att3,tbu);

elast21 = regress(ed_att1,ebu);

elast22 = regress(ed_att2,ebu);

elast23 = regress(ed_att3,ebu);



elast = [ elast11 elast12 elast13 elast21 elast22 elast23 ]; % elasticities are estimated here as I want to simulate an identical plan 
                                            % (I don't want the *structure* of the plan to vary across states => differences 
                                            % of effects are only due to differences in marginal effects
                                            


numatt = length(elast);

% declare and normalize Gamma
fulltof = [ epl brr ];
tof = fulltof;
tof(del,:) = [];

epl(del,:) = [];
brr(del,:) = [];
pmr(del,:) = [];
pmr = normalize(pmr);

for jj = 1:size(tof,2)
factorXem(:,jj) = normalize(tof(:,jj));
end


factordata = factorXem;


%% Baseline

close all

Idata = [ factordata(:,1) factordata(:,2) factordata(:,1).*factordata(:,2) ]; %interaction variables

hominterc = 0; % if = 1  inclusion of heterogeneous intercept

demean = 1; % if = 1 fixed effects

hor = 4; %lenght of IRF

lag = 3; %number of lags of X

parametric = 0; % = 0 to bootstraped CI, =1 for parametrics

restr =[ 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0]; % Declare with Xs need to be interacted [=1 interacted, =0 not interacted]

centered = 0; % = 0 median, = 1 mean

pct = 10; % = 10 90% CI

block = 3; %size of the block bootstrap

shock = [ 1 2 3 4 5 6 7 8 ]; % Variables that compose the plan

shock2 = [ 1 2 ]; % how many shocks

rep = 1000; % number of bootstrap replications

nint=size(Idata,2);   % number of interaction terms
neq=size(y,2);   % number of ebndogenous variables
nvar=size(x,2); %number of explanatori variables
nnlin=sum(restr); %number of interactions


%% Calculate the columns to suppress in the various lags to include only the plans announced at t-i implemented at t

toc = [ ]; % columns to suppress in the first lag
tc = isempty(toc);
if tc == 0
supp = [];
futsim = [];
for j = 1:size(toc,1)
for i = 1:size(toc,2)
    pos = toc(j,i);
    if pos == 0
        supph = zeros(1,1+nint);
    else
        suma = 1:nint+1;
        supph = suma+((pos-1)*(1+nint)+j*(nvar-nnlin + nnlin*(1+nint)));
    end
    supp = [supp supph];
    clear pos mul supph
end
end
futsim = supp(find(supp~=0));
else
futsim=[];
end


clear toc tc supph mul pos util futsimh posi

%% values of the percentiles
   
values=[prctile(Idata(:,1),80) prctile(Idata(:,2),50) prctile(Idata(:,3),50) 
        prctile(Idata(:,1),20) prctile(Idata(:,2),50) prctile(Idata(:,3),50) 
   	    prctile(Idata(:,1),50) prctile(Idata(:,2),80) prctile(Idata(:,3),50) 
	    prctile(Idata(:,1),50) prctile(Idata(:,2),20) prctile(Idata(:,3),50) ]; % Defining values at which interaction terms are to be evaluated for the IRF
    values=values';

%% Estimation and bootstrap

for gamma = 1:hor+1
    
    [beta,sterr,errors,ysel,Xsel]=estimatePLP(y,x,I,lag,Idata,restr,demean,gamma,cont,sum(restr));

%% Bootstraped CI

[BETAMAT ERRORMAT]=IPLPboot(parametric,Idata,ysel,Xsel,errors,beta,rep,I,lag,block);

ERRORMATG(:,gamma) = ERRORMAT;
BETAMATG(:,gamma) = BETAMAT;

end

%% Simulation

for boot = 1:rep

[DIFF,SIMUL] = simulationLP(BETAMATG(boot,:),elast,gamma,values,neq,lag,x,restr,shock2,numatt,futsim);


MATSIMUL(:,:,:,:,boot) = SIMUL(:,:,:,:);
MATDIFF(:,:,:,:,boot) = DIFF(:,:,:,:);

end

%% Generate Median and CI
[MEDIAN, STD] = distr(MATSIMUL,size(values,2),neq,shock2,centered,pct,gamma);
[DIFFMED, DIFFSTD] = distr(MATDIFF,size(values,2)/2,neq,shock2,centered,pct,gamma);


%% Generating Pictures

close all force

dopicsmain(DIFFMED,DIFFSTD,MEDIAN,STD,shock2,values,8,ynms)



%% END AUTOMATIC PART