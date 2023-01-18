function impcompareoutout(MA,STD,line,nvar,numpic,shock,Names)

%MA is period x Number of IRs matrix and stands for the IRS that should be
%comparded in a row 
%STD is a period x 2 matrix Numbber of IRs   and contains the
%confidence intervals 

%IRS are supposed to be for a single variable and a single schock from different
%estimations.

%variable is a string   

% GENERATE IRF 
%-------------------------------------------------------------------

% Number of periods and IRS
[nHorizon noIR]=size(MA);
noIR = noIR-(noIR/3);
nHorizon = nHorizon-1;
count = 0;
% Number of periods

Highcolor     =[0 0 1];      HighbandFillColor  =[.2 .6 1];
Lowcolor      =[.9 .1 .1];   LowbandFillColor   =[.9 .6 .4];
Diffcolor     =[0 0.7 0.2];    DiffbandFillColor  =[.1 .6 .1];

max1=max(max(max(max(STD(:,:,:)))),max(max(max(STD(:,:,:)))));
min1=min(min(min(min(STD(:,:,:)))),min(min(min(STD(:,:,:)))));
maximum=max1*sign(max1)*1.1+max1*(1-sign(max1))*0.9;
minimum=min1*(1-sign(min1))*1.1+min1*(sign(min1))*0.9;
% PLOT IRF
%-------------------------------------------------------------------
for i=1:noIR
    if mod(i,2) == 1
    intercept = zeros(size(MA));
    subplot(nvar/2,noIR*2,(line+numpic-1))
    hold on
    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
    [STD(1:nHorizon+1,1,i+count)' fliplr(STD(1:nHorizon+1,2,i+count)')],...
    HighbandFillColor,'EdgeColor',HighbandFillColor,'LineStyle','-.','facealpha',.3);

    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
    [STD(1:nHorizon+1,1,i+1+count)' fliplr(STD(1:nHorizon+1,2,i+1+count)')],...
    LowbandFillColor,'EdgeColor',LowbandFillColor,'LineStyle','--','facealpha',.3);

    %irfs
    p0=plot(0:nHorizon,MA(:,i+count), '-.','LineWidth',1.5,'color',Highcolor);
    p1=plot(0:nHorizon,MA(:,i+1+count),  '--','LineWidth',1.5,'color',Lowcolor);

    %zero line
    plot(0:nHorizon,zeros(size(1:nHorizon+1)),'k')

    hold off; axis tight
    ylim([(min1+0.1*min1) (max1+0.1*max1)])        
    xlim([0 nHorizon]);
    set(gca,'XTick',0:1:nHorizon,'XTickLabel',cellstr(num2str((0:1:nHorizon)')),'Layer','top','FontSize',8)
    
    ax.YTickMode = 'manual';
    

        ax = gca;
        ticks=ax.YTick;%retrieve current ticks
        ax.YTickLabel=round(ticks,2)*100;%multiply  

               
    ylabel(Names.Properties.VariableNames(numpic),'FontSize',10)

        [lh,icons]=legend([p0 p1],{['High']; ['Low']},'FontSize',8,'Location','best', 'Orientation', 'horizontal','NumColumnsMode','manual');
        lh.NumColumns = 2;
        set(lh,'box','off')
        p1 = icons(1).Position;
        icons(1).Position = [0.3 p1(2) p1(3)];
        icons(3).XData = [0.2 0.27];
        p2 = icons(2).Position;
        icons(2).Position = [0.53 p2(2) p2(3)];
        icons(5).XData = [0.43 0.5];
    else
     count=count+1;
    intercept = zeros(size(MA));
    subplot(nvar/2,noIR*2,(line+numpic))
    hold on
    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
    [STD(1:nHorizon+1,1,i+count)' fliplr(STD(1:nHorizon+1,2,i+count)')],...
    DiffbandFillColor,'EdgeColor',DiffbandFillColor,'LineStyle','-.','facealpha',.3);

    %irfs
    p0=plot(0:nHorizon,MA(:,i+count), '-.','LineWidth',1.5,'color',Diffcolor);

    %zero line
    plot(0:nHorizon,zeros(size(1:nHorizon+1)),'k')

    hold off; axis tight
    ylim([(min1+0.1*min1) (max1+0.1*max1)])        
    xlim([0 nHorizon]);
    set(gca,'XTick',0:1:nHorizon,'XTickLabel',cellstr(num2str((0:1:nHorizon)')),'Layer','top','FontSize',8)
    
    ax.YTickMode = 'manual';
    

    ax = gca;
    ticks=ax.YTick;%retrieve current ticks
    ax.YTickLabel=round(ticks,2)*100;%multiply  

    end
end
end
   
