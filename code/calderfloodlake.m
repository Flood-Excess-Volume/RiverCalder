clear;
%%

% values as of 29/5/19 for RRA review

FEV = 1.65;  % M cubic metres
depth = 2; % m
Lside = sqrt(FEV*10^6/depth);
xmax = 950;
Vr = [0.44 0.88];   % Reservoir:  M cubic metres 100 percent
Vnfm = [0.07 0.14]; % Leaky dams: M cubic metres 100 percent
% Vb = [0.025*FEV 0.05*FEV]; % Trees: M cubic metres 100 percent
Vb = [0.03*FEV 0.06*FEV]; % Trees: M cubic metres 100 percent
Lxv = [Vr, Vnfm, Vb]/FEV;
nfm = mean(Lxv(3:4));
res = mean(Lxv(1:2));
tree = mean(Lxv(5:6));

lower_fev = [Lxv(1) Lxv(3) Lxv(5)];
upper_fev = [Lxv(2) Lxv(4) Lxv(6)];

lower_fev_sum = cumsum(lower_fev);
upper_fev_sum = cumsum(upper_fev);

alph = 0.4;
fs = 16;
%%

figure(101);
Lx = Lside;
Ly = Lside;

xr = Lx*[0,res,res,0,0];
yr = Ly*[0,0,1,1,0];
plot(xr,yr,'-b','linewidth',2); hold on;
patch(xr,yr,[0 0 0.99]); alpha(alph); hold on;

xn = Lx*[res,res+nfm,res+nfm,res,res];
yn = Ly*[0,0,1,1,0];
plot(xn,yn,'-r','linewidth',2); hold on;
patch(xn,yn,[0.99 0 0]); alpha(alph); hold on;

xt = Lx*[res+nfm,res+nfm+tree,res+nfm+tree, res+nfm,res+nfm];
yt = Ly*[0,0,1,1,0];
plot(xt,yt,'-g','linewidth',2); hold on;
patch(xt,yt,[0 0.99 0]); alpha(alph); hold on;

arr = annotation('doublearrow');
arr.Parent = gca;
arr.X = Lx*[Lxv(1) Lxv(2)];
arr.Y = [0.8*Ly 0.8*Ly];
arr.Color = 'blue';
arr.Head1Style = 'vback3';
arr.Head2Style = 'vback3';
arr.Head1Length = 5;
arr.Head2Length = 5;
arr.LineWidth = 2; 

arn = annotation('doublearrow');
arn.Parent = gca;
arn.X = Lx*[res+Lxv(3) res+Lxv(4)];
arn.Y = [0.6*Ly 0.6*Ly];
arn.Color = 'red';
arn.Head1Style = 'vback3';
arn.Head2Style = 'vback3';
arn.Head1Length = 5;
arn.Head2Length = 5;
arn.LineWidth = 2; 

art = annotation('doublearrow');
art.Parent = gca;
art.X = Lx*[res+nfm+Lxv(5) res++nfm+Lxv(6)];
art.Y = [0.4*Ly 0.4*Ly];
art.Color = 'green';
art.Head1Style = 'vback3';
art.Head2Style = 'vback3';
art.Head1Length = 5;
art.Head2Length = 5;
art.LineWidth = 2; 

arm = annotation('doublearrow');
arm.Parent = gca;
arm.X = Lx*[0.3344 0.6678];
arm.Y = [0.2*Ly 0.2*Ly];
arm.Color = 'black';
arm.Head1Style = 'vback3';
arm.Head2Style = 'vback3';
arm.Head1Length = 5;
arm.Head2Length = 5;
arm.LineWidth = 2; 

xlabel('Sidelength (m)','fontsize',fs);
ylabel('Sidelength (m)','fontsize',fs);
axis([0 Lside 0 Lside 0 depth]);
title(sprintf('Flood-excess lake: FEV $$\\approx %d$$m$$^2$$ x $$2$$m $$ \\approx %.3f$$Mm$$^3$$',...
        round(Lside,0),FEV),'Interpreter','latex','fontsize',20);
% box on


%%
% gtext('53.3%','fontsize',20);
% gtext('8.48%','fontsize',20);
% gtext('5%','fontsize',20);
% gtext('26.7%','fontsize',20);
% gtext('4.24%','fontsize',20);
% gtext('2.5%','fontsize',20);
% gtext('reservoirs','fontsize',20);
% gtext('NFM','fontsize',20);
% gtext('trees','fontsize',24);
% gtext('mean mitigation','fontsize',20);
% gtext('50.11%','fontsize',20);
% gtext('66.81%','fontsize',20);
% gtext('33.41%','fontsize',20);
% gtext('$\pounds$30M at $\pounds$[0.56,1.13]M$/1\%$','Interpreter','latex','fontsize',20);
% gtext('$\pounds$5.38M at $\pounds$[0.63,1.27]M$/1\%$','Interpreter','latex','fontsize',20);
% gtext('$\pounds$5M at $\pounds$[1,2]M$/1\%$','Interpreter','latex','fontsize',20);
% gtext('$\pounds$40.38M at $\pounds$0.804M$/1\%$','Interpreter','latex','fontsize',20);

%%
figure(103);
ax1 = axes('Position',[0.11 0.11 0.75 0.75]);

%reservoir
xr = [0,lower_fev_sum(1),upper_fev_sum(1),0,0];
yr = [0,0,1,1,0];
% plot(xr,yr,'-b','linewidth',2); hold on;
patch(xr,yr,[0 0 0.75]); 
alpha(alph); 
hold on;

%nfm
xn = [lower_fev(1),lower_fev_sum(2),upper_fev_sum(2),upper_fev_sum(1),lower_fev_sum(1)];
yn = [0,0,1,1,0];
% plot(xn,yn,'-r','linewidth',2); hold on;
patch(xn,yn,[0.25 0 0]); 
alpha(alph); 
hold on;

%trees
xt = [lower_fev_sum(2),lower_fev_sum(3),upper_fev_sum(3),upper_fev_sum(2),lower_fev_sum(2)];
yt = [0,0,1,1,0];
% plot(xt,yt,'-g','linewidth',2); hold on;
patch(xt,yt,[0 0.75 0]); 
alpha(alph); 
hold on;

%total
xtot = [0,lower_fev_sum(3),upper_fev_sum(3),0,0];
ytot = [0,0,1,1,0];
plot(xtot,ytot,'-k','linewidth',3); hold on;


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

ax1.XTick =  [0 lower_fev_sum 1];
ax1.YTick = [ ];
ax1.XTickLabel = round(100*[0 lower_fev_sum 1],1);
ax1.XLim = [0 1];
ax1.YLim = [0 1];
ax1.XTickLabelRotation = 60;
ax1.XLabel.String = '% of FEV mitigated';
ax1.FontSize = fs-2;
ax1.TickLabelInterpreter = 'latex';
ax1.TickDir ='both';


ax2.XTick = [0 upper_fev_sum 1];
ax2.YTick = [ ];
ax2.XTickLabel = round(100*[0 upper_fev_sum 1],1);
ax2.XLim = [0 1];
ax2.YLim = [0 1];
ax2.XTickLabelRotation = 60;
% ax2.XLabel.String = '% of FEV mitigated';
ax2.FontSize = fs-2;
ax2.TickLabelInterpreter = 'latex';
ax2.TickDir ='both';


pos = get(gca,'Position');

arr = annotation('doublearrow');
arr.X = 1.025*[pos(1)+pos(3) pos(1)+pos(3)];
arr.Y = [pos(2) pos(2)+pos(4)];
arr.Color = 'black';
arr.Head1Style = 'plain';
arr.Head2Style = 'plain';
arr.Head1Length = 5;
arr.Head2Length = 5;
arr.LineWidth = 2; 

best = text(1.07,1,'Best case','HorizontalAlignment', 'right','Fontsize',fs);
set(best,'Rotation',90);
worst = text(1.07,0,'Worst case','HorizontalAlignment', 'left','Fontsize',fs);
set(worst,'Rotation',90);

tres = annotation('textbox','String',{'RESERVOIRS','Mean FEV: $40\%$','Cost: $\pounds 30$M','Value: $\pounds$[0.56, 1.13]M$/1\%$'},...
    'Interpreter','latex');
tres.Color = [0 0 0.75];
tres.LineStyle = 'none';
tres.FontSize = fs;
tres.Position = [1.025*pos(1),pos(2)+0.6*pos(4),0.2,0.2];
tres.FitBoxToText = 'on';

tnfm = annotation('textbox','String',{'NFM: LEAKY DAMS','Mean FEV: $6.36\%$','Cost: $\pounds 5.38$M','Value: $\pounds$[0.63,1.27]M$/1\%$'},...
    'Interpreter','latex');
tnfm.Color = [0.25 0 0];
tnfm.LineStyle = 'none';
tnfm.FontSize = fs;
tnfm.Position = [0.55,pos(2)+0.3*pos(4),0.2,0.2];
tnfm.FitBoxToText = 'on';

% anfm = annotation('arrow');
% anfm.Parent = gca;
% anfm.X = Lx*[Lxv(1) Lxv(2)];
% anfm.Y = [0.8*Ly 0.8*Ly];
% anfm.Color = [0.75 0 0];
% anfm.HeadStyle = 'vback3';
% anfm.HeadLength = 5;
% anfm.LineWidth = 2; 

ttre = annotation('textbox','String',{'NFM: TREES','Mean FEV: $4.5\%$','Cost: $\pounds 6$M','Value: $\pounds$[1,2]M$/1\%$'},...
    'Interpreter','latex');
ttre.Color = [0 0.5 0];
ttre.LineStyle = 'none';
ttre.FontSize = fs;
ttre.Position = [0.6,pos(2)+0.55*pos(4),0.2,0.2];
ttre.FitBoxToText = 'on';

ttot = annotation('textbox','String',{'TOTAL','Mean FEV: 50.86\%','Cost: $\pounds 41.38$M','Value: $\pounds$[0.61,1.22]M$/1\%$'},...
    'Interpreter','latex');
ttot.Color = [0 0 0];
ttot.LineWidth = 3;
ttot.FontSize = fs;
ttot.Position = [0.5,pos(2),0.2,0.2];
ttot.FitBoxToText = 'on';


title(sprintf('Flood-excess lake: FEV $$\\approx (%d^2 \\times 2)$$m$$^3 \\approx %.3f$$Mm$$^3$$',...
        round(Lside,0),FEV),'Interpreter','latex','fontsize',20);


%% to export: x by y cm?
