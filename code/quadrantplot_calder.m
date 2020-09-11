%%% Script carries out the following:
%
% > loads and extracts relevant flow-data for Mytholmroyd
% > calculates rating curve from empirical formula
% > computes FEV and duration etc.
% > plots data using 2 subroutines
% <plot3panel.m> 
% <plotFEVhT.m>
% <plot_h_year.m>
% <plot_ratingcurve.m>

% TK, August 2018 - adapted from OB's <flowdatafloodsmy.m>

clear;

%% DATA: load + extract
% Flood data analysis:
%
% Mytholmroyd Calder River excel files: 22 to 32181
%
% Next is same for Aire and Calder:
armleyexc = [22, 32181];
nlenarmleyfi = armleyexc(2)-armleyexc(1)+1;
armleymeimaartdayspmo = [31,30,31,31,30,31,30,31,31,28,31];
day15min = 24*4;
daysfromay = 31+30+31+31+30+31+30.0-1.0;
no15minarmley = (sum(armleymeimaartdayspmo))*day15min;
fprintf('number 15min intervals May 2015 to March 2016 Armley/Calder flood gauge: %g %g\n',no15minarmley,nlenarmleyfi);
ndecarmley(1:2) = [(sum(armleymeimaartdayspmo(1:7))+24), (sum(armleymeimaartdayspmo(1:7))+29)]*day15min;
%

%%
fid = fopen('Mytholmroyd F1204 Flow 15 min May 15 to Mar 16.csv'); %
dataflow = textscan(fid, '%s %f %s %s %s %s %s %s %s', 'Delimiter', ',', 'HeaderLines', 21);
fclose(fid);

fid2 = fopen('Mytholmroyd F1204 Stage 15 min May 15 to Mar 16.csv'); %
datastage = textscan(fid2, '%s %f %s %f %f %s %f %s %s %s %s', 'Delimiter', ',', 'HeaderLines', 21);
fclose(fid2);

% extract Q data        
timestampQ = datetime(dataflow{1},'InputFormat','dd/MM/yyyy HH:mm:SS');
nq = length(dataflow{2});
timeQ = zeros(1,nq);
dischargeQ = zeros(1,nq);
for ii = 1:nq
    timeQ(ii) = 0.25*ii/24; % Times 0.25 gives hours, 1/24 gives day
    dischargeQ(ii) = dataflow{2}(ii);
end

% extract h data
timestamph = datetime(datastage{1},'InputFormat','dd/MM/yyyy HH:mm');
nh = length(datastage{2});
timeh = zeros(1,nh);
stage = zeros(1,nh);
for ii = 1:nh
    timeh(ii) = 0.25*ii/24;
    stage(ii) = datastage{2}(ii);
end

%% Rating curve
arml = [0., 2.107, 3.088]; % lower stage limit
armu = [2.107, 3.088, 5.8]; % upper stage limit
se = [0.849, 0.136, 0.136];

%rc coeffs
Crc = [8.459, 21.5, 2.086];
brc = [2.239,1.37,2.515];
arc = [0.342,0.826,-0.856];

Ns = 10000;
[armlevc,ic] = max(stage);
rcstage = arml(1):(armlevc-arml(1))/Ns:armlevc;
[armlevc,ic] = max(rcstage);

rcdischarge = rcstage; % Q for rating curve
rcdischargeL = rcstage;
rcdischargeU = rcstage;

discharge = stage;
dischargeL = stage;
dischargeU = stage;

%%
% for ii = 1:na
%     if (stage(ii) < armu(1)) & (stage(ii) > arml(1))
%         discharge(ii) = Crc(1)*(stage(ii)-arc(1))^brc(1);
%     else
%         if (stage(ii) < armu(2))
%             discharge(ii) = Crc(2)*(stage(ii)-arc(2))^brc(2);
%         else
%             discharge(ii) = Crc(3)*(stage(ii)-arc(3))^brc(3);
%         end
%     end
% end
% %
% for ii = 1:Ns+1
%     if (rcstage(ii) < armu(1)) & (rcstage(ii) > arml(1))
%         rcdischarge(ii) = real(Crc(1)*(rcstage(ii)-arc(1))^brc(1));
%     else
%         if (stage(ii) < armu(2))
%             rcdischarge(ii) = Crc(2)*(rcstage(ii)-arc(2))^brc(2);
%         else
%             rcdischarge(ii) = Crc(3)*(rcstage(ii)-arc(3))^brc(3);
%         end
%     end
% end

% % Q = Q(h)
for ii = 1:Ns+1
    
    if (rcstage(ii) <= armu(1)) && (rcstage(ii) >= arml(1))
        
        rcdischarge(ii) = real(Crc(1)*(rcstage(ii)-arc(1))^brc(1));
        rcdischargeL(ii) = (1.0-se(1))*rcdischarge(ii); % -SE
        rcdischargeU(ii) = (1.0+se(1))*rcdischarge(ii); % +SE
        
    elseif (rcstage(ii) <= armu(2)) && (rcstage(ii) > arml(2))
        
        rcdischarge(ii) = Crc(2)*(rcstage(ii)-arc(2))^brc(2);
        rcdischargeL(ii) = (1.0-se(2))*rcdischarge(ii); % -SE
        rcdischargeU(ii) = (1.0+se(2))*rcdischarge(ii); % +SE
        
    elseif (rcstage(ii) > armu(2))
        
        rcdischarge(ii) = Crc(3)*(rcstage(ii)-arc(3))^brc(3);
        rcdischargeL(ii) = (1.0-se(3))*rcdischarge(ii); % -SE
        rcdischargeU(ii) = (1.0+se(3))*rcdischarge(ii); % +SE
        
    end
    
end

% Q = Q(h(t)) = Q(t)
for ii = 1:nh
    
    if (stage(ii) <= armu(1)) && (stage(ii) >= arml(1))
        
        discharge(ii) = Crc(1)*(stage(ii)-arc(1))^brc(1);
        dischargeL(ii) = (1.0-se(1))*discharge(ii); % -SE
        dischargeU(ii) = (1.0+se(1))*discharge(ii); % +SE
        
    elseif (stage(ii) <= armu(2)) && (stage(ii) > arml(2))
        
        discharge(ii) = Crc(2)*(stage(ii)-arc(2))^brc(2);
        dischargeL(ii) = (1.0-se(2))*discharge(ii); % -SE
        dischargeU(ii) = (1.0+se(2))*discharge(ii); % +SE
        
    elseif (stage(ii) > armu(2))
        
        discharge(ii) = Crc(3)*(stage(ii)-arc(3))^brc(3);
        dischargeL(ii) = (1.0-se(3))*discharge(ii); % -SE
        dischargeU(ii) = (1.0+se(3))*discharge(ii); % +SE
        
    end
    
end

%% calculate FEV etc

armin = 4.0;
armax = 5.7;
Na = 20;
da = 0.1; % (armax-armin)/Na;
armstep = armin:da:armax;
Na = size(armstep,2);
nle = size(stage(ndecarmley(1):ndecarmley(2)),2);
arnle = zeros(1,nle);

for nna = 1:Na
    armcrit = armstep(nna); % m
    [harmax,iarmax] = max(stage(ndecarmley(1):ndecarmley(2)));
    [harmup,iarmup] = min(abs(stage(ndecarmley(1):ndecarmley(1)+iarmax)-armcrit));
    [harmdo,iarmdo] = min(abs(stage(ndecarmley(1)+iarmax:ndecarmley(2))-armcrit));
    iarmdo = iarmax+iarmdo;
    %
    %
    flowarmcrit = discharge(ndecarmley(1)+iarmup-1);
    armexcessvol = 15*60*sum( discharge(ndecarmley(1)+iarmup-1:ndecarmley(1)+iarmdo-1)-flowarmcrit ); % 15min
    armexcvol(nna) = armexcessvol;
    %
    flowarmcrit = Crc(3)*(armcrit-arc(3))^brc(3);
    posarm = max(discharge(ndecarmley(1):ndecarmley(2))-flowarmcrit,arnle);
    armexcessvol = 15*60*sum(posarm);
    flowarmcritvec(nna) = Crc(3)*(armcrit-arc(3))^brc(3);
    Tfvec(nna) = 15*60*nnz(posarm);
    armexcvol(nna) = armexcessvol;
end

%%
armcrit = 4.5; % hT [m]

QT = Crc(3)*(armcrit-arc(3))^brc(3);
QTminus = QT*(1-se(3));
QTplus = QT*(1+se(3));

posarm = max(discharge(ndecarmley(1):ndecarmley(2))-QT,arnle);
posarmL = max(dischargeL(ndecarmley(1):ndecarmley(2))-QTplus,arnle);
posarmU = max(dischargeU(ndecarmley(1):ndecarmley(2))-QTminus,arnle);

% FEVs
armexcessvol = 15*60*sum(posarm);
armexcessvolL = 15*60*sum(posarmL); 
armexcessvolU = 15*60*sum(posarmU);

% flood durations
Tf = 15*60*nnz(posarm);
TfL = 15*60*nnz(posarmL);
TfU = 15*60*nnz(posarmU);

% flowarmcrit = Crc(3)*(armcrit-arc(3))^brc(3);
% posarm = max(discharge(ndecarmley(1):ndecarmley(2))-flowarmcrit,arnle);
% armexcessvol = 15*60*sum(posarm)
% Tf = 15*60*nnz(posarm)


Qm = QT+armexcessvol/Tf; 
hm = (Qm/Crc(3))^(1/brc(3))+arc(3);

%
[harmax,iarmax] = max(stage(ndecarmley(1):ndecarmley(2)));
[harmup,iarmup] = min(abs(stage(ndecarmley(1):ndecarmley(1)+iarmax)-armcrit));
[harmdo,iarmdo] = min(abs(stage(ndecarmley(1)+iarmax:ndecarmley(2))-armcrit));
iarmdo = iarmax+iarmdo;
%
% flowarmcrit = myflow2(ndecarmley(1)+iarmup);
% armexcessvol = 15*60*sum(myflow2(ndecarmley(1)+iarmup:ndecarmley(1)+iarmdo)-flowarmcrit); % 15min
armeymean = 0.5*(max(stage)+armcrit);
armexcessvolest = max(discharge)*(armeymean/max(stage))*(max(stage)-armcrit)/max(stage)*(iarmdo-iarmup)*0.25*3600;
armexcessvolest2 = max(discharge)*(armeymean/max(stage))*(max(stage)-armeymean)/max(stage)*(iarmdo-iarmup)*0.25*3600;

%


%% data for plotting
t = timeh(ndecarmley(1):ndecarmley(2))-daysfromay; %time
h = stage(ndecarmley(1):ndecarmley(2)); %h data
q = dischargeQ(ndecarmley(1):ndecarmley(2)); %q data

q2 = discharge(ndecarmley(1):ndecarmley(2));
q2L = dischargeL(ndecarmley(1):ndecarmley(2));
q2U = dischargeU(ndecarmley(1):ndecarmley(2));

t = t(1:4*day15min+1);
h = h(1:4*day15min+1);
q = q(1:4*day15min+1);
q2 = q2(1:4*day15min+1);
q2L = q2L(1:4*day15min+1);
q2U = q2U(1:4*day15min+1);

hrc = rcstage; % h rating curve
qrc = rcdischarge; % q rating curve

%% Tf as a function of QT for error calc
here = find(Tfvec == Tf);
here = here(2);
dTfdQT = (Tfvec(here+1) - Tfvec(here-1))/(flowarmcritvec(here+1) - flowarmcritvec(here-1));
x = linspace(flowarmcritvec(here-2),flowarmcritvec(here+2),11);
y = dTfdQT*(x - QT) + Tf;
plot(flowarmcritvec, Tfvec,'k'); hold on;
% plot([QT, QT],[Tf,Tfvec(end)],'k:');
% plot([flowarmcritvec(1),QT],[Tf, Tf],'k:');
plot(x,y,'r','linewidth',2);
plot(QT,Tf,'or','linewidth',2);
text(1.01*QT,1.01*Tf,'$(Q_T, T_f)$','fontsize',16, 'HorizontalAlignment', 'left','Interpreter','latex');
% text(0.9*QT,0.9*Tf,'$$\frac{\partial Q_T}{\partial T_f)$$','fontsize',16, 'HorizontalAlignment', 'left','Interpreter','latex');
xlabel('$Q_T$ [m$^3$/s]','fontsize',16,'Interpreter','latex');
ylabel('$T_f$ [s]','fontsize',16,'Interpreter','latex');

%% error tests
sd = max(se);
sd = 0.136;
dt = 15*60;
qk = q(q>QT);
q2k = q2(q2>QT);

Ve = armexcessvol;
dVedTf = Ve/Tf;

% both ignoring Tf contributions
errVe2 = dt^2*sd^2*sum(qk.^2) + Tf^2*sd^2*QT^2; %applying rc N+1 times (Qk for k=1,...,N and QT)
errVm2 = Tf^2*sd^2*Qm^2 + Tf^2*sd^2*QT^2; % applying rc 1+1 times (Qm and QT)

% both including Tf contributions
errVe2 = dt^2*sd^2*sum(qk.^2) + (Tf^2 + dVedTf^2*dTfdQT^2)*sd^2*QT^2; %applying rc N+1 times (Qk for k=1,...,N and QT)
errVm2 = Tf^2*sd^2*Qm^2 + (Tf^2 + dVedTf^2*dTfdQT^2)*sd^2*QT^2; % applying rc 1+1 times (Qm and QT)

errVe = sqrt(errVe2);
errVm = sqrt(errVm2);

errVefrac = errVe/Ve
errVmfrac = errVm/Ve

%% plotting routines

% 3 panel with h(t), Q(h), Q(t)
plot3panel;
plot3panelerr;

% FEV and sidelength as a function of threshold
plotFEVhT;

% h(t) for whole year
plot_h_year;

% rating curve and discharge with errors
plot_ratingcurve;