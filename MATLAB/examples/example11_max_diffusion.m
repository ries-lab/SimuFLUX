%% moving, bleaching fluorophore and tracking with Abberior sequence
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path

savepath='/Users/ries/Library/CloudStorage/OneDrive-Personal/Projects/SimulFLUX/Figures/exports for figures/diffusion';
%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"SimSequencefile")
    sim=SimSequencefile;
end
fname='Tracking_2D.json';
fname2='PSFvectorial2D.json'; %use a PSF that is defined via a json file
sim.loadsequence(fname,fname2);
sim.sequence.locLimit=100;
sim.makepatterns;

defaultsequence=sim.sequence;
defaultdeadtimes=sim.deadtimes;
defaultestimators=sim.estimators;

posfig=[562,541,1450,410];

%conditions to change
% sim.sequence.damping=1;
%iterative estimator:
% sim.estimators(4).function='est_qLSQiter2D';
% sim.estimators(4).par{3}=false; % no ccr check
% sim.sequence.Itr(4).phtLimit=20;
% sim.sequence.Itr(4).patDwellTime=0.1e-3;
% sim.sequence.Itr(4).patGeoFactor=L/360;
 % L=itrs(itr).patGeoFactor*360; %nm;
% sim.deadtimes.point=0; sim.deadtimes.estimator=0;

laserpower=1;
updatetime=0.0025; %ms
repetitions=20;
maxerr=100;

%% plot example track
figure(256); clf
D=3.;
sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
fl=FlMoveBleach;
fl.photonbudget=inf;
fl.brightness=1000*laserpower;
fl.makediffusion(D,updatetime)
sim.fluorophores=fl;
out=sim.runSequence(repetitions=1,resetfluorophores=true);
err=sqrt((out.loc.xnm-out.loc.xfl1).^2+(out.loc.ynm-out.loc.yfl1).^2);
indlost=find(err<maxerr, 1,'last');
subplot(2,2,1)
hold off; plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl1)%;hold on; plot(out.loc.loccounter,out.loc.xgalvo)
xlabel('time (localizations)'); ylabel('x position (nm)'); 
hold on; plot(out.loc.loccounter, out.loc.ynm);hold on; plot(out.loc.loccounter,out.loc.yfl1)%;hold on; plot(out.loc.loccounter,out.loc.ygalvo)
plot([indlost,indlost],[min(out.loc.ynm),max(out.loc.ynm)])
xlabel('time (localizations)'); ylabel('x position (nm)'); 
legend('x estimated position','x fluorophore position','y estimated position','y fluorophore position')

subplot(2,2,2)
hold off; plot(out.loc.xnm(1:end-1), out.loc.ynm(1:end-1));
tfl=0:min(diff(out.loc.time)):out.loc.time(end);

hold on; plot(interp1(fl.pos(:,1),fl.pos(:,2),tfl),interp1(fl.pos(:,1),fl.pos(:,3),tfl))
xlabel('x position (nm)'); ylabel('y position (nm)'); 
axis equal


%% investigate what happens if fluorophore gets lost
D=3.;
numtracks=300;
lenwin=10;
averagejump=zeros(lenwin,1);
averageerr=zeros(lenwin,1);
nj=0;
for k=1:numtracks
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
    fl=FlMoveBleach;
    fl.photonbudget=inf;
    fl.brightness=1000*laserpower;
    fl.makediffusion(D,updatetime)
    sim.fluorophores=fl;
    out=sim.runSequence(repetitions=1,resetfluorophores=true);
    
    err=sqrt((out.loc.xnm-out.loc.xfl1).^2+(out.loc.ynm-out.loc.yfl1).^2);
    jx=vertcat(diff(out.loc.xfl1),0);
    jy=vertcat(diff(out.loc.yfl1),0);
    % jx=smoothdata(jx,'movmean',5);
    % jy=smoothdata(jy,'movmean',5);
    jump=sqrt(jx.^2+jy.^2);
    
    indlost=find(err<maxerr, 1,'last');
    try
        averagejump=averagejump+jump(indlost-lenwin+4:indlost+3);
        averageerr=averageerr+err(indlost-lenwin+4:indlost+3);
        nj=nj+1;
    end
end

subplot(2,2,3)
hold off
nx=(1:length(averageerr))'-7;
plot(nx,averagejump/nj); hold on; plot(nx,averageerr/nj);
xlabel('time (localization)')
ylabel('jump, localization error (nm)')
legend('average jump', 'average error')

%% MSD analysis
tlast=0;
D=3.;
% plot example track
msdwin=10;msdwinr=1000;
msdmf=zeros(msdwin,1);msdfl=msdmf;msdflr=zeros(msdwinr,1);
imsd=0; ir=0;
repetitionsmsd=100;
dtmsd=[];dtmsdr=[];
for k=1:repetitionsmsd
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
    fl=FlMoveBleach;
    fl.photonbudget=inf;
    fl.brightness=1000*laserpower;
    fl.makediffusion(D,updatetime)
    sim.fluorophores=fl;
    out=sim.runSequence(repetitions=1,resetfluorophores=true);
    err=sqrt((out.loc.xnm-out.loc.xfl1).^2+(out.loc.ynm-out.loc.yfl1).^2);
    indlost=find(err<maxerr, 1,'last');
    
    [msdhfr,dtmsdr]=getmsd(fl.pos(:,2),fl.pos(:,3),fl.pos(:,1),msdwinr,fl.pos(2,1)-fl.pos(1,1));
    msdflr=msdflr+msdhfr;
    ir=ir+1;
    
    if indlost>50
        time=out.loc.time(1:indlost);
        if isempty(dtmsd)
            dtmsd=min(diff(time));
        end
        xfl=interp1(fl.pos(:,1),fl.pos(:,2),time,"nearest");
        yfl=interp1(fl.pos(:,1),fl.pos(:,3),time,"nearest");
        % [msdhf,dt]=getmsd(out.loc.xfl1(1:indlost),out.loc.yfl1(1:indlost),time,msdwin);
         [msdhf,dtmsd]=getmsd(xfl,yfl,time,msdwin,dtmsd);
        
        [msdhm,dtmsd]=getmsd(out.loc.xnm(1:indlost),out.loc.ynm(1:indlost),time,msdwin,dtmsd);
        if ~(any(isnan(msdhm)|any(isnan(msdfl))))
            msdfl=msdfl+msdhf;
            msdmf=msdmf+msdhm;
            imsd=imsd+1;
        end
    end
end

msdfl=msdfl/imsd;
msdmf=msdmf/imsd;
msdflr=msdflr/ir;

%%
subplot(2,2,4);
hold off
[Dflr,offflr]=msdfit(msdflr,dtmsdr,'m');
[Dfl,offfl]=msdfit(msdfl,dtmsd,'r');
[Dmf,offmf,msdt]=msdfit(msdmf,dtmsd,'b');
xlim([0,msdt(end)/1e3])


legend("fluorophore all","D="+string(Dflr)+"µm2/s","fluorophore tracked","D="+string(Dfl)+"µm2/s","Minflux","D="+string(Dmf)+"µm2/s")

savefig(gcf,[savepath filesep 'tracks_abort_msd.fig'])

%% investigate maximum diffusion coefficient for various conditions
sim.sequence=defaultsequence; %reset
sim.makepatterns;
Dtest=[0 0.1 0.2 0.25:0.25:7];
Dtest=0:1;

laserpowers=[0.05 0.1 0.2:0.2:0.8 1:0.5:5];
clear Dmaxl efoD0l rmseDmaxl
f=figure(301);clf
for k=1:length(laserpowers)
    figure(f)
    [Dmaxl(k),efoD0l(k),rmseDmaxl(k)]=maxDiffusion(sim, laserpowers(k),repetitions, Dtest,string(laserpowers(k)))
end
title('laserpowers')
f.Position=posfig;
savefig(f,[savepath filesep 'laserpowers.fig'])
saveas(f,[savepath filesep 'laserpowers.eps'],'epsc')

figure(258); clf
subplot(3,4,1)
plot(laserpowers,Dmaxl)
xlabel('laserpower (a.u.)')
ylabel('max Diffusion D µm^2/s')
title('laserpower')
text(0,0.35,num2str(efoD0l','%2.0f'))

subplot(3,4,5)
plot(laserpowers,rmseDmaxl)
xlabel('laserpower (a.u.)')
ylabel('rmse at Dmax (nm)')
title('laserpower')

subplot(3,4,9)
plot(laserpowers,efoD0l)
ylabel('efo (kHz)')
xlabel('laser power (a.u.)')
title('laserpower')

subplot(3,4,10)
plot(efoD0l,rmseDmaxl)
xlabel('efo (kHz)')
ylabel('rmse at Dmax (nm)')
title('laserpower')
% text(0,0.35,num2str(efoDmaxl','%2.0f'))

%% redo with no dead times and faster pattern dwell times
sim.sequence=defaultsequence; %reset
sim.deadtimes.point=0; sim.deadtimes.estimator=0;
sim.sequence.Itr(4).patDwellTime=2e-5; %20 us
sim.makepatterns()
Dtestf=[0 0.1 0.2 0.3 0.5:0.5:20];
Dtestf=0:1
clear Dmaxlf efoD0lf rmseDmaxlf
f=figure(302); clf
for k=1:length(laserpowers)
    figure(f);
    [Dmaxlf(k),efoD0lf(k),rmseDmaxlf(k)]=maxDiffusion(sim, laserpowers(k),repetitions, Dtestf,string(laserpowers(k)))
end
sim.deadtimes=defaultdeadtimes; %back to default;

f.Position=posfig;
title('laserpower fast')
% legend(string(laserpowers'))
savefig(f,[savepath filesep 'laserpowersfast.fig'])
saveas(f,[savepath filesep 'laserpowersfast.eps'],'epsc')

figure(258)
subplot(3,4,1)
hold on
plot(laserpowers,Dmaxlf)
legend('standard','no deadtime')

subplot(3,4,5)
hold on
plot(laserpowers,rmseDmaxlf)
legend('standard','no deadtime')

subplot(3,4,9)
hold on
plot(laserpowers,efoD0lf)


subplot(3,4,10)
hold on
plot(efoD0lf,rmseDmaxlf)




%% L scan pattern size
sim.sequence=defaultsequence; %reset
clear DmaxL efoD0L rmseDmaxL
Ldefault=sim.sequence.Itr(4).patGeoFactor*360; % 100
Ls=[40 75 100 150 200 250];
laserpowerLnorm=laserpower*Ldefault^2./Ls.^2;
f=figure(303); clf
for k=1:length(Ls)
    figure(f);
    sim.sequence.Itr(4).patGeoFactor=Ls(k)/360;
    sim.makepatterns;
    
    [DmaxL(k),efoD0L(k),rmseDmaxL(k)]=maxDiffusion(sim, laserpowerLnorm(k),repetitions, Dtest,string(Ls(k)))
end
title('L')
subplot(1,3,2)
text(0,0.35,num2str(laserpowerLnorm','%2.2f'))
f.Position=posfig;
savefig(f,[savepath filesep 'Ls.fig'])
saveas(f,[savepath filesep 'Ls.eps'],'epsc')

figure(258)
subplot(3,4,2)
plot(Ls,DmaxL)
xlabel('Scan pattern size L (nm)')
ylabel('max Diffusion D µm^2/s')
title('L')

subplot(3,4,6)
plot(Ls,rmseDmaxL)
xlabel('Scan pattern size L (nm)')
ylabel('rmse at Dmax (nm)')
title('L')

%% photon limit
sim.sequence=defaultsequence; %reset
clear Dmaxpl efoD0pl rmseDmaxpl
photlim=[5 10 20 30 50 100 200];
f=figure(304);clf
for k=1:length(photlim)
    figure(f)
    sim.sequence.Itr(4).phtLimit=photlim(k);
    sim.makepatterns;
    [Dmaxpl(k),efoD0pl(k),rmseDmaxpl(k)]=maxDiffusion(sim, laserpower,repetitions, Dtest,string(photlim(k)))
end
title('photon limit')
f.Position=posfig;
savefig(f,[savepath filesep 'photlim.fig'])
saveas(f,[savepath filesep 'photlim.eps'],'epsc')

figure(258)
subplot(3,4,3)
plot(photlim,Dmaxpl)
xlabel('photon limit')
ylabel('max Diffusion D µm^2/s')
title('photon limit')

subplot(3,4,7)
plot(photlim,rmseDmaxpl)
xlabel('photon limit')
ylabel('rmse at Dmax (nm)')
title('photon limit')



%% dwell time
sim.sequence=defaultsequence; %reset
dwelltimes=[10 20 50 100 200 500 1000]; %us
clear Dmaxdt efoD0dt rmseDmaxdt
f=figure(305);clf
for k=1:length(dwelltimes)
    figure(f)
    sim.sequence.Itr(4).patDwellTime=dwelltimes(k)*1e-6;
    sim.makepatterns();
    [Dmaxdt(k),efoD0dt(k),rmseDmaxdt(k)]=maxDiffusion(sim, laserpower,repetitions, Dtest,string(dwelltimes(k)))
end
title('pattern dwell time')
f.Position=posfig;
savefig(f,[savepath filesep 'dwelltimes.fig'])
saveas(f,[savepath filesep 'dwelltimes.eps'],'epsc')

figure(258)
subplot(3,4,4)
plot(dwelltimes,Dmaxdt)
xlabel('pattern dwell time µs')
ylabel('max Diffusion D µm^2/s')
title('pattern dwell time')

subplot(3,4,8)
plot(dwelltimes,rmseDmaxdt)
xlabel('pattern dwell time µs')
ylabel('rmse at Dmax (nm)')
title('pattern dwell time')

%% redo with fast instrument
sim.sequence=defaultsequence; %reset
sim.deadtimes.point=0; sim.deadtimes.estimator=0;
f=figure(306);clf
clear Dmaxdtf efoD0dtf rmseDmaxdtf
for k=1:length(dwelltimes)
    figure(f)
    sim.sequence.Itr(4).patDwellTime=dwelltimes(k)*1e-6;
    sim.makepatterns();
    [Dmaxdtf(k),efoD0dtf(k),rmseDmaxdtf(k)]=maxDiffusion(sim, laserpower,repetitions, Dtest,string(dwelltimes(k)))
end
sim.deadtimes=defaultdeadtimes; %back to default;
title('pattern dwell time fast')
f.Position=posfig;
savefig(f,[savepath filesep 'dwelltimesfast.fig'])
saveas(f,[savepath filesep 'dwelltimesfast.eps'],'epsc')

figure(258)
subplot(3,4,4)
hold on
plot(dwelltimes,Dmaxdtf)
xlabel('pattern dwell time µs')
ylabel('max Diffusion D µm^2/s')
title('pattern dwell time')
legend('normal','fast')

subplot(3,4,8)
hold on
plot(dwelltimes,rmseDmaxdtf)
xlabel('pattern dwell time µs')
ylabel('rmse at Dmax (nm)')
title('pattern dwell time')



%% damping factor
sim.sequence=defaultsequence; %reset
dampingfs=[0 0.2 0.5 1 2 5 10];
f=figure(307);clf
clear Dmaxdamp efoD0damp rmseDmaxdamp
for k=1:length(dampingfs)
    figure(f)
    sim.sequence.damping=dampingfs(k);
    sim.makepatterns();
    [Dmaxdamp(k),efoD0damp(k),rmseDmaxdamp(k)]=maxDiffusion(sim, laserpower,repetitions, Dtest,string(dampingfs(k)))
end
title('damping factor')
f.Position=posfig;
savefig(f,[savepath filesep 'dampingfactor.fig'])
saveas(f,[savepath filesep 'dampingfactor.eps'],'epsc')

figure(258)
subplot(3,4,11)
hold on
plot(dampingfs,Dmaxdamp)
xlabel('pattern dwell time µs')
ylabel('max Diffusion D µm^2/s')
title('damping factor')

subplot(3,4,12)
hold on
plot(dampingfs,rmseDmaxdamp)
xlabel('pattern dwell time µs')
ylabel('rmse at Dmax (nm)')
title('damping factor')

savefig(gcf,[savepath filesep 'D_scan_vars.fig'])
%% compare different scan patterns
sim.sequence=defaultsequence; %reset
sim.makepatterns();
repetitions=100;
f=figure(310); clf
% reference
[Dmaxp6,efop6]=maxDiffusion(sim, laserpower,repetitions,[],"hexagon")

% 3 instead of 6 probing positions
sim.sequence.Itr(4).Mode.pattern='triangle';
sim.makepatterns;
[Dmaxp3,efop3]=maxDiffusion(sim, laserpower,repetitions,[],"triangle")

% 4 instead of 6 probing positions
sim.sequence.Itr(4).Mode.pattern='square';
sim.makepatterns;
[Dmaxp4,efop4]=maxDiffusion(sim, laserpower,repetitions,[],"square")

title("scan patterns")
f.Position=posfig;
savefig(f,[savepath filesep 'scanpatterns2.fig'])
saveas(f,[savepath filesep 'scanpatterns2.eps'],'epsc')

%%
f=figure(312); clf
repetitions=100;
% standard
sim.sequence=defaultsequence;
sim.makepatterns;
[Dmaxp6,efop6]=maxDiffusion(sim, laserpower,repetitions,[],"standard")

%optimized

Lm=150;Ldefault=102;
sim.sequence=defaultsequence;
sim.sequence.Itr(4).patDwellTime=5e-5; %50 us
sim.sequence.Itr(4).Mode.pattern='square';
sim.sequence.Itr(4).patGeoFactor=Lm/360;
laserpowerLnorm=laserpower*Ldefault^2./Lm.^2;
sim.sequence.Itr(4).phtLimit=10;
sim.makepatterns;
[Dmaxopt,efoopt]=maxDiffusion(sim, laserpowerLnorm,repetitions,[],"optimized")
title("D_s="+Dmaxp6+", D_o="+Dmaxopt)

f.Position=posfig;
savefig(f,[savepath filesep 'optim2.fig'])
saveas(f,[savepath filesep 'optim2.eps'],'epsc')
%no dead times
% also decrease pattern time

%%
% 
% sim.sequence=defaultsequence;
% sim.deadtimes.point=0; sim.deadtimes.estimator=0;
% sim.sequence.Itr(4).patDwellTime=2e-5; %20 us
% sim.sequence.Itr(4).Mode.pattern='hexagon';
% sim.deadtimes.point=0; sim.deadtimes.estimator=0;
% sim.makepatterns;
% [Dmaxnod,efonod]=maxDiffusion(sim, laserpower,repetitions,[],"fast hex")
% sim.deadtimes=defaultdeadtimes; %back to default;

%%
%better estimator
figure(314); clf
% standard
sim.sequence=defaultsequence;
sim.makepatterns;
[Dmaxp6b,efop6b]=maxDiffusion(sim, laserpower,repetitions,[],"standard")

sim.sequence=defaultsequence;
sim.estimators(4).function='est_qLSQiter2D';
sim.estimators(4).par{3}=false; % no ccr check
sim.estimators(4).par(1)=[];
sim.makepatterns;

[D0est,efoest]=maxDiffusion(sim, laserpower,repetitions,[],"iterLSQ")
sim.estimators=defaultestimators;



%% helper functions
function [Dmax, efoD0,rmseDmax]=maxDiffusion(sim, laserpower,repetitions,Ds,dname,fig)
if nargin<4 || isempty(Ds)
    Ds=[0 0.1 0.2 0.25:0.25:8];
end
if nargin<5 || isempty(dname)
    dname="data";
end

if nargin<6
    fig=gcf;
end
rmse=zeros(length(Ds),repetitions,3);efo=zeros(length(Ds),repetitions); efostart=efo;
indlosta=zeros(length(Ds),repetitions);

maxerr=100;

for d=1:length(Ds)
    for k=1:repetitions
        sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;

        fl=FlMoveBleach;
        fl.photonbudget=inf;
        fl.brightness=1000*laserpower;
        
        D=Ds(d); %um^2/s
        fl.makediffusion(D,0.003)
        sim.fluorophores=fl;
        out=sim.runSequence(repetitions=1,resetfluorophores=true);
        err=sqrt((out.loc.xnm-out.loc.xfl1).^2+(out.loc.ynm-out.loc.yfl1).^2);
        indlost=find(err<maxerr, 1,'last');
        if isempty(indlost)
            indlost=0;
            rmse(d,k,:)=NaN*[1 1 1];
            efo(d,k)=NaN;
        else
            filter=out.loc.itr==max(out.loc.itr);
            filter(indlost+1:end)=false;
            sr=sim.summarize_results(out,display=false,filter=filter);
            rmse(d,k,:)=sr.rmse;
            efo(d,k)=mean(out.loc.efo(filter),'omitnan');
        end
        indlosta(d,k)=indlost;
        % filter(max(2,length(filter)-10):end)=false;
        % efostart(d,k)=mean(out.loc.efo(filter),'omitnan'); %take out 10 last locs where fluorophore might be lost
    end
end

% rmsem=squeeze(mean(rmse,2,'omitnan'));
% indconverged=rmse<min(rmse(1,1,1)*3,75);

isconverged=indlosta>sim.sequence.locLimit-10; %safe margin
fconverged=squeeze(sum(isconverged,2))/repetitions;
% fconverged=mean(fconverged(:,1:2),2);

startind=find(fconverged>0.5,1,'last');
if isempty(startind)
    startind=1;
end
startD=Ds(startind);
fitfun=fittype(@(D,alpha,x) .5-erf(alpha*(x-D))/2 ); %sigmoidal
fp=fit(Ds',fconverged,fitfun,StartPoint=[startD,1]);
Dmax=fp.D;

indconv=find(Ds<=Dmax,1,'last');


% rmsec=rmse; rmsec(~isconverged)=NaN;
rmsem=squeeze(mean(rmse,2,'omitnan'));

efoDmax=mean(efo(indconv,:),'omitnan');
efoD0=mean(efo(1,:),'omitnan');
rmseDmax=mean(rmsem(indconv,1:2),'omitnan');
if any(size(rmseDmax)~=1)
    rmseDmax=NaN;
end
figure(fig)
subplot(1,3,1)

hp=plot(Ds,fconverged,DisplayName=dname);
hold on
ylabel("fraction tracked")
xlabel('Diffusion coefficient um^2/s')
title("Dmax: "+ string(Dmax) +" µm2/s")
plot(Ds,fp(Ds),'Color',hp.Color,LineWidth=1, DisplayName=dname+"fit")
legend

subplot(1,3,2)

plot(Ds,mean(efo,2,'omitnan'),'Color',hp.Color,DisplayName=dname)
hold on
xlabel('Diffusion coefficient um^2/s')
ylabel("efo kHz")
legend

subplot(1,3,3)

plot(Ds,mean(rmsem(:,1:2),2,'omitnan'),'Color',hp.Color,DisplayName=dname)
hold on
ylabel("RMSE of converged (nm)")
xlabel('Diffusion coefficient um^2/s')
legend

% Dmax, efoD0,rmseDmax
end


function [msd,dt]=getmsd(x,y,time,maxd,dt)
% maxp=max(time)/dt;
timeint=min(time):dt:max(time);
intm="nearest";
intm="linear";
xint=interp1(time,x,timeint,intm);
yint=interp1(time,y,timeint,intm);
% maxd=10;
msd=zeros(maxd,1);

for di=1:maxd
    msd(di)=mean((xint(1:end-di)-xint(1+di:end)).^2+(yint(1:end-di)-yint(1+di:end)).^2);
end
end

function [Dfit, off,msdt]=msdfit(msd,dt,linec)
msdt=(0:dt:length(msd)*dt)';
plot(msdt(2:end)/1e3,msd/1e6,linec)
fp=fit(msdt(2:end),msd(1:end),"poly1");
hold on;plot(msdt/1e3,fp(msdt)/1e6,linec+"--")
Dfit=fp.p1/1000/4; %um^2/s
off=sqrt(fp.p2);
% title("Dfit: "+string(Dfit)+" µm2/s, offset: " + string(fp.p2)+ " nm, " + "sqrt(off): " +string(off))
xlabel('time(s)')
ylabel('MSD (µm^2)')
end