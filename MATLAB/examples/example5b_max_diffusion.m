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
% L=150;
sim.sequence.locLimit=100;
% sim.sequence.Itr(4)=[];
% sim.sequence.damping=1;

%iterative estimator:
% sim.estimators(4).function='est_qLSQiter2D';
% sim.estimators(4).par{3}=false; % no ccr check

% sim.sequence.Itr(4).phtLimit=20;
% sim.sequence.Itr(4).patDwellTime=0.1e-3;
% sim.sequence.Itr(4).patGeoFactor=L/360;
 % L=itrs(itr).patGeoFactor*360; %nm;
% sim.deadtimes.point=0; sim.deadtimes.estelimator=0;
sim.makepatterns;
laserpower=0.6;
updatetime=0.01; %ms

maxerr=100;

%% plot example track
figure(256); clf
D=3;
sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
fl=FlMoveBleach;
fl.photonbudget=inf;
fl.brightness=1000*laserpower;
fl.makediffusion(D,updatetime)
sim.fluorophores=fl;
out=sim.runSequence(repetitions=1,resetfluorophores=true);
err=sqrt((out.loc.xnm-out.loc.xfl1).^2+(out.loc.ynm-out.loc.yfl1).^2);
indlost=find(err<maxerr, 1,'last');
subplot(2,3,1)
hold off; plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl1)%;hold on; plot(out.loc.loccounter,out.loc.xgalvo)
xlabel('time (localizations)'); ylabel('x position (nm)'); legend('estimated position','fluorophore position')
plot([indlost,indlost],[min(out.loc.xnm),max(out.loc.xnm)])
subplot(2,3,2)
hold off; plot(out.loc.loccounter, out.loc.ynm);hold on; plot(out.loc.loccounter,out.loc.yfl1)%;hold on; plot(out.loc.loccounter,out.loc.ygalvo)
plot([indlost,indlost],[min(out.loc.ynm),max(out.loc.ynm)])
xlabel('time (localizations)'); ylabel('x position (nm)'); 
subplot(2,3,3)
hold off; plot(out.loc.xnm(1:end-1), out.loc.ynm(1:end-1));hold on; plot(out.loc.xfl1,out.loc.yfl1)%;hold on; plot(out.loc.loccounter,out.loc.ygalvo)
xlabel('x position (nm)'); ylabel('y position (nm)'); axis equal


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

subplot(2,3,4)
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
msdwin=10;
msdmf=zeros(msdwin,1);msdfl=msdmf;
imsd=0;
repetitionsmsd=100;
dtmsd=[];
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
subplot(2,3,5);
hold off
[Dfl,offfl]=msdfit(msdfl,dtmsd,'r');
[Dmf,offmf]=msdfit(msdmf,dtmsd,'b');
legend("fluorophore","D="+string(Dfl)+"µm2/s","Minflux","D="+string(Dmf)+"µm2/s")

savefig(gcf,[savepath filesep 'tracks_abort_msd.fig'])

%% investigate maximum diffusion coefficient for various conditions
% change default settings
fname='Tracking_2D.json';
fname2='PSFvectorial2D.json'; %use a PSF that is defined via a json file
sim.loadsequence(fname,fname2);
sim.sequence.locLimit=100;

repetitions=20;
Dtest=0:0.25:9;

sim.makepatterns;

figure(258); clf
laserpowers=[0.05 0.1 0.2:0.2:0.8 1:0.5:3];
clear Dmaxl efoDmaxl rmseDmaxl
figure(301);clf
for k=1:length(laserpowers)
    [Dmaxl(k),efoDmaxl(k),rmseDmaxl(k)]=maxDiffusion(sim, laserpowers(k),repetitions, Dtest)
    drawnow
end
title('laserpowers')
savefig(gcf,[savepath filesep 'laserpowers.fig'])

figure(258)
subplot(3,4,1)
plot(laserpowers,Dmaxl)
xlabel('laserpower (a.u.)')
ylabel('max Diffusion D µm^2/s')
title('laserpower')
text(0,0.35,num2str(efoDmaxl','%2.0f'))

subplot(3,4,5)
plot(laserpowers,rmseDmaxl)
xlabel('laserpower (a.u.)')
ylabel('rmse at Dmax (nm)')
title('laserpower')

subplot(3,4,9)
plot(efoDmaxl,Dmaxl)
xlabel('efo (kHz)')
ylabel('max Diffusion D µm^2/s')
title('laserpower')
text(0,0.35,num2str(efoDmaxl','%2.0f'))

subplot(3,4,10)
plot(efoDmaxl,rmseDmaxl)
xlabel('efo (kHz)')
ylabel('rmse at Dmax (nm)')
title('laserpower')
% text(0,0.35,num2str(efoDmaxl','%2.0f'))

%redo with no dead times and faster pattern dwell times
deadtimesold=sim.deadtimes;
patterntimeold=sim.sequence.Itr(4).patDwellTime;
sim.deadtimes.point=0; sim.deadtimes.estelimator=0;
sim.sequence.Itr(4).patDwellTime=2e-5; %20 us
sim.makepatterns();
clear Dmaxlf efoDmaxlf rmseDmaxlf
figure(302); clf
for k=1:length(laserpowers)
    [Dmaxlf(k),efoDmaxlf(k),rmseDmaxlf(k)]=maxDiffusion(sim, laserpowers(k),repetitions, Dtest)
end
sim.deadtimes=deadtimesold; %revert to standard
sim.sequence.Itr(4).patDwellTime=patterntimeold;
sim.makepatterns();
title('laserpower fast')
legend(string(laserpowers'))
savefig(gcf,[savepath filesep 'laserpowersfast.fig'])

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
plot(efoDmaxlf,Dmaxlf)


subplot(3,4,10)
hold on
plot(efoDmaxlf,rmseDmaxlf)


Dtest=0:0.25:6;
% L scan pattern size
clear DmaxL efoDmaxL rmseDmaxL
Ls=[40 75 150 200 250 300];
patGeoFactorOld=sim.sequence.Itr(4).patGeoFactor;
figure(303); clf
for k=1:length(Ls)
    sim.sequence.Itr(4).patGeoFactor=Ls(k)/360;
    [DmaxL(k),efoDmaxL(k),rmseDmaxL(k)]=maxDiffusion(sim, laserpower,repetitions, Dtest)
end
title('L')
legend(string(Ls'))
savefig(gcf,[savepath filesep 'Ls.fig'])

sim.sequence.Itr(4).patGeoFactor=patGeoFactorOld;

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


clear Dmaxpl efoDmaxpl rmseDmaxl
photlim=[5 10 20 30 50 100];
figure(304);clf
for k=1:length(photlim)
    sim.sequence.Itr(4).phtLimit=photlim(k);
    [Dmaxpl(k),efoDmaxpl(k),rmseDmaxpl(k)]=maxDiffusion(sim, laserpower,repetitions, Dtest)
end
title('photon limit')
legend(string(photlim'))
savefig(gcf,[savepath filesep 'photlim.fig'])

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



% dwell time
dwelltimes=[10 20 50 100 200 500 1000]; %us
% dwelltimes=10;

patterntimeold=sim.sequence.Itr(4).patDwellTime;
clear Dmaxdt efoDmaxdt rmseDmaxdt
figure(305);clf
for k=1:length(dwelltimes)
    sim.sequence.Itr(4).patDwellTime=dwelltimes(k)*1e-6;
    sim.makepatterns();
    [Dmaxdt(k),efoDmaxdt(k),rmseDmaxdt(k)]=maxDiffusion(sim, laserpower,repetitions, Dtest)
end
title('pattern dwell time')
 legend(string(dwelltimes'))
 savefig(gcf,[savepath filesep 'dwelltimes.fig'])
%redo with fast instrument
sim.deadtimes.point=0; sim.deadtimes.estimator=0;

figure(306);clf
clear Dmaxdtf efoDmaxdtf rmseDmaxdtf
for k=1:length(dwelltimes)
    sim.sequence.Itr(4).patDwellTime=dwelltimes(k)*1e-6;
    sim.makepatterns();
    [Dmaxdtf(k),efoDmaxdtf(k),rmseDmaxdtf(k)]=maxDiffusion(sim, laserpower,repetitions, Dtest)
end
title('pattern dwell time fast')
 legend(string(dwelltimes'))
 savefig(gcf,[savepath filesep 'dwelltimesfast.fig'])

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
sim.deadtimes=deadtimesold;
sim.sequence.Itr(4).patDwellTime=patterntimeold;
sim.makepatterns();

savefig(gcf,[savepath filesep 'D_scan_vars.fig'])
%% compare different conditions
repititions=10;
figure(259)
clf
% reference
[Dmaxp,efop]=maxDiffusion(sim, laserpower,repetitions)

% 3 instead of 6 probing positions
sim.sequence.Itr(4).Mode.pattern='triangle';
sim.makepatterns;
[Dmaxp,efop]=maxDiffusion(sim, laserpower,repetitions)

% 4 instead of 6 probing positions
sim.sequence.Itr(4).Mode.pattern='square';
sim.makepatterns;
[Dmaxp,efop]=maxDiffusion(sim, laserpower,repetitions)

%no dead times
% also decrease pattern time
patterntimeold=sim.sequence.Itr(4).patDwellTime;
sim.deadtimes.point=0; sim.deadtimes.estelimator=0;
sim.sequence.Itr(4).patDwellTime=2e-5; %20 us
sim.sequence.Itr(4).Mode.pattern='hexagon';
deadtimesold=sim.deadtimes;
sim.deadtimes.point=0; sim.deadtimes.estelimator=0;
sim.makepatterns;
[Dmaxnod,efonod]=maxDiffusion(sim, laserpower,repetitions)
sim.deadtimes=deadtimesold;
sim.sequence.Itr(4).patDwellTime=patterntimeold;
sim.makepatterns;
%
%better estimator
simestold=sim.estimators(4);
sim.estimators(4).function='est_qLSQiter2D';
sim.estimators(4).par{3}=false; % no ccr check
sim.estimators(4).par(1)=[];

[Dmaxest,efoest]=maxDiffusion(sim, laserpower,repetitions)
sim.estimators(4)=simestold;

legend('reference','triangle','fast','square','iterative LSQ')
savefig(gcf,[savepath filesep 'compare_conditions.fig'])

%% helper functions
function [Dmax, efoDmax,rmseDmax]=maxDiffusion(sim, laserpower,repetitions,Ds)
if nargin<4
    Ds=0:0.25:8;
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
        fl.makediffusion(D,0.01)
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
fitfun=fittype(@(D,alpha,x) .5-erf(alpha*(x-D))/2 ); %sigmoidal
fp=fit(Ds',fconverged,fitfun,StartPoint=[2,1]);
Dmax=fp.D;

indconv=find(Ds<=Dmax,1,'last');


% rmsec=rmse; rmsec(~isconverged)=NaN;
rmsem=squeeze(mean(rmse,2,'omitnan'));

efoDmax=mean(efo(indconv,:),'omitnan');
rmseDmax=mean(rmsem(indconv,1:2),'omitnan');
subplot(1,3,1)

hp=plot(Ds,fconverged);
hold on
ylabel("fraction tracked")
xlabel('Diffusion coefficient um^2/s')
title("Dmax: "+ string(Dmax) +" µm2/s")
plot(Ds,fp(Ds),'Color',hp.Color)
subplot(1,3,2)

plot(Ds,mean(efo,2,'omitnan'))
hold on
xlabel('Diffusion coefficient um^2/s')
ylabel("efo kHz")

subplot(1,3,3)

plot(Ds,mean(rmsem(:,1:2),2,'omitnan'))
hold on
ylabel("RMSE of converged (nm)")
xlabel('Diffusion coefficient um^2/s')
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

function [Dfit, off]=msdfit(msd,dt,linec)
msdt=(0:dt:length(msd)*dt)';
plot(msdt(2:end),msd,linec)
fp=fit(msdt(2:end),msd(1:end),"poly1");
hold on;plot(msdt,fp(msdt),linec+"--")
Dfit=fp.p1/1000/4; %um^2/s
off=sqrt(fp.p2);
% title("Dfit: "+string(Dfit)+" µm2/s, offset: " + string(fp.p2)+ " nm, " + "sqrt(off): " +string(off))
end