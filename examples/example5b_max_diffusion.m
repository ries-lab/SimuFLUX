%% moving, bleaching fluorophore and tracking with Abberior sequence
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"SimSequencefile")
    sim=SimSequencefile;
end
fname='Tracking_2D.json';
fname2='PSFvectorial2D.json'; %use a PSF that is defined via a json file
sim.loadsequence(fname,fname2);
sim.sequence.locLimit=200;
% sim.sequence.Itr(4)=[];
% sim.sequence.damping=1;
sim.sequence.Itr(4).phtLimit=20;
% sim.sequence.Itr(4).patDwellTime=0.1e-3;
sim.sequence.Itr(4).patGeoFactor=0.5;
% sim.deadtimes.point=0; sim.deadtimes.estelimator=0;
sim.makepatterns;
laserpower=.4;
%% make diffusing, bleaching fluorophores
Ds=0.5:0.1:1.5;
repetitions=5;
rmse=zeros(length(Ds),repetitions,3);efo=zeros(length(Ds),repetitions);

for d=1:length(Ds)
    for k=1:repetitions
        sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;

        fl=FlMoveBleach;
        fl.photonbudget=inf;
        fl.brightness=1000*laserpower;
        updatetime=0.01; %ms
        D=Ds(d); %um^2/s
        fl.makediffusion(D,updatetime)
        sim.fluorophores=fl;
        out=sim.runSequence(repetitions=1,resetfluorophores=true);
        filter=out.loc.itr==max(out.loc.itr);
        sr=sim.summarize_results(out,display=false,filter=filter);
        rmse(d,k,:)=sr.rmse;
        efo(d,k)=mean(out.loc.efo(filter));
    end
end

rmsem=squeeze(mean(rmse,2));
indconverged=rmse<rmse(1,1,1)*3;
fconverged=squeeze(sum(indconverged,2))/repetitions;
rmsec=rmse; rmsec(~indconverged)=NaN;
rmsecm=mean(rmsec,2,'omitnan');

figure(255)
subplot(1,3,1)
plot(Ds,mean(fconverged(:,1:2),2))
hold on
ylabel("fraction tracked")
xlabel('Diffusion coefficient um^2/s')

subplot(1,3,2)
plot(Ds,mean(efo,2))
hold on
xlabel('Diffusion coefficient um^2/s')
ylabel("efo kHz")

subplot(1,3,3)
plot(Ds,mean(rmsecm(:,1:2),2))
hold on
ylabel("RMSE of converged (nm)")
xlabel('Diffusion coefficient um^2/s')
