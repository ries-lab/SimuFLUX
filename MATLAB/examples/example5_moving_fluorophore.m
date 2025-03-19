%% moving, bleaching fluorophore and tracking with Abberior sequence
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"SimSequencefile")
    sim=SimSequencefile;
end

fname='Tracking_2D.json';
fname2='PSFvectorial2D.json'; %use a PSF that is defined via a json file
sim.loadsequence(fname,fname2);
sim.makepatterns;

%% make diffusing, bleaching fluorophores
figure(250)
tiledlayout(1,2,"TileSpacing","tight")
nexttile

sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;sim.background=0;
fl=FlMoveBleach;
fl.photonbudget=100000;
updatetime=0.01; %ms
D=0.5; %um^2/s
fl.makediffusion(D,updatetime)
%diffusion coefficient, update time args.startpos,dim, numpoints, buondarybox
sim.fluorophores=fl;
out=sim.runSequence("repetitions",1);
sim.plotpositions(out,xvalues="time");
nexttile
plot(out.loc.xnm,out.loc.ynm,'k',out.loc.xfl1,out.loc.yfl1,'r',out.loc.xgalvo+out.loc.xeod,out.loc.ygalvo+out.loc.yeod,'g')
legend('estimated','fluorophore','galvo+EOD')
axis equal
xlabel("x (nm)")
ylabel("y (nm)")
title("diffusion")

%% make stepping fluorophore
% fig. 1
sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
fl2=FlMoveBleach;
fl2.photonbudget=5000;
fl2.brightness=200;
updatetime=0.01; %us
stepsize=16; %nm
dwelltime=28; %ms
fl2.makesteps(stepsize,dwelltime, updatetime,angle=0,startpos=[50,0,0])
      % args.startpos, dim, numpoints,angle (degree);

sim.fluorophores=fl2;sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
out=sim.runSequence("repetitions",1);

nexttile
% figure(99)
sim.plotpositions(out,xvalues="time");
title("stepping")

%% instabilities: vibrations
fl3=FlMoving(brightness=50000); %collect more photons
fl3.posmode='function';

figure(261)
tiledlayout(1,4); nexttile
frequencies=[0.05 0.1, 0.2, 0.5, 1, 2, 5, 10]; %kHz

% frequencies=1;
amplitude=5; %nm

stdx=0*frequencies;stdxrel=stdx;
for k=1:length(frequencies)
    posfl=[0 0 0];
    fl3.posfunction={@(t) amplitude*sin(frequencies(k)*t)+posfl(1), @(t) 0*t+posfl(2), @(t) 0*t+posfl(3)};    
    sim.fluorophores=fl3;sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
    out=sim.runSequence("repetitions",1);
    filter=out.loc.itr==max(out.loc.itr)&out.loc.vld==1;
    sr=sim.summarize_results(out,display=false,filter=filter);
    stdx(k)=sr.stdraw(1);
    stdxrel(k)=sr.stdraw(1)/sr.sCRB(1);
end
semilogx(frequencies,stdxrel)
xlabel('frequncy (kHz)')
ylabel('std(x)/sCRB(x)')
title("vibrations: frequency")

frequency=0.1;
amplitudes=0:2:10;
stdxa=0*amplitudes;stdxrela=stdxa;
for k=1:length(amplitudes)
    posfl=[0 0 0];
    fl3.posfunction={@(t) amplitudes(k)*sin(frequency*t)+posfl(1), @(t) 0*t+posfl(2), @(t) 0*t+posfl(3)};    
    sim.fluorophores=fl3;sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
    out=sim.runSequence("repetitions",1);
    filter=out.loc.itr==max(out.loc.itr)&out.loc.vld==1;
    sr=sim.summarize_results(out,display=false,filter=filter);
    stdxa(k)=sr.stdraw(1);
    stdxrela(k)=sr.stdraw(1)/sr.sCRB(1);
end
nexttile
plot(amplitudes,stdxrela)
xlabel('amplitude (nm)')
ylabel('std(x)/sCRB(x)')
title("vibrations: amplitude")