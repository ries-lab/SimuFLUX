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

%% make diffusing, bleaching
sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
fl=FlMoveBleach;
fl.photonbudget=100000;
updatetime=0.01; %ms
D=0.5; %um^2/s
fl.makediffusion(D,updatetime)
%diffusion coefficient, update time args.startpos,dim, numpoints, buondarybox
sim.fluorophores=fl;


out=sim.runSequence("repetitions",1);

sim.plotpositions(out,figure=211,xvalues="time");

%% make stepping fluorophore
sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
fl2=FlMoveBleach;
fl2.photonbudget=10000;
updatetime=0.01; %us
stepsize=16; %nm
dwelltime=5; %ms
fl2.makesteps(stepsize,dwelltime, updatetime,angle=45)
      % args.startpos, dim, numpoints,angle (degree);

sim.fluorophores=fl2;sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
out=sim.runSequence("repetitions",1);
sim.plotpositions(out,figure=212);

