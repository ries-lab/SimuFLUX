%% moving, bleaching fluorophore and tracking with Abberior sequence
%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"Sim_sequencefile")
    sim=Sim_sequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
end
fname='Tracking_2D.json';
fname2='PSFvectorial2D.json'; %use a PSF that is defined via a json file
sim.loadsequence(fname,fname2);
sim.makepatterns;

%% make diffusing, bleaching
fl=Fl_move_bleach;
fl.photonbudget=10000;
updatetime=10; %us
D=.1; %um^2/s
fl.makediffusion(D,updatetime)
%diffusion coefficient, update time args.startpos,dim, numpoints, buondarybox
sim.fluorophores=fl;


out=sim.runSequence("repetitions",1);

sim.plotpositions(out)

%% make stepping fluorophore
fl2=Fl_move_bleach;
fl2.photonbudget=10000;
updatetime=10; %us
stepsize=16; %nm
dwelltime=5000; %us
fl2.makesteps(stepsize,dwelltime, updatetime,angle=45)
      % args.startpos, dim, numpoints,angle (degree);

sim.fluorophores=fl2;sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
out=sim.runSequence("repetitions",1);
sim.plotpositions(out)

