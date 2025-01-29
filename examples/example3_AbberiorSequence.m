% Run an abberior sequence file and see how EOD an Galvo move
if ~exist('sim','var') || ~isa(sim,"Sim_sequencefile")
    sim=Sim_sequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];
end
laserpower=50;
fl=Fl_bleach; %define a bleaching fluorophore
fl.photonbudget=10000;
fl.pos=[200 50 0];
fl.brightness=100*laserpower; %kHz if excited at the center of a Gaussian beam
sim.fluorophores=fl;


% fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit100__locLim50_.json';
fname='Tracking_2D.json';
fname2='PSFvectorial2D.json'; %use a PSF that is defined via a json file
sim.loadsequence(fname,fname2);
sim.makepatterns;
out=sim.runSequence("repetitions",1);

figure(88);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl1);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
plot(out.loc.loccounter,out.loc.xeod)
xlabel('time (itr)')
ylabel('x position(nm)')
legend('estimated', 'fluorophore','xgalvo','EOD')

sim.displayresults(out) %display summary of simulation
% sim.displayresults(patternkey, out)