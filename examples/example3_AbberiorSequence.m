% Run an abberior sequence file and see how EOD an Galvo move
if ~exist('sim','var') || ~isa(sim,"Sim_sequencefile")
    sim=Sim_sequencefile;
end
fl=Fl_static; %define a static fluorophore
fl.pos=[200 50 0];
fl.brightness=10000; %kHz if excited at the center of a Gaussian beam
sim.fluorophores=fl;
sim.posgalvo=[0 0 0];

fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit100__locLim50_.json';
fname2='PSFvectorialJonas.json'; %use a PSF that is defined via a json file
sim.loadsequence(fname,fname2);
sim.makepatterns;
out=sim.runSequence;

figure(88);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl1);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
plot(out.loc.loccounter,out.loc.xeod)
xlabel('time (itr)')
ylabel('x position(nm)')
legend('estimated', 'fluorophore','xgalvo','EOD')

sim.displayresults(out) %display summary of simulation
% sim.displayresults(patternkey, out)