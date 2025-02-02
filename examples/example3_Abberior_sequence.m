% Run an abberior sequence file and see how EOD an Galvo move
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist('sim','var') || ~isa(sim,"SimSequencefile")
    sim=SimSequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];
end
laserpower=50;
fl=FlBleach; %define a bleaching fluorophore
fl.photonbudget=10000;
fl.pos=[200 50 0];
fl.brightness=200*laserpower; %kHz 
sim.fluorophores=fl;


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
indf=out.loc.itr==max(out.loc.itr) & out.loc.vld==1;
sim.summarize_results(out,filter=indf); %display summary of simulation
