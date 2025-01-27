% Example getting started
% Define programmatically a simple MINFLUX experiment and look at the
% results

fl=Fl_static; %define a static fluorophore
fl.pos=[10 0 0];
fl.brightness=1000; %kHz if excited at the center of a Gaussian beam

psf_donut=PSF_donut2D; %here you define a PSF. In this case, an analytical 2D donut PSF

sim=Sim_Simulator(fl); %make a simulator and attach fluorophore

numberOfLocalizations=1000;

%define scan pattern
L=75; %size of scan pattern
orbitpoints=6; %number of probing points in orbit
probecenter=true; %should we also probe the center?
laserpower=5; %relative, increases brightness
pointdwelltime=100; %measurement time in each point
repetitions=2; %how often to repeat the pattern scan
sim.definePattern("donut", psf_donut, makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=dwelltime,laserpower=laserpower,repetitions=repetitions)
% sim.definePattern(key, PSF_object, arguments...)

%we need an estimator. Define as component
fwhm=360;% size of the donut, needed for proper estimation. 
sim.defineComponent("estgauss","estimator",@estimators,parameters={"donut2D",sim.patterns("donut").pos,L,fwhm},dim=1:2);
% sim.defineComponent(key,type (estimator),function handle of estimator function,parameters);

%sequence: 
seq={"donut","estgauss"};
out=sim.runSequence(seq);
% out.loc: localizations
% out.fluorophores: position of fluorophores
% out.raw: photon measurements

sim.displayresults("donut", out) %display summary of simulation
% sim.displayresults(patternkey, out)