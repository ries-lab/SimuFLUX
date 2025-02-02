% Example getting started
% Define programmatically a simple MINFLUX experiment and look at the
% results
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path

fl=FlStatic; %define a static fluorophore
fl.pos=[10 0 0];
fl.brightness=1000; %kHz if excited at the center of a Gaussian beam

psf_donut=PsfDonut2D; %here you define a PSF. In this case, an analytical 2D donut PSF

sim=Simulator(fl); %make a simulator and attach fluorophore

numberOfLocalizations=1000;

%define scan pattern
L=75; %size of scan pattern
orbitpoints=6; %number of probing points in orbit
probecenter=true; %should we also probe the center?
laserpower=5; %relative, increases brightness
pointdwelltime=0.1; %ms, measurement time in each point
repetitions=2; %how often to repeat the pattern scan
sim.definePattern("donut", psf_donut, makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
% sim.definePattern(key, PSF_object, arguments...)

%we need an estimator. Define as component
fwhm=360;% size of the donut, needed for proper estimation. 
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,fwhm},dim=1:2);
% sim.defineComponent(key,type (estimator),function handle of estimator function,parameters);

%sequence: 
seq={"donut","estdonut"};
out=sim.runSequence(seq,maxlocs=numberOfLocalizations);
%args: maxlocs: number of times the pattern is scanned (one trace)
%       repetitions: number of times the sequence is repeated
% out.loc: localizations
% out.fluorophores: position of fluorophores
% out.raw: photon measurements

sim. summarize_results(out); %display summary of simulation
