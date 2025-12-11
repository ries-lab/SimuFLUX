% Example to get started
% Define programmatically a simple MINFLUX experiment and look at the
% results
addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath')))))); %add all folders to serach path

fl=FlStatic; %define a static fluorophore
fl.pos=[100 20 0];
fl.brightness=100000; %kHz if excited at the center of a Gaussian beam

if ~exist("psf_vec","var") 
    psf_vec=PsfVectorial; 
end
psf_vec.setpinhole("AU",1);

sim=Simulator(fluorophores=fl); %make a simulator and attach fluorophore

numberOfLocalizations=1000;

%define scan pattern
L=288; %size of scan pattern
sigma=160;
sigma=110; %tophat
orbitpoints=6; %number of probing points in orbit. 
probecenter=false; %should we also probe the center?
pinholeorbit=true;
laserpower=1; %relative, increases brightness
pointdwelltime=1/(orbitpoints+probecenter); %ms, measurement time in each point
repetitions=2; %how often to repeat the pattern scan
sim.definePattern("pho_do", psf_vec, phasemask="tophat", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,...
    laserpower=laserpower,repetitions=repetitions,pinholeorbit=pinholeorbit)
% sim.definePattern(key, PSF_object, arguments...)

%we need an estimator. Define as component
sim.defineComponent("estdonut","estimator",@est_qLSQiter2D,parameters={L,probecenter},dim=1:2);
sim.defineComponent("estgauss","estimator",@est_pinholeorbit,parameters={"patternpos", L,sigma,probecenter},dim=1:2);
% sim.defineComponent("estdonut","estimator",@est_donutLSQ1_2D,parameters={"patternpos",L,360},dim=1:2);

%sequence: 
seq={"pho_do","estgauss"};

out=sim.runSequence(seq,maxlocs=numberOfLocalizations);

sim.summarize_results(out); %display summary of simulation
xs=0:20:200;
figure(22); clf
sim.scan_fov(seq,xs, tag="tophat",ax1="pos",linestyle='r',maxlocs=numberOfLocalizations);
plot(xs,xs,xs,0*xs)
