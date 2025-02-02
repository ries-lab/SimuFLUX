% Investigate the effect of background
% Causes: imperfect zero, autofluorescence background, nearby fluorophore
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);
sim=Simulator;

numberOfLocalizations=1000;
L=75;
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=4, ...
    probecenter=true,orbitL=L,laserpower=100)
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
seq={"donut","estdonut"};
%% no background
fl1=FlStatic;
fl1.pos=[20 0 0];
psf_vec.zerooffset=0;
sim.fluorophores=fl1;
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp("no background")
sim.summarize_results(out);

%% imperfect zero
psf_vec.zerooffset=0.2;
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp("zerooffset: ")
sim.summarize_results(out);

%% Autofluorescence background
psf_vec.zerooffset=0.0;
sim.background=50;
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp("fluorescence background: ")
sim.summarize_results(out);

%% nearby fluorophore
fl2=FlStatic;
fl2.pos=[50 100 600];

fc=FlCollection;
fc.add([fl1, fl2])
psf_vec.zerooffset=0;
sim.fluorophores=fc; %make a simulator and attach fluorophore

out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp("nearby fluorophore: ")
sim.summarize_results(out);

