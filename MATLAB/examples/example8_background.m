% Investigate the effect of background
% Causes: imperfect zero, autofluorescence background, nearby fluorophore
%XXX background handling of estimator wrong, this makes this not well
%interpretable. Use bg calibration and background in estimator
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);
estimatebackground=true;
sim=Simulator;
laserpower=10;
numberOfLocalizations=1000;
L=75;
probecenter=true;
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=4, ...
    probecenter=probecenter,orbitL=L,laserpower=laserpower,pointdwelltime=0.1)
sim.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
sim.defineComponent("bg","background",@backgroundsubtractor,parameters={"background_estimated"});
seq={"donut","bg","estdonut"};
% bgf=sim.patterns('donut').backgroundfac(1); %background used for simulation
fltestpos=[10 0 0];
fl1=FlStatic;
fl1.pos=fltestpos;
sim.fluorophores=fl1;
figure(281);
tiledlayout(1,3); 
xs=0:5:65;
%% no background
sim.background_estimated=0;
fl1.pos=fltestpos;
psf_vec.zerooffset=0;
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp("no background")
sim.summarize_results(out);

nexttile(1)
hold off
sim.scan_fov(seq,xs, tag="no bg",ax1=["std"],linestyle='k');
hold on
nexttile(3)
hold off
sim.scan_fov(seq,xs, tag="no bg",ax1=["bias"],linestyle='k');
hold on
nexttile(2)
hold off
sim.scan_fov(seq,xs, tag="no bg",ax1=["rmse"],linestyle='k');
hold on
%% imperfect zero
psf_vec.zerooffset=0.01;
% estimate background by scanning centered fluorophore
fl1.pos=[0 0 0];
sim.background_estimated=0;
outtest=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
bgphot=mean(outtest.raw(:,end));
bgzeroscan=bgphot/outtest.par{1}.pointdwelltime(1);
bgestdirect=fl1.brightness*psf_vec.zerooffset*laserpower*psf_vec.PSFs("pinhole633    0    0").interp([0 0 0]); %illustration how zero offset is converted into background
sim.background_estimated=bgzeroscan*estimatebackground;
fl1.pos=fltestpos;
% sim.background_estimated=psf_vec.zerooffset*fl.brightness*bgf; %in general, the GT background is not known but needs to be calibrated 
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp("zerooffset: ")
sim.summarize_results(out);

nexttile(1)
sim.scan_fov(seq,xs, tag="imperfect zero",ax1=["std"],linestyle='r--');
nexttile(3)
sim.scan_fov(seq,xs, tag="imperfect zero",ax1=["bias"],linestyle='r--');
nexttile(2)
sim.scan_fov(seq,xs, tag="imperfect zero",ax1=["rmse"],linestyle='r--');

%% Autofluorescence background
psf_vec.zerooffset=0.0;
sim.background=30;
sim.background_estimated=sim.background*laserpower*estimatebackground; %in general, the GT background is not known but needs to be calibrated 
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp("fluorescence background: ")
sim.summarize_results(out);


nexttile(1)
sim.scan_fov(seq,xs, tag="auto fluorescence",ax1=["std"],linestyle='b:');
xlim([0 max(xs)])
nexttile(3)
sim.scan_fov(seq,xs, tag="auto fluorescence",ax1=["bias"],linestyle='b:');
xlim([0 max(xs)])
nexttile(2)
sim.scan_fov(seq,xs, tag="auto fluorescence",ax1=["rmse"],linestyle='b:');
xlim([0 max(xs)])

%% nearby fluorophore
fl2=FlStatic;
fl2.pos=[50 50 400];

fc=FlCollection;
fc.add({fl1, fl2})
psf_vec.zerooffset=0;
sim.fluorophores=fc; %make a simulator and attach fluorophore
sim.background_estimated=0; %we don't know bg
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp("nearby fluorophore: ")
sim.summarize_results(out);

