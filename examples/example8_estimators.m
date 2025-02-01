% Compare different estimators with and without backgroud
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);
fl=FlStatic;
fl.pos=[20 0 0];
sim=Simulator(fl);

numberOfLocalizations=1000;
L=75;
probecenter=true;
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=4, ...
    probecenter=probecenter,orbitL=L,laserpower=100)


%% no background
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
seq={"donut","estdonut"};
psf_vec.zerooffset=0;
xcoords=0:5:100;
figure(293); statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true);


% out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
% disp("est_donut2d")
% sim.displayresults(out);

