% Compare different estimators with and without backgroud
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);

fl=FlStatic;
fl.brightness=10000; %very bright to look at bias
fl.pos=[20 0 0];
sim=Simulator(fluorophores=fl);

numberOfLocalizations=100;
L=75;
orbitpoints=4;
laserpower=100;
xcoords=0:2:L;
probecenter=true;
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,laserpower=laserpower)
% sim.calculateCRB("donut")
% sim.calculateCRBscan("donut")

ax1v="pos"; 
%% no background, simple estimator
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
seq={"donut","estdonut"};
psf_vec.zerooffset=0;

figure(270); 
tiledlayout("TileSpacing","tight"); nexttile; hold off 
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=true,tag="simple est");

%iterative
sim.defineComponent("estiter","estimator",@est_quad2Diter,parameters={L,probecenter,30},dim=1:2);
seq={"donut","estiter"};
% sim.fluorophores.pos=[30 20 0];
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true, ax1=ax1v,clearfigure=false,tag="iterative est");
% out=sim.runSequence(seq,maxlocs=1);
% sim.summarize_results(out);

% no background, direct estimator
sim.defineComponent("estdirect","estimator",@est_quadraticdirect1D,parameters={L},dim=1);
seq={"donut","estdirect"};
hold off
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false, tag="direct est",title="comparison of estimators");

if ax1v=="pos"
    plot([0 L],[0 L],'k--')
end

%% explore impact of background on estimator
%no background
fl.brightness=1000;
fl.pos=[30 30 0];
sim.background=0;
sim.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
seq={"donut","estdonut"};
psf_vec.zerooffset=0;
% xcoords=0:5:50;
nexttile; hold off
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,tag="no bg");
% out=sim.runSequence(seq,maxlocs=1);
% sim.summarize_results(out);

%background, 
sim.background=50;
% xcoords=0:2:100;
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="bg");

%background subtracted,
% bgf=sim.patterns('donut').backgroundfac(1); %background used for simulation
sim.background_estimated=sim.background*laserpower; %in general, the GT background is not known but needs to be calibrated 

%sim.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
sim.defineComponent("bg","background",@backgroundsubtractor,parameters={"background_estimated"});
seq={"donut","bg","estdonut"};
psf_vec.zerooffset=0;
% xcoords=0:2:100;
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="bg est",title="estimator bias from background");


sim.defineComponent("estimatorbg","estimator",@est_quadraticdirectbg1D,parameters={L},dim=1);
seq={"donut","estimatorbg"};
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="bg fit par",title="estimator bias from background");



if ax1v=="pos"
    plot([0 L],[0 L],'k--')
end

% out=sim.runSequence(seq,maxlocs=1);
% sim.summarize_results(out);

% sim.background=30;
% sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360,"bg_phot"},dim=1:2);
% seq={"donut","estdonut"};
% psf_vec.zerooffset=0;
% % xcoords=0:2:100;
% figure(294); statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax2=["biasrel"],clearfigure=false);
% 



% out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
% disp("est_donut2d")
% sim.summarize_results(out);

