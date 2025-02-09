% Compare different estimators with and without backgroud
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);

fl=FlStatic;
fl.brightness=1000000;
fl.pos=[20 0 0];
sim=Simulator(fl);

numberOfLocalizations=100;
L=75;
xcoords=0:2:2*L;
probecenter=true;
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=4, ...
    probecenter=probecenter,orbitL=L,laserpower=100)
% sim.calculateCRB("donut")
% sim.calculateCRBscan("donut")

ax1v="pos"; ax2v="rmse";
%% no background, simple estimator
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
seq={"donut","estdonut"};
psf_vec.zerooffset=0;

figure(293); statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax2=ax2v, ax1=ax1v,clearfigure=true);

%% iterative
sim.defineComponent("estiter","estimator",@est_quad2Diter4points,parameters={L,15},dim=1:2);
seq={"donut","estiter"};
% sim.fluorophores.pos=[30 20 0];
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax2=ax2v, ax1=ax1v,clearfigure=false);
% out=sim.runSequence(seq,maxlocs=1);
% sim.summarize_results(out);
%% no background, debiased estimator
sim.defineComponent("estmLMS","estimator",@est_donut2d_mLMS,parameters={sim.patterns("donut").pos,L,360,[1.3 3.6]},dim=1:2);
seq={"donut","estmLMS"};
hold off
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax2=ax2v, ax1=ax1v,clearfigure=false);

%% no background, direct estimator
sim.defineComponent("estdirect","estimator",@est_quadraticdirect1D,parameters={L},dim=1);
seq={"donut","estdirect"};
hold off
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax2=ax2v, ax1=ax1v,clearfigure=false);

%% no background, quadratic approx
sim.defineComponent("estq","estimator",@est_quadraticND,parameters={sim.patterns("donut").pos},dim=1:2);
seq={"donut","estq"};
hold off
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax2=ax2v, ax1=ax1v,clearfigure=false);

%% no background, simple estimator
sim.background=0;
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360,"bg_phot"},dim=1:2);
seq={"donut","estdonut"};
psf_vec.zerooffset=0;
xcoords=0:5:50;
figure(294); statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax2=["biasrel"]);

%background, simple estimator
sim.background=30;
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360,0},dim=1:2);
seq={"donut","estdonut"};
psf_vec.zerooffset=0;
% xcoords=0:2:100;
figure(294); statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax2=["biasrel"],clearfigure=false);

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

