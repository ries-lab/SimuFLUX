% Compare different estimators with and without backgroud
% Figure 2
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);

fl=FlStatic;
fl.brightness=10000; %very bright to look at bias
fl.pos=[0 0 0];
sim=Simulator(fluorophores=fl);

numberOfLocalizations=1000;
L=75;
orbitpoints=4;
laserpower=100;
xcoords=0:5:L;
probecenter=true;
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,laserpower=laserpower)

ax1v="pos"; 
% ax1v="rmse"; 
%% no background, simple estimator
sim.background=0;sim.background_estimated=0;
sim.defineComponent("estdonut","estimator",@est_donutLSQ1_2D,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
seq={"donut","estdonut"};
psf_vec.zerooffset=0;

figure(270); 
tiledlayout(1,2,"TileSpacing","tight"); nexttile(1); hold off 
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=true,tag="simple est");
phottxt="photons. Simple: "+ string(statout.phot(1));
%iterative
sim.defineComponent("estiter","estimator",@est_qLSQiter2D,parameters={L,probecenter,20},dim=1:2);
seq={"donut","estiter"};
% sim.fluorophores.pos=[30 20 0];
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true, ax1=ax1v,clearfigure=false,tag="iterative est");
phottxt=phottxt+", iter: "+ string(statout.phot(1));


sim.defineComponent("estimator2p","estimator",@est_qDirectFitBg1D,parameters={L,probecenter},dim=1);
seq={"donut","estimator2p"};
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="direct eq");
phottxt=phottxt+", direct: "+ string(statout.phot(1));

sim.defineComponent("estimatorMLE","estimator",@est_qMLE1D,parameters={L},dim=1);
seq={"donut","estimatorMLE"};
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="MLE");
phottxt=phottxt+", MLE: "+ string(statout.phot(1));

LG=L*2;
sim.definePattern("gauss", psf_vec, phasemask="flat", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=LG,laserpower=laserpower)
sim.defineComponent("estimatorGauss","estimator",@est_GaussLSQ1_2D,parameters={sim.patterns("gauss").pos,LG,110,probecenter},dim=1);
seq={"gauss","estimatorGauss"};
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="Gauss");
phottxt=phottxt+", Gauss: "+ string(statout.phot(1));

title(phottxt)

if ax1v=="pos"
    plot([0 L],[0 L],'k--')
end
ax=gca; ax.YLim(1)=0;ax.YLim(2)=65;
ax.XLim(end)=L*0.75;


%% explore impact of background on estimator
figure(270)
%no background
fl.brightness=2000;
fl.pos=[0 0 0];
sim.background=0;
sim.defineComponent("estdonut","estimator",@est_qLSQiter2D,parameters={L,probecenter},dim=1:2);
seq={"donut","estdonut"};
psf_vec.zerooffset=0;
% xcoords=0:5:50;
nexttile(2); hold off
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,tag="no bg");
photons=statout.phot(1)
% out=sim.runSequence(seq,maxlocs=1);
% sim.summarize_results(out);

%background, 
sim.background=fl.brightness/20;
% xcoords=0:2:100;
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="bg");
photonsbg=statout.phot(1)
%background subtracted,
% bgf=sim.patterns('donut').backgroundfac(1); %background used for simulation
sim.background_estimated=sim.background*laserpower; %in general, the GT background is not known but needs to be calibrated 

%sim.defineComponent("estdonut","estimator",@est_qLSQiter2D,parameters={L,probecenter},dim=1:2);
sim.defineComponent("bg","background",@backgroundsubtractor,parameters={"background_estimated"});
seq={"donut","bg","estdonut"};
psf_vec.zerooffset=0;
% xcoords=0:2:100;
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="bg est",title="estimator bias from background");

sim.background_estimated=0;
sim.defineComponent("estimatorbg","estimator",@est_qDirectFitBg1D,parameters={L},dim=1);
seq={"donut","estimatorbg"};
statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=ax1v,clearfigure=false,tag="bg fit par",title="estimator bias from background");


if ax1v=="pos"
    plot([0 L],[0 L],'k--')
end
ax=gca; ax.YLim(1)=0; ax.YLim(2)=65;
ax.XLim(end)=L*0.75;


%% convergence of the iterative estimator
figure(271); clf

% sim.fluorophores.pos=[30 20 0];
xcoords=0:5:L;
iters=[1,2,3,4, 5, 6,7, 8, 10,  15, 20, 30, 50];
% iters=50
for k=1:length(iters)
    sim.defineComponent("estiter","estimator",@est_qLSQiter2D,parameters={L,probecenter,iters(k)},dim=1:2);
    seq={"donut","estiter"};
    statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true, ax1=ax1v,clearfigure=false,tag="iterations: "+string(iters(k)));
end
statout.phot
plot([0 L],[0 L],'k--')
ax=gca; ax.YLim(1)=0;


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
% %% Fig, 1
% figure(271)
% seq={"donut","estdonut"};
% psf_vec.zerooffset=0;
% sim.background=0;
% xcoords=0:5:50;
% numberOfLocalizations=10000;
% statout=sim.scan_fov(seq,xcoords,"maxlocs",numberOfLocalizations,"display",true,ax1=["rmse","sCRB"],clearfigure=true,tag="simple est");
% ylim([0 12])