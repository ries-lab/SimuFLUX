% Simulate several fluorophores in a fluorophore collection
%% two fluorophores
fl1=Fl_static;
fl2=Fl_static;
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PSF_vectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);

fc=Flcollection;
fc.add([fl1, fl2])

sim=Sim_Simulator(fc); %make a simulator and attach fluorophore

numberOfLocalizations=1000;

L=75;
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=4, ...
    probecenter=false,orbitL=L,laserpower=100)
sim.defineComponent("estdonut","estimator",@estimators,parameters={"donut2D",sim.patterns("donut").pos,L,360},dim=1:2);


seq={"donut","estdonut"};
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
stats=sim.displayresults(out);

x=0:20:750;
biasx=zeros(length(x),3); stdx=biasx; rmsex=biasx;
for k=1:length(x)
    fl2.pos(1)=x(k);
    out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
    stats=sim.displayresults(out);
    biasx(k,:)=stats.bias;
    stdx(k,:)=stats.std;
    rmsex(k,:)=stats.rmse;
end

figure(145)
yyaxis right
plot(x,biasx(:,1))
ylabel('bias (nm)')
yyaxis left
plot(x,rmsex(:,1),x,stdx(:,1))
legend('rmse','std','bias')
ylabel('rmse / std (nm)')
xlabel('separation x (nm)')


%% Imaging of blinking fluorophores with Abberior sequence
%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"Sim_sequencefile")
    sim=Sim_sequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
end
fname='Imaging_2D.json';
sim.loadsequence(fname); %only sequence file, then simple gauss and donut PSFs are used (fast)
sim.makepatterns;
sim.makescoutingpattern([-100 -100; 400 250 ]) %for imaging

% make fake NPCs
clear posnpc;
R=50;dphi=pi/4;
phi=0:dphi:2*pi-dphi; posnpc(:,1)=R*cos(phi); posnpc(:,2)= R*sin(phi);

% make a fluorophore collection with blinking fluorophores
fc=Flcollection_blinking;

%set parameterst for caged fluorophore, PAFP or similar
laserpower=5;
switchpar.brightness=100*laserpower;
switchpar.toffsmlm=3*1e6; %on-switching time in ns
switchpar.photonbudget=5000;
switchpar.tonsmlm=1e8; % ns stays on, only bleached
switchpar.activations=5; %re activations
switchpar.starton=0; %fluorophores start in random on / off state, determined by tonsmlm, toffsmlm
fc.setpar(switchpar)

%add fake NPCs
fc.add(posnpc);
% add more NPCs at positions dpos
dpos=[250 50 ]; fc.addstatic(posnpc+dpos);
dpos=[50 150 ]; fc.addstatic(posnpc+dpos);

sim.fluorophores=fc;
out=sim.scoutingSequence(maxrep=1000);

%plot results
vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
vldcfr=vld & out.loc.cfr<0.1;
figure(89); hold off;
plot(sim.scoutingcoordinates(:,1),sim.scoutingcoordinates(:,2),'k+')
hold on
plot(out.loc.xnm(vld),out.loc.ynm(vld),'r.')
plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
posfl=vertcat(fc.flall.pos);
plot(posfl(:,1),posfl(:,2),'mo')
axis equal
legend('scouting','last itr vld','last itr vld +cfr', 'fluorophore')