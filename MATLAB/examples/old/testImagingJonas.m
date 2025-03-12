% testImagingJonas



% testImagingJonas
fl=bleachingFluorophore;fl.pos=[0 0 0];
fl.photonbudget=10000;
fl2=bleachingFluorophore;fl2.pos=[600 0 0];fl2.photonbudget=10000;
fc=Fluorophorecollection;
fc.add([fl fl2]);

sim=AbberiorSimulator(fl);
fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit100__locLim50_.json';
sim.loadsequence(fname)

sim.makepatterns
sim.makescoutingpattern([-200 -200; 1200 1200 ])
% out=sim.runSequence;
out=sim.scoutingSequence(maxrep=1);
% out=sim.repeatSequence('',1,5000);
% figure(88);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
% plot(out.loc.loccounter,out.loc.xeod)
% xlabel('time (itr)')
% legend('xest', 'xfl','xgalvo','EOD')
vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
vldcfr=vld & out.loc.cfr<0.3;
figure(89); hold off;
plot(out.loc.xnm(~vldcfr),out.loc.ynm(~vldcfr),'r.')
hold on
plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'b+')
axis equal
asfdasd
%%

clear pos
R=50; phi=0:pi/4:2*pi; pos(:,1)=R*cos(phi); pos(:,2)= R*sin(phi);
fc=blinkingFluorophorecollection;

switchpar.brightness=1000;
switchpar.toffsmlm=100000;
switchpar.photonbudget=1000;
switchpar.tonsmlm=1e8; %stays on, only bleached
switchpar.activations=1;
switchpar.starton=0;
fc.setpar(switchpar)

% fc.setpar('brightness',100)

fc.add(pos);

dpos=[500 0 ];
fc.addstatic(pos+dpos);

dpos=[0 1000 ];
fc.addstatic(pos+dpos);

% fl.brightness=1000;
% fl.photonbudget=20000;
psfdonut=PSFMF_donut2D;
sim=AbberiorSimulator(fc);

numlocs=12;
% fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit1000_.json';
fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit100__locLim50_.json';
sim.loadsequence(fname)
% sim.sequence_json.damping=.1;
sim.makepatterns
sim.makescoutingpattern([-200 -200; 1200 1200 ])
% out=sim.runSequence;
out=sim.scoutingSequence(maxrep=100);
% out=sim.repeatSequence('',1,5000);
% figure(88);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
% plot(out.loc.loccounter,out.loc.xeod)
% xlabel('time (itr)')
% legend('xest', 'xfl','xgalvo','EOD')
vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
vldcfr=vld & out.loc.cfr<0.3;
figure(89); hold off;
plot(out.loc.xnm(vld),out.loc.ynm(vld),'r.')
hold on
plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'b+')
axis equal