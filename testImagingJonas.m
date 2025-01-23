% testImagingJonas
clear pos
R=50; phi=0:pi/4:2*pi; pos(:,1)=R*cos(phi); pos(:,2)= R*sin(phi);
fc=Fluorophorecollection;
fc.switchpar.brightnes=1000;
fc.switchpar.toffsmlm=100000;
fc.switchpar.photonbudget=500000;
fc.addstatic(pos);

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
vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) & out.loc.cfr<0.2;
figure(88); plot(out.loc.xnm(vld),out.loc.ynm(vld),'+')
axis equal