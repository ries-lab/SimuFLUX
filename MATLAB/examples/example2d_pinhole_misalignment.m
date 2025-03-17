if ~exist("psf_vecpp","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vecph=PsfVectorial; %simple 2D donut PSF
end
pinholepos=0:100:500; %um

stdx=zeros(length(pinholepos),3);crb1=stdx;biasx=stdx;rmsex=stdx;phot=zeros(length(pinholepos),1);
fl=FlStatic(brightness=1000); %define a static fluorophore
fl.pos=[0 0 0];
sim=Simulator(fluorophores=fl); 

L=75; %size of scan pattern
probecenter=true; %should we also probe the center?
orbitpoints=6;
laserpower=5; %relative, increases brightness
pointdwelltime=0.1; % ms, measurement time in each point
repetitions=1; %how often to repeat the pattern scan
sim.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
clear psfall
for k=1:length(pinholepos)
    psf_vecph.setpinhole("AU",1,"offset",[pinholepos(k) 0]);
    sim.definePattern("ph_misaligned", psf_vecph, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
        probecenter=probecenter,orbitL=L,pointdwelltime=0.1,laserpower=laserpower,repetitions=repetitions);
       [psfall(:,:,:,k),gridv]=psf_vecph.imagestack("vortex");
  
        seq={"ph_misaligned","estdonut"};
        out=sim.runSequence(seq,"maxlocs",10000);
        sr=sim.summarize_results(out);
        stdx(k,:)=sr.std;
        crb1(k,:)=sr.sCRB1;
        biasx(k,:)=sr.bias;
        rmsex(k,:)=sr.rmse;
        phot(k)=sr.phot;
end
%%
zind=3;
stdxrel=stdx./crb1(1,:,:).*sqrt(phot); %normalized to perfectly aligned phaseplate and photon numbers
rmserel=rmsex./crb1(1,:,:).*sqrt(phot); %normalized to perfectly aligned phaseplate and photon numbers
crbrel=crb1./crb1(1,:,:);
crb=crb1./sqrt(phot);


figure(228) %ED
subplot(4,3,[1 4])
hold off

    plot(pinholepos,stdxrel(:,1)','k');
    hold on
    plot(pinholepos,stdxrel(:,2)','k--')

% plot(phaseplateposmm,stdx(:,zind,1),'b',phaseplateposmm,stdx(:,zind,2),'b--',phaseplateposmm,crb(:,zind,1),'r',phaseplateposmm,crb(:,zind,2),'r--');
% plot(phaseplateposmm,stdxrel(:,zind,1),'b',phaseplateposmm,stdxrel(:,zind,2),'b--',phaseplateposmm,crbrel(:,zind,1),'r',phaseplateposmm,crbrel(:,zind,2),'r--');
% plot(phaseplateposmm,stdxrel(:,:,1),phaseplateposmm,stdxrel(:,:,2),'--');


xlabel('misalignment of pinhole (nm)')
ylabel('std / CRB aligned')
title("Standard deviation")
legend("x","y")
ylim([0.9 1.1])

colors=lines(8);

subplot(4,3,[3 6])

    plot(pinholepos,biasx(:,1),'k');
    hold on
    plot(pinholepos,biasx(:,2),'k--')

xlabel('misalignment of pinhole (nm)')
ylabel('bias (nm)')
ylim([-1 1])

title("Bias")

subplot(4,3,[2 5])
hold off

    plot(pinholepos,rmserel(:,1)','k') 
    hold on
    plot(pinholepos,rmserel(:,2)','k--')


% plot(phaseplateposmm,rmserel(:,:,1),phaseplateposmm,rmserel(:,:,2),'--');
xlabel('misalignment of pinhole (nm)')
ylabel('rmse / CRB aligned')
title("Root mean square error (rmse)")
ylim([0.9 1.1])
%%
clear indz

indz0=find(gridv{3}>=0,1,"first");
for k=1:length(pinholepos)
    if k==1
        psfx=psfall(:,:,indz0,k);
    else
    psfx=horzcat(psfx,psfall(:,:,indz0,k));
    end
end
subplot(4,3,10:12)
imagesc(psfx)
axis equal off
title('misalignment')

% plot(phaseplateposmm,stdxrel(:,:,1),phaseplateposmm,stdxrel(:,:,2))
% disp("misaligned phase plate:")
% sim.summarize_results(out);
% psfab=psf_vec2.imagestack("vortex");
% imx(horzcat(psf0,psfab),'Parent',figure(221),'Title',"Misaligned phase plate"); %compare the two PSFs
% 