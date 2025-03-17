addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist("psf_vecpp","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vecpp=PsfVectorial; %simple 2D donut PSF
end
psf_vecpp.setpinhole("AU",1);
phaseplateposmm=0:0.1:0.5; %mm
zpos=[-200 -100 0 100 200];
phaseplateposrel=phaseplateposmm/2.5; %pupil diameter assumed to be 5 mm; 
stdx=zeros(length(phaseplateposrel),length(zpos),3);crb1=stdx;biasx=stdx;rmsex=stdx;phot=zeros(length(phaseplateposrel),length(zpos),1);
fl=FlStatic(brightness=1000); %define a static fluorophore
fl.pos=[10 0 0];
sim=Simulator(fluorophores=fl); 

L=75; %size of scan pattern
probecenter=true; %should we also probe the center?
orbitpoints=6;
laserpower=5; %relative, increases brightness
pointdwelltime=0.1; % ms, measurement time in each point
repetitions=2; %how often to repeat the pattern scan
sim.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
clear psfall
for k=1:length(phaseplateposrel)
    sys_mis.maskshift=[phaseplateposrel(k),0]; % radius of pupil function is 1
    psf_vecpp.setpar(sys_mis)
    sim.definePattern("donut_misaligned", psf_vecpp, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
        probecenter=probecenter,orbitL=L,pointdwelltime=0.1,laserpower=laserpower,repetitions=1);
    [psfall(:,:,:,k),gridv]=psf_vecpp.imagestack("vortex");
    for z=1:length(zpos)
        sim.fluorophores.pos=[0 0 zpos(z)];
        seq={"donut_misaligned","estdonut"};
        out=sim.runSequence(seq);
        sr=sim.summarize_results(out);
        stdx(k,z,:)=sr.std;
        crb1(k,z,:)=sr.sCRB1;
        biasx(k,z,:)=sr.bias;
        rmsex(k,z,:)=sr.rmse;
        phot(k,z)=sr.phot;
    end
end
%%
zind=3;
stdxrel=stdx./crb1(1,:,:).*sqrt(phot); %normalized to perfectly aligned phaseplate and photon numbers
rmserel=rmsex./crb1(1,:,:).*sqrt(phot); %normalized to perfectly aligned phaseplate and photon numbers
crbrel=crb1./crb1(1,:,:);
crb=crb1./sqrt(phot);
figure(229) %Fig. 1
plot(phaseplateposmm,stdx(:,zind,1),'b',phaseplateposmm,crb(:,zind,1),'b--',phaseplateposmm,rmsex(:,zind,1),'r',phaseplateposmm,-biasx(:,zind,1),'r--');
xlabel('misalignment of phase plate (mm)')
ylabel('std, CRB, RMSE (nm)')
legend('std','CRB','rmse','bias')
ylim([0 40])

figure(226) %ED
subplot(4,3,[1 4])
hold off
for k=1:size(stdxrel,2)
    plot(phaseplateposmm,stdxrel(:,k,1)','Color',colors(k,:));
    hold on
    plot(phaseplateposmm,stdxrel(:,k,2)','--','Color',colors(k,:))
end

% plot(phaseplateposmm,stdx(:,zind,1),'b',phaseplateposmm,stdx(:,zind,2),'b--',phaseplateposmm,crb(:,zind,1),'r',phaseplateposmm,crb(:,zind,2),'r--');
% plot(phaseplateposmm,stdxrel(:,zind,1),'b',phaseplateposmm,stdxrel(:,zind,2),'b--',phaseplateposmm,crbrel(:,zind,1),'r',phaseplateposmm,crbrel(:,zind,2),'r--');
% plot(phaseplateposmm,stdxrel(:,:,1),phaseplateposmm,stdxrel(:,:,2),'--');


xlabel('misalignment of phase plate (mm)')
ylabel('std / CRB aligned')
title("Standard deviation")
legend("x","y")

colors=lines(8);

subplot(4,3,[3 6])
hold off
for k=1:size(biasx,1)
    plot(zpos,biasx(k,:,1),'Color',colors(k,:));
    hold on
    plot(zpos,biasx(k,:,2),'--','Color',colors(k,:))
end
xlabel('z position (nm)')
ylabel('bias (nm)')
pn=vertcat(phaseplateposmm,phaseplateposmm);
legend(string(pn(:)))
title("Bias")

subplot(4,3,[2 5])
hold off
for k=1:size(rmserel,2)
    plot(phaseplateposmm,rmserel(:,k,1)','Color',colors(k,:)) 
    hold on
    plot(phaseplateposmm,rmserel(:,k,2)','--','Color',colors(k,:))
end

% plot(phaseplateposmm,rmserel(:,:,1),phaseplateposmm,rmserel(:,:,2),'--');
xlabel('misalignment of phase plate (mm)')
ylabel('rmse / CRB aligned')
title("Root mean square error (rmse)")
zn=vertcat(zpos,zpos);
legend(string(zn(:)))
%%
clear indz

for z=1:length(zpos)
    indz(z)=find(gridv{3}>=zpos(z),1,"first");
    if z==1
        psfz=psfall(:,:,indz,end);
    else
    psfz=horzcat(psfz,psfall(:,:,indz(z),end));
    end
end
subplot(4,3,7:9)
imagesc(psfz)
axis equal off
title('z pos')

indz0=find(gridv{3}>=0,1,"first");
for k=1:length(phaseplateposmm)
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