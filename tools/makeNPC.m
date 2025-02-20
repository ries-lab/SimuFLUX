function posnpc=makeNPC(args)
arguments
    args.pos=[0,0,0];
    args.R=50;
    args.copynumber=32;
    args.dz=50;
    args.rotation=0; % in plane rotation
    args.twistangle=pi/32;
    args.shiftangle=pi/16;
    args.dR=3;
end
dphi=pi/4;
phi=0:dphi:2*pi-dphi-args.rotation; 
v0=0*phi;
switch args.copynumber
    case 8
        phiall=phi;
        zall=v0;
        Rall=v0+args.R;
    case 16
        phiall=horzcat(phi,-phi+args.twistangle);
        zall=horzcat(v0+args.dz/2,v0-args.dz/2);
        Rall=horzcat(v0,v0)+args.R;
    case 32
        phiall=horzcat(phi,phi+args.shiftangle,-phi+args.twistangle,-phi+args.twistangle-args.shiftangle);
        zall=horzcat(v0+args.dz/2,v0+args.dz/2,v0-args.dz/2,v0-args.dz/2);
        Rall=horzcat(v0+args.R+args.dR,v0+args.R,v0+args.R+args.dR,v0+args.R);
    otherwise
        disp('NPC copy number: 8, 16 or 32')
        pos=[];
        return
end

posnpc(:,1)= Rall.*cos(phiall); posnpc(:,2)= Rall.*sin(phiall);
posnpc(:,3)= zall;
posnpc=posnpc+args.pos;

