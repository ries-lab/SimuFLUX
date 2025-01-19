classdef PSFMF_realistic<PSFMF
    properties
        sigmaz=0; %don't do z sectioning.
    end
    methods
        function [io,iaf]=intensity(obj, flposrel ,phasepattern, L)
            key=phasepattern+L;
            psfint=obj.PSFs(key);
            sigmaz=obj.sigmaz;
            if sigmaz>0
                zfac=exp(-flposrel(:,3).^2/2/sigmaz^2);
            else 
                zfac=1;
            end
            io=psfint(flposrel)*zfac+obj.zerooffset;          
        end
        function calculatePSFs(obj,phasepattern,Lxs)
            key=phasepattern+Lxs;
            if obj.PSFs.isConfigured && obj.PSFs.isKey(key)
                return
            end
            fprintf(key + ", ")
            addpath('/PSF_simulation/library');
            [sys,opt,out]=systemparameters;
            intmethod='cubic'; %linear,cubic?
            Lx=str2double(Lxs);
            %convert L into phase
            dxdphi=1.4; %nm/degree, from LSA paper Deguchi
            dzdphi=-3.6; %nm/degree
            

            % interpolant:
            nx=(-opt.radiusCanvas:opt.pixSize:opt.radiusCanvas)*1e9; %to nm
            nz=(-opt.depthCanvas:out.dz:opt.depthCanvas)*1e9; %to nm
            [X,Y,Z] = ndgrid(nx,nx,nz);

            %Gaussian for normalization
            sys.Ei = {'circular'};
            outc=out;
            outc.dz=opt.depthCanvas; % Axial range. 
            outc=effField(sys,outc, opt);      
            outc=effIntensity(sys,outc);
            zmid=ceil(size(outc.I,4)/2);
            % normfact=sum(sum(outc.I(1,:,:,zmid)));%normalized to integral =1
            normfact=max(max(outc.I(1,:,:,zmid)));%normalized to max =1
            switch phasepattern
                case {'x','y'}
                    sys.delshift = deg2rad(Lx/dxdphi);
                    sys.Ei = {'halfmoon', 'linear'};   
                    out=effField(sys,out, opt);      
                    out=effIntensity(sys,out);
                    PSF=squeeze(out.I(:,:,:,:))/normfact;
                    PSFx = griddedInterpolant(X,Y,Z,permute(PSF,[2,1,3]),intmethod);
                    PSFy = griddedInterpolant(X,Y,Z,PSF,intmethod);     
                    obj.PSFs("x"+Lxs)=PSFx; 
                    obj.PSFs("y"+Lxs)=PSFy; 
                case 'tophat'
                    sys.delshift = deg2rad(Lx/dzdphi);
                    sys.Ei = {'pishift', 'circular'};  
                    out=effField(sys,out, opt);      
                    out=effIntensity(sys,out);
                    PSF=squeeze(out.I(:,:,:,:))/normfact;
                    PSFdonut = griddedInterpolant(X,Y,Z,PSF,intmethod);
                    obj.PSFs(key)=PSFdonut;   
                case 'vortex'
                    sys.Ei = { 'phaseramp',  'circular'};    
                    out=effField(sys,out, opt);      
                    out=effIntensity(sys,out);
                    PSF=squeeze(out.I(:,:,:,:))/normfact;
                    PSFdonut = griddedInterpolant(X,Y,Z,PSF,intmethod);
                    obj.PSFs(key)=PSFdonut;
                otherwise
                    xxx
            end
        end
        function setpinhole(obj,args)
            arguments
                obj 
                args.lambda=600; %nm
                args.AU =1; %in airy units
                args.diameter=[] %nm
                args.NA=1.5;
                args.refractiveIndex=1.5; %1.33 water, %1.5 oil
            end
            if isempty(args.diameter)
                args.diameter=args.AU*1.22*args.lambda/args.NA;
            end
            n=args.refractiveIndex;
            fwhm2=(0.88*args.lambda/(n-sqrt(n^2-args.NA^2)))^2+(sqrt(2)*n*args.diameter/args.NA)^2;
            obj.sigmaz=sqrt(fwhm2)/2.35;
        end
        function savePSF(obj,name)
            PSFinterpolant=obj.PSFinterpolant;
            save(name,'PSFinterpolant')
        end
        function loadPSF(obj,name)
            load(name,'PSFinterpolant');
            obj.PSFinterpolant=PSFinterpolant;
        end
        function plotprofile(obj,key)
            psf=obj.PSFs(key);
            inten=psf.Values;
            xaxv=(psf.GridVectors{1});
            figure(33)
            plot(xaxv, inten(:,round(end/2),round(end/2)))
        end
    end
end

function [sys,opt,out]=systemparameters
sys = [];
out = [];    
opt.Et = 0; % Display Et during the calculation.
opt.Ef = 0;    % Display Ef during the calculation.
opt.calbar = 0;    % Display calculation progress bar.
opt.mem = 50000;
opt.pixSize = 5e-9; % Although the data size will be huge, pixel size of 2e-9 is recommended for fine calculation.
opt.radiusCanvas = 0.7e-6; % Lateral range.
opt.depthCanvas = 0.8e-6; % Axial range. 
opt.polAngle = 0.0*pi; % Linear polarization angle.
opt.phaseImage = false; % Show the phase image at the back focal plane.
opt.intImage = false;

% Setting specific parameters
[sys,out]=effInit_oil_exc(sys,out,opt);        % Assigning initial parameters. 
sys.NA=1.35;
sys.nm=1.406;
sys.ns=1.406;
sys.Mt=100;
sys.wa=7e-3;
sys.Na=200; 
sys.phaseImage = opt.phaseImage;
sys.intImage = opt.intImage;
sys.wa=7e-3;
sys.rz = 0.8e-6;% Axial range.
sys.pl = 0; % Angle of the linear polarization.
end