classdef PsfVectorial<Psf
    properties
        parameters
        PSFs=dictionary;
        PSFph=[];
        normfactgauss
        pinholepar
        sigma
        fwhm
    end
    methods
        function obj=PsfVectorial(varargin)
            obj@Psf(varargin{:});
            addpath([fileparts(mfilename('fullpath')) '/PSF_simulation/library']);
            obj.parameters=obj.loadparameters('/settings/defaultsystemparameters_vectorialPSF.m');
        end
        function setpar(obj,varargin)
            if nargin >2 %structure passed on
                spar=struct(varargin{:});
            else
                spar=varargin{1};
            end
            % overwrite any of the parameter fields if they exist,
            % otherwise add to addpar
            fn=fieldnames(spar);
            fnpar=fieldnames(obj.parameters);
            recalculate=false;
            for k=1:length(fn)
                assigned=false;
                for l=1:length(fnpar)
                    if isfield(obj.parameters.(fnpar{l}),fn{k})
                        if ~isequal(obj.parameters.(fnpar{l}).(fn{k}),spar.(fn{k}))
                            obj.parameters.(fnpar{l}).(fn{k})=spar.(fn{k});
                            recalculate=true;
                        end
                        assigned=true;
                    end
                end
                if ~assigned
                    obj.parameters.addpar.(fn{k})=spar.(fn{k});
                    % recalculate=true;
                end 
            end
            if recalculate 
                obj.PSFs=dictionary;%delete pre-calculated PSFs
                obj.normfactgauss=[];
                if ~isempty(obj.pinholepar)
                    pp=struct2pairs(obj.pinholepar);
                    obj.setpinhole(pp{:})
                end
            end
        end
        function [idet,phfac]=intensity(obj, flpos ,patternpos, phasepattern, L)
           
            if nargin<5
                key=phasepattern;
            else
                key=phasepattern+L;
            end
            %flpos with respect to optical axis ([0,0] of PSF coordinate
            %system)
            % patternpos: position by EOD only of excitation pattern, no
            % descanning. Pinhole is with respect to [0,0] (but can be
            % shifted in pinhole definition).
            
            psfint=obj.PSFs(key);
            phkey=obj.PSFph;
            flposrel=flpos-patternpos;
            iexc=psfint.interp(flposrel)+obj.zerooffset; 
            iexc=max(iexc,0);
            if isempty(phkey) 
                phfac=ones(size(iexc));
            else
                psfph=obj.PSFs(phkey);
                phfac=psfph.interp(flpos);
            end
            idet=iexc.*phfac;
        end
        function key=calculatePSFs(obj,phasepattern,Lxs,forcecalculation)
            if nargin<4
                forcecalculation=false;
            end
            if isnumeric(Lxs)
                Lx=Lxs;
                Lxs=num2str(Lxs);
            else
                Lx=str2num(Lxs);
            end
            key=string(phasepattern)+Lxs;
            if obj.PSFs.isConfigured && obj.PSFs.isKey(key) && ~forcecalculation
                return
            end
            fprintf(key + ", ")
            intmethod='cubic'; %linear,cubic?
            extraolation_method='nearest'; 

            opt=obj.parameters.opt;
            out=obj.parameters.addpar;
             sys=obj.parameters.sys;
            % interpolant:
            % nx=(-opt.radiusCanvas:out.dr:opt.radiusCanvas)*1e9; %to nm
            % nz=(-opt.depthCanvas:out.dz:opt.depthCanvas)*1e9; %to nm
            % [X,Y,Z] = ndgrid(nx,nx,nz);
            
           
            % sigma from wavelength, NA
            fwhm=0.51*sys.loem/sys.NA*1e9;
            obj.sigma=fwhm/2.35;
            
            if isempty(obj.normfactgauss)
                sys.Ei = {'circular'};
                outc=effField(sys,out, opt);      
                outc=effIntensity(sys,outc);
                zmid=ceil(size(outc.I,4)/2);
                % normfact=sum(sum(outc.I(1,:,:,zmid)));%normalized to integral =1
                normfact=max(max(outc.I(1,:,:,zmid)));%normalized to max =1
                obj.normfactgauss=normfact;
            else
                normfact=obj.normfactgauss;
            end
                
            % %fit Gaussian:
            % vg=squeeze(outc.I(1,:,ceil(end/2),zmid))';
            % x=1:length(vg); x=x'-mean(x); x=x*out.dr*1e9;
            % fitp=fit(x,vg,'gauss1');
            % obj.sigma=fitp.c1/sqrt(2);

            switch phasepattern
                case 'flat'
                    sys.Ei = {'circular'};    
                    sys=addzernikeaberrations(sys,out);
                    out=effField(sys,out, opt);      
                    out=effIntensity(sys,out);
                    PSF=squeeze(out.I(:,:,:,:))/normfact;
                    PSF = obj.beadsize(PSF, sys.beadradius);
                    [PSF,PSFdonut.normalization]=normpsf(PSF);
                    [X,Y,Z]=meshgrid4PSF(PSF,out.dr,out.dz);
                    PSFdonut.interp = griddedInterpolant(X,Y,Z,PSF,intmethod, extraolation_method);              
                    obj.PSFs(key)=PSFdonut;
                case {'halfmoonx','halfmoony'}
                    %convert L into phase
                    dxdphi=1.4; %nm/degree, from LSA paper Deguchi
                    % dzdphi=-3.6; %nm/degree
                    sys.delshift = deg2rad(Lx/dxdphi);
                    sys.Ei = {'halfmoon', 'linear'};  
                    sys=addzernikeaberrations(sys,out);
                    out=effField(sys,out, opt);      
                    out=effIntensity(sys,out);
                    PSF=squeeze(out.I(:,:,:,:))/normfact;
                    PSF = obj.beadsize(PSF, sys.beadradius);
                    [PSF,normalization]=normpsf(PSF);
                    [X,Y,Z]=meshgrid4PSF(PSF,out.dr,out.dz);
                    PSFx.interp = griddedInterpolant(X,Y,Z,permute(PSF,[2,1,3]),intmethod, extraolation_method);
                    PSFy.interp = griddedInterpolant(X,Y,Z,PSF,intmethod, extraolation_method);     
                    PSFx.normalization=normalization;
                    PSFy.normalization=normalization;
                    obj.PSFs("halfmoonx"+Lxs) =PSFx; 
                    obj.PSFs("halfmoony"+Lxs)=PSFy; 
                case 'pinhole' %here
                    sys.Ei = {'circular'};
                    sys=addzernikeaberrations(sys,out);
                    phdiameter=Lx(1);
                    if length(Lx)==3
                        phpos=[Lx(2) Lx(3)];
                    else
                        phpos=[0 0];
                    end
                    out=effField(sys,out, opt);      
                    out=effIntensity(sys,out);
                    psfg=squeeze(out.I);
                    % psfgmaxnorm=psfg/normfact;
                    norm=sum(sum(psfg(:,:,round(end/2))));
                    psfg=psfg/norm;
                    pixelsize=opt.pixSize*1e9;
                    nx=1:size(psfg,1);
                    nx=(nx-mean(nx))*pixelsize;
                    [Xk,Yk]=meshgrid(nx);
                    kernel=double((Xk-phpos(1)).^2+(Yk-phpos(2)).^2<(phdiameter/2)^2);
                    psfph=conv2fft(psfg,kernel);
                    [psfph,PSFdonut.normalization]=normpsf(psfph);
                    [X,Y,Z]=meshgrid4PSF(psfph,out.dr,out.dz);
                    PSFdonut.interp = griddedInterpolant(X,Y,Z,psfph,intmethod,extraolation_method);
                    obj.PSFs(key)=PSFdonut;   
                case 'tophat'
                    dzdphi=-3.6; %nm/degree
                    sys.delshift = deg2rad(Lx/dzdphi);
                    sys.Ei = {'pishift', 'circular'};  
                    sys=addzernikeaberrations(sys,out);
                    out=effField(sys,out, opt);      
                    out=effIntensity(sys,out);
                    PSF=squeeze(out.I(:,:,:,:))/normfact;
                    PSF = obj.beadsize(PSF, sys.beadradius);
                    [PSF,PSFdonut.normalization]=normpsf(PSF);
                    [X,Y,Z]=meshgrid4PSF(PSF,out.dr,out.dz);
                    PSFdonut.interp = griddedInterpolant(X,Y,Z,PSF,intmethod, extraolation_method);
                    obj.PSFs(key)=PSFdonut;   
                case 'vortex'
                    sys.Ei = { 'phaseramp',  'circular'};  
                    sys=addzernikeaberrations(sys,out);
                    out=effField(sys,out, opt);      
                    out=effIntensity(sys,out);
                    PSF=squeeze(out.I(:,:,:,:))/normfact;
                    PSF = obj.beadsize(PSF, sys.beadradius);
                    [PSF,PSFdonut.normalization]=normpsf(PSF);
                    [X,Y,Z]=meshgrid4PSF(PSF,out.dr,out.dz);
                    PSFdonut.interp = griddedInterpolant(X,Y,Z,PSF,intmethod, extraolation_method);              
                    obj.PSFs(key)=PSFdonut;
                otherwise
                    warning(phasepattern+" PSF name not defined")
            end
        end
        function psfo=beadsize(obj,psf,R)
            if R==0
                psfo=psf;
                return
            end
            dr=obj.parameters.addpar.dr; dz=obj.parameters.addpar.dz;
            nx=1:size(psf,1);
            nx=(nx-mean(nx))*dr;
            ny=1:size(psf,2);
            ny=(ny-mean(ny))*dr;
            nz=1:size(psf,3);
            nz=(nz-mean(nz))*dz;
            [X,Y,Z]=meshgrid(nx,ny,nz);
            circpsf=double((X).^2+(Y).^2+(Z).^2<=R^2);
            circpsf=circpsf/sum(circpsf(:));
            psfo=convnfft(psf,circpsf);
            % psfo=ifftshift(ifftshift(ifftshift(psfo,1),2),3); %inconsitent use of fftshift, but ok because bead is symmetric
        end
        function setpinhole(obj, args)
            arguments
                obj 
                args.lambda=600; %nm
                args.AU =1; %in airy units
                args.diameter=[] %nm
                args.offset=[0,0];
            end
            if length(args.offset)~=2
                error('PSF_vectorial.setpinhole: args.offset needs to have two entries (x,y)')
            end
            if isempty(args.diameter)
                % [sys,opt,out]=systemparameters;
                args.diameter=round(args.AU*1.22*obj.parameters.sys.loem*1e9/obj.parameters.sys.NA); %single nm accuracy should be sufficient
            end

            %sdetermine if PH is already calculated
            phdefined=~isempty(obj.PSFph) && ~isempty(obj.pinholepar);
            notchanged=phdefined && (args.diameter == obj.pinholepar.diameter...
                    && all(args.offset == obj.pinholepar.offset)); 

            calculated=notchanged && ((obj.PSFs.isConfigured) && isKey(obj.PSFs,"pinhole"+num2str([args.diameter,args.offset])));

            % if isempty(obj.PSFph)||(~isempty(obj.pinholepar) && (args.diameter ~= obj.pinholepar.diameter...
            %         || any(args.offset ~= obj.pinholepar.offset)) ) || (~isempty(obj.PSFs) && ~isKey(obj.PSFs,"pinhole"+num2str([args.diameter,args.offset])))
            if ~calculated
                obj.PSFph=obj.calculatePSFs('pinhole',[args.diameter,args.offset],1);
            end
            obj.pinholepar.offset=args.offset;
            obj.pinholepar.diameter=args.diameter;
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
        function showpsf(obj,key)
            psf=obj.PSFs(key);
            inten=psf.Values;
            if isempty(obj.PSFph)
                psfph=1;
            else
                psfphx=obj.PSFs(obj.PSFph);
                psfph=psfphx.Values;
            end
            psftot=inten.*psfph;
            imx(psftot);
        end
        function par=loadparameters(obj,filen)
            [fpath,fname,ext]=fileparts(filen);
            switch ext
                case '.m'
                    if exist(fpath,"dir")
                        addpath(fpath)
                    end
                    [sys,opt,out]=eval(fname);
                    [sys2,out]=effInit_oil_exc(sys,out,opt);
                    sys=copyfields(sys,sys2);
                    par.sys=sys;
                    par.opt=opt;
                    par.addpar=out;
                case '.json'
                    disp('json not implemented yet')
            end
        end
        % function out=fwhm(obj)
        %      out=obj.sigma*2.35*1.2; %comes from comparison with simple donut
        % end
        function loadexperimental(obj,key,fname)
            [~,fnh,ext]=fileparts(fname);
            if contains(fname,"_3dcal") && ext==".mat" %3Dcal.mat
                l=load(fname);
                PSFvol=l.SXY(1).PSF{1};
                px=l.parameters.pixelsize{1}(1)*1000; %nm
                pz=l.parameters.dz*1000;
                spsf=size(PSFvol);
                nx=1:spsf(1);nx=nx-mean(nx);nx=nx*px;
                ny=1:spsf(2);ny=ny-mean(ny);ny=ny*px;
                nz=1:spsf(3);nz=nz-mean(nz);nz=nz*pz;
                [X,Y,Z] = ndgrid(nx,nx,nz);
                intmethod='cubic'; %linear,cubic?
                PSFexp = griddedInterpolant(X,Y,Z,PSFvol,intmethod);              
                obj.PSFs(key+"0")=PSFexp;
            

            else
                disp('not implemented')
            end
            %uiPSF
            %tiffstack
        end
        function vout=imagestack(obj,key,args)
            arguments
                obj
                key
                args.show=false;
            end
            if ~isKey(obj.PSFs,key)
                key=key+"0";
            end
            psf=obj.PSFs(key);
            vout=psf.interp.Values;
            if ~isempty(obj.PSFph)
                phpsf=obj.PSFs(obj.PSFph);
                vout=vout.*phpsf.interp.Values  ; 
            end

            if args.show
                imx(double(vout))
            end
        end
        function nf=normfactor(obj,phasemask,number)
            nf=obj.PSFs(string(phasemask)+number).normalization;
        end
    end
end

function sys=addzernikeaberrations(sys,addpar)
if isfield(sys,'Zr') && ~isempty(sys.Zr)
    sys.Ei{end+1}='zernike';
elseif isfield(addpar,'Zr') && ~isempty(addpar.Zr)
    sys.Ei{end+1}='zernike';
    sys.Zr=addpar.Zr;
end
end


function [PSF,normf]=normpsf(PSF)
normf=max(PSF(:));
PSF=PSF/normf;
end

function [X,Y,Z]=meshgrid4PSF(PSF,dr,dz)
nx=((1:size(PSF,1))-mean((1:size(PSF,1))))*dr*1e9;
nz=((1:size(PSF,3))-mean((1:size(PSF,3))))*dz*1e9;
[X,Y,Z]=ndgrid(nx,nx,nz);
end