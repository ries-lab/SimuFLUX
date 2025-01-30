classdef Sim_Simulator<handle
    properties
        patterns=dictionary;
        sequences=dictionary;
        fluorophores
        background=0; %kHz from AF, does not count towards photon budget
        posgalvo=[0 0 0]; %nm, position of pattern center
        posEOD=[0 0 0]; %nm, not descanned, with respect to posgalvo.
        time=0;
    end
    methods
        function obj=Sim_Simulator(fl)
            if nargin>0    
                obj.fluorophores=fl;
            end
        end
        function defineComponent(obj,key,type,functionhandle,args)
            %define estimators, recenterers or other functions to be called
            %within a pattern repeat
            arguments 
                obj
                key
                type
                functionhandle
                args.parameters={};
                args.dim=1:3;
            end
            component.functionhandle=functionhandle;
            component.type=string(type);
            component.dim=args.dim;
            component.parameters=args.parameters;
            obj.patterns(key)=component;
        end
        function definePattern(obj,key,psf,args)
            %defines a pattern based on a PSF that is passed on and other
            %parameters
            arguments 
                obj
                key
                psf
                args.phasemask=""; %parameter that defines the shape of the phase mask
                args.zeropos=0; %position of the zero when calculating PSFs (e.g., PhaseFLUX)
                args.patternpos=[0, 0, 0];% list of 3D positions where the PSF is moved to. Use either zeropos or pos
                args.makepattern=[];
                args.orbitpoints=4;
                args.orbitL=100;
                args.probecenter=true; 
                args.orbitorder=[]; %change the order of measurement points
                args.pointdwelltime=10; %us
                args.laserpower=1; %usually we use relative, but can also be absolute.
                args.repetitions=1;
                args.dim=1:2;
            end
            pattern.repetitions=args.repetitions;
            %make a list of either positions of the PSF or of zeropositions
            %to be calculated.
            pos=args.patternpos;
            zeropos=args.zeropos;
            if ~isempty(args.makepattern)
                posh=obj.makeorbitpattern(args.makepattern,args.orbitpoints,args.probecenter,args.orbitorder);
                pos=posh.*args.orbitL/2;
            end
            if size(pos,1)==1 && length(zeropos)>1
                pos=repmat(pos,length(zeropos),1);
            end
            if size(pos,1)>1 && length(zeropos)<=1
                zeropos=repmat(zeropos,1,size(pos,1));
            end
            pattern.pos=pos;
            
            %create values for all measurement points of the pattern.
            for k=1:size(pos,1)
                pattern.psf(k)=psf;
                pattern.phasemask(k)=string(args.phasemask);
                pattern.zeropos(k)=string((zeropos(k)));
                pattern.psf(k).calculatePSFs(args.phasemask,pattern.zeropos(k));
                pattern.pointdwelltime(k)=args.pointdwelltime;
                pattern.laserpower(k)=args.laserpower;
            end
            pattern.L=args.orbitL;
            pattern.dim=args.dim;
            pattern.type="pattern";
            obj.patterns(key)=pattern;
        end
       

        function out=patternscan(obj,key)
            % scans a defined pattern
            pattern=obj.patterns(key);
            fluorophores=obj.fluorophores;
            numpoints=length(pattern.zeropos);
            intall=zeros(numpoints,1);
            repetitions=pattern.repetitions;
            timestart=obj.time;
           
            timep=0;
            posgalvo=obj.posgalvo;
            posEOD=obj.posEOD;
            flpos=zeros(fluorophores.numberOfFluorophores,3);
            flintall=zeros(fluorophores.numberOfFluorophores,1);
            flproperties=fluorophores.getproperties; %for performance
            time=obj.time;
            for r=1:repetitions
                for k=1:numpoints
                    timep=timep+time; %for calculating average time point
                    [flposh,isactive]=fluorophores.position(obj.time,flproperties);
                    flposrel=flposh-posgalvo;
                    [intensityh,pinholehfac]=pattern.psf(k).intensity(flposrel,pattern.pos(k,:)+posEOD,pattern.phasemask(k),pattern.zeropos(k));
                    intensityh=intensityh*pattern.laserpower(k);
                    flint=fluorophores.intensity(intensityh,pattern.pointdwelltime(k),time,pinholehfac,flproperties);
                    intensity=sum(flint);
                    flpos(isactive,:)=flpos(isactive,:)+flposh;
                    flintall(isactive,:)=flintall(isactive,:)+flint;
                    time=time+pattern.pointdwelltime(k);
                    intall(k)=intall(k)+intensity+obj.background; %sum over repetitions, fluorophores
                end
                fluorophores.updateonoff(time);
            end
            out.phot=poissrnd(intall); %later: fl.tophot(intenall): adds bg, multiplies with brightness, does 
            out.photbg=obj.background*sum(pattern.pointdwelltime);
            out.flpos=flpos/numpoints/repetitions;
            out.flint=flintall;
            out.averagetime=timep/numpoints/repetitions;
            out.patterntotaltime=time-timestart;
            out.counter=1;
            out.repetitions=repetitions;
            obj.time=time;
            
        end

        function out=runSequenceintern(obj,seq,maxlocalizations)
            numseq=length(seq);
            numpat=0;
            for s=numseq:-1:1
                pat=obj.patterns(seq{s});
                if strcmp(pat.type,"pattern")
                    ls(s)=length(pat.zeropos);
                    numpat=numpat+1;
                end
            end
            photch=zeros(maxlocalizations*numpat,max(ls))-1;
            loc=initlocs(maxlocalizations*numpat,{'xnm','ynm','znm','xfl1',...
                'yfl1','zfl1','phot','time','abortcondition', 'xgalvo', 'ygalvo', 'zgalvo',...
                'xeod','yeod','zeod'});
   
            bleached=false;
            loccounter=0;
            for k=1:maxlocalizations
                if obj.fluorophores.remainingphotons<1
                    bleached=true;
                    break
                end
                posgalvo_beforecenter=obj.posgalvo;
                xest=[0,0,0];
                for s=1:numseq
                    component=obj.patterns(seq{s});
                    switch component.type
                        case "pattern"
                            scanout=obj.patternscan(seq{s});
                            loccounter=loccounter+1; % every localization gets a new entry, as in abberior
                            loc.phot(loccounter)=loc.phot(loccounter)+sum(scanout.phot);
                            photch(loccounter,1:length(scanout.phot))=scanout.phot;
                            loc.time(loccounter)=loc.time(loccounter)+scanout.averagetime;
                            flpos=scanout.flpos(1,:);
                            loc.xfl1(loccounter)=flpos(1);loc.yfl1(loccounter)=flpos(2);loc.zfl1(loccounter)=flpos(3);
                            fluorophores.pos(loccounter,1:size(scanout.flpos,1),:)=scanout.flpos;
                            fluorophores.int(loccounter,1:size(scanout.flpos,1))=scanout.flint;
                        case "estimator"
                            % par=replaceinlist(component.parameters,'patternpos',patternpos,'L',L,'probecenter',probecenter);
                            xesth=component.functionhandle(scanout.phot,component.parameters{:});
                            if length(xesth)==3
                                xest(component.dim)=xesth(component.dim);
                            else
                                xest(component.dim)=xesth;
                            end
                            xesttot=xest+posgalvo_beforecenter;
                            loc.xnm(loccounter)=xesttot(1);loc.ynm(loccounter)=xesttot(2);loc.znm(loccounter)=xesttot(3);
                            loc.xgalvo(loccounter,1)=obj.posgalvo(1);loc.ygalvo(loccounter,1)=obj.posgalvo(2);loc.zgalvo(loccounter,1)=obj.posgalvo(3);
                            loc.xeod(loccounter,1)=obj.posEOD(1);loc.yeod(loccounter,1)=obj.posEOD(2);loc.zeod(loccounter,1)=obj.posEOD(3);
                        
                        case "positionupdater"
                            [posgalvo,posEOD]=component.functionhandle(xest,obj.posgalvo,obj.posEOD,component.parameters{:});
                            obj.posgalvo(component.dim)=posgalvo(component.dim);obj.posEOD(component.dim)=posEOD(component.dim);
                    end
                end
            end
            loc=removeempty(loc,loccounter);
            loc.abortcondition=zeros(size(loc.phot));loc.abortcondition(end)=1+2*bleached;
            out.loc=loc;
            out.raw=photch(1:loccounter,:);
            out.fluorophores=fluorophores;
        end
        function out=runSequence(obj,seq,args)
            arguments
                obj
                seq
                args.maxlocs=1000;
                args.repetitions=1;
            end
            out=[];
            for k=1:args.repetitions
                obj.fluorophores.reset(obj.time);
                out2=obj.runSequenceintern(seq,args.maxlocs);
                out=obj.addout(out,out2);
            end
            out.sequence=seq;
        end
        function out1=addout(obj,out1,out2) %helper function
            if isempty(out2)
                return
            end
            if isempty(out1)
                out1=out2;
                out1.raw=out2.raw;
                out1.loc.rep=1+0*out1.loc.xnm;
                out1.fluorophores.pos=out2.fluorophores.pos;
                out1.fluorophores.int=out2.fluorophores.int;
            else 
                sr=size(out2.raw);
                if length(sr)==2 %Abberior
                    out1.raw(end+1:end+sr(1),1:sr(2))=out2.raw;
                else
                    out1.raw(end+1,:,:,:)=out2.raw;
                end
                out1.fluorophores.pos=cat(1,out1.fluorophores.pos,out2.fluorophores.pos);
                out1.fluorophores.int=cat(1,out1.fluorophores.int,out2.fluorophores.int);
                out2.loc.rep=size(out1.raw,1)+0*out2.loc.xnm;
                out1.loc=appendstruct(out1.loc,out2.loc);
            end
        end

        function pospattern=makeorbitpattern(obj,patterntype, orbitpoints,usecenter,orbitorder,dim)
            switch patterntype
                case "orbitscan"
                    dphi=2*pi/orbitpoints;
                    phi=(0:dphi:2*pi-dphi)';
                    x=cos(phi);
                    y=sin(phi);
                    pospattern=horzcat(x,y,0*phi);
                    if usecenter
                        pospattern(end+1,:)=[0,0,0];
                    end
                    if nargin>3 && ~isempty(orbitorder)
                        pospattern=pospattern(orbitorder,:);
                    end
                case "zscan"
                    pospattern(1,3)=-1;pospattern(2,3)=1;
                    if usecenter
                        pospattern(3,3)=0;
                    end
            end
        end
        function lp=locprec(obj,photons,L)
            lp=L./sqrt(8*(photons));
        end
        function locprec=calculateCRB(obj,patternnames,args)
            arguments
                obj
                patternnames
                args.dim=1:3;
                args.position=obj.fluorophores.position(0);
            end
            dim=args.dim;
            flpos=args.position(1,:);
            % flpos=obj.fluorophores(1).pos; %the main one
            brightness=obj.fluorophores.brightness;
            bg=obj.background/brightness(1);
            
            ih=0;
            % make x,y, z fisher matrix
            eps=1;
            IFisher=zeros(length(dim));
                pattern=obj.patterns(patternnames);
                for k=length(pattern.zeropos):-1:1
                    for coord=1:length(dim)
                        dposa=[0 0 0]; dposa(dim(coord))=eps/2; dposa2=dposa; dposa2(dim(coord))=-eps/2;
                        dpdc(dim(coord))=(pi(dposa)-pi(dposa2))/eps;
                    end
                    for coord=1:length(dim)
                        for coord2=1:length(dim)
                            % IFisher(dim(coord),dim(coord2))=IFisher(dim(coord),dim(coord2))+dpdc(dim(coord))*dpdc(dim(coord2))/(pi([0 0 0])+1e-5);
                             IFisher((coord),(coord2))=IFisher((coord),(coord2))+dpdc(dim(coord))*dpdc(dim(coord2))/(pi([0 0 0])+1e-5);
                        end
                    end
                end
                crlb=(inv(IFisher));
                locprech=diag(sqrt(crlb))';
                locprec=[0 0 0];
                locprec(dim)=locprech;%(dim);

            function iho=pi(dpos)
                    ih=pattern.psf(k).intensity(flpos-pattern.pos(k,:)-obj.posgalvo,dpos,pattern.phasemask(k),pattern.zeropos(k))+bg;
                    
                    ihm=0;
                    for m=length(pattern.zeropos):-1:1
                         ihm=ihm+pattern.psf(m).intensity(flpos-pattern.pos(m,:)-obj.posgalvo,dpos,pattern.phasemask(m),pattern.zeropos(m))+bg;
                    end
                    iho=ih/ihm;
            end
        end
        function stats=displayresults(obj,out)%, lpcrb,L)
          
            photraw=out.raw;
            photraw(photraw==-1)=NaN;
            if length(size(photraw))>3
                photch=squeeze(mean(mean(photraw,1,'omitnan'),2,'omitnan'));
            else
                photch=squeeze(mean(photraw,1,'omitnan'));  
            end
            xest=horzcat(out.loc.xnm,out.loc.ynm,out.loc.znm);
            flpos=horzcat(out.loc.xfl1,out.loc.yfl1,out.loc.zfl1);
            phot=out.loc.phot;
            
            % sigmaCRB=obj.calculateCRB(patternkey,dim=pattern.dim);

            seq=out.sequence;
            lp=[0,0,0]; sigmaCRB=[0,0,0];
            for k=1:length(seq)
                pattern=obj.patterns(seq{k});
                if pattern.type=="pattern"
                    sh=obj.calculateCRB(seq{k},dim=pattern.dim);
                    sigmaCRB(pattern.dim)=sh(pattern.dim);
                    Lh=obj.locprec(mean(phot),pattern.L);
                    lp(pattern.dim)=Lh;
                end
            end

             ff1='%1.1f,';
            disp([ 'photch: ', num2str(photch(:)',ff1),...
                ' mean(phot): ', num2str(mean(phot),ff1)])
             ff='%1.2f,';
            disp(['std:  ', num2str(std(xest,'omitnan'),ff),...
                ' rmse: ', num2str(rmse(xest,flpos,'omitnan'),ff),...
                ' pos: ', num2str(mean(xest,'omitnan'),ff),...
                ' bias: ', num2str(mean(xest-flpos,'omitnan'),ff)]);
            disp(['locp: ', num2str(lp,ff),...
                ' sCRB: ', num2str(sigmaCRB/sqrt(mean(phot)),ff),...
                ' sCRB*sqrt(phot): ', num2str(sigmaCRB,ff1)])

            stats.photch=photch(:);
            stats.phot=mean(phot);
            stats.std=std(xest,'omitnan');
            stats.rmse=rmse(xest,flpos,'omitnan');
            stats.bias=mean(xest-flpos,'omitnan');
            stats.locp=lp;
            stats.sCRB=sigmaCRB/sqrt(mean(phot));
            stats.sCRB1=sigmaCRB;

        end
    end
end


function locs=initlocs(maxlocalizations,fields)
for k=1:length(fields)
    locs.(fields{k})=zeros(maxlocalizations,1);
end
end

function loco=removeempty(loc,loccounter)
fields=fieldnames(loc);
for k=1:length(fields)
    loco.(fields{k})=loc.(fields{k})(1:loccounter);
end
end