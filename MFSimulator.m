classdef MFSimulator<handle
    properties
        patterns=dictionary;
        sequences=dictionary;
        fluorophores
        background=0; %kHz from AF, does not count towards photon budget
        % dwelltime=10; % us
        posgalvo=[0 0 0]; %nm, position of pattern center
        posEOD=[0 0 0]; %nm, not descanned, with respect to posgalvo.
        time=0;
    end
    methods
        function obj=MFSimulator(fl)
            if nargin>0    
                obj.fluorophores=fl;
            end
        end
        function definePattern(obj,patternname,psf,args)
            arguments 
                obj
                patternname
                psf
                args.psfpar="";
                args.zeropos=0; %position of the zero when calculating PSFs (e.g., PhaseFLUX)
                args.pos=[0, 0, 0];
                args.makepattern=[];
                args.orbitpoints=4;
                args.orbitL=100;
                args.probecenter=true; 
                args.orbitorder=[];
                args.pointdwelltime=10; %us
                args.laserpower=1; %usually we use relative, but can also be absolute.
            end
            pos=args.pos;
            zeropos=args.zeropos;
            if strcmp(args.makepattern,'orbitscan')
                pos=obj.makeorbitpattern(args.orbitpoints,args.probecenter,args.orbitorder)*args.orbitL/2;
            end
            if size(pos,1)==1 && length(zeropos)>1
                pos=repmat(pos,length(zeropos),1);
            end
            if size(pos,1)>1 && length(zeropos)<=1
                zeropos=repmat(zeropos,1,size(pos,1));
            end
            pattern.pos=pos;
            
            for k=1:size(pos,1)
                pattern.psf(k)=psf;
                pattern.psfpar(k)=string(args.psfpar);
                pattern.zeropos(k)=string((zeropos(k)));
                pattern.psf(k).calculatePSFs(args.psfpar,pattern.zeropos(k));
                pattern.pointdwelltime(k)=args.pointdwelltime;
                pattern.laserpower(k)=args.laserpower;
            end
            obj.patterns(patternname)=pattern;
        end
        
        function combinepatterns(obj,newname, allpatterns)
            patternnew=obj.patterns(allpatterns(1));
            for k=2:length(allpatterns)
                pattern2=obj.patterns(allpatterns(k));
                patternnew=appendstruct(patternnew,pattern2);
            end
            obj.patterns(newname)=patternnew;
        end

        function out=patternscan(obj,patternname)
            pattern=obj.patterns(patternname);
            fl=obj.fluorophores;
            numpoints=length(pattern.zeropos);
            intall=zeros(numpoints,1);
            % flpos=zeros(length(fl),3);
            % flint=zeros(length(fl),1);
            
            timestart=obj.time;
           
            timep=0;
            posgalvo=obj.posgalvo;
            poseod=obj.posEOD;
            flpos=zeros(fl.numberOfFluorophores,3);
            flintall=zeros(fl.numberOfFluorophores,1);
            for k=1:numpoints
                timep=timep+obj.time;
                [flposh,isactive]=fl.position(obj.time);
                flposrel=flposh-posgalvo;
                [ih,phfac]=pattern.psf(k).intensity(flposrel,pattern.pos(k,:)+poseod,pattern.psfpar(k),pattern.zeropos(k));
                ih=ih*pattern.laserpower(k);
                flint=fl.intensity(ih,pattern.pointdwelltime(k),phfac);
                inten=sum(flint);
                flpos(isactive,:)=flpos(isactive,:)+flposh;
                flintall(isactive,:)=flintall(isactive,:)+flint;
                obj.time=obj.time+pattern.pointdwelltime(k);
                intall(k)=inten+obj.background;
            end
            out.phot=poissrnd(intall); %later: fl.tophot(intenall): adds bg, multiplies with brightness, does 
            out.photbg=obj.background*sum(pattern.pointdwelltime);
            out.flpos=flpos/numpoints;
            out.flint=flintall;
            out.time=timep/numpoints;
            out.measuretime=obj.time-timestart;
            out.counter=1;
            fl.updateonoff(obj.time)
        end

        function out=patternrepeat(obj,patternname,reps)
            out=obj.patternscan(patternname);
            for rep=2:reps
                out2=obj.patternscan(patternname);
                out=sumstruct(out,out2);
            end
            out.flpos=out.flpos/out.counter;
            out.time=out.time/out.counter;
            out.counter=1; %already normalized
        end

        function defineSequence(obj,key,sequencelist)
            obj.sequences(key)={sequencelist};
        end

        function out=runSequence(obj,key,maxlocalizations)
            seq=obj.sequences(key);seq=seq{1};
            numseq=size(seq,1);
            for s=numseq:-1:1
                pat=obj.patterns(seq{s,1});
                ls(s)=length(pat.zeropos);
            end
            photch=zeros(maxlocalizations,numseq,max(ls))-1;

            xest=zeros(maxlocalizations,3);
            flpos=zeros(maxlocalizations,3);
            photall=zeros(maxlocalizations,1);
            timep=zeros(maxlocalizations,1);
            bleached=false;
            for k=1:maxlocalizations
                if obj.fluorophores.remainingphotons<1
                    bleached=true;
                    break
                end
                pospattern_beforecenter=obj.posgalvo;
                
                for s=1:numseq
                    scanout=obj.patternrepeat(seq{s,1},seq{s,2});
                    xh=seq{s,3}(scanout.phot);
                    xest(k,seq{s,4})=xh(seq{s,4});
                    flpos(k,seq{s,4})=scanout.flpos(seq{s,4});
                    if seq{s,6} %recenter
                        seq{s,5}(xest(k,:));
                    end
                    photall(k)=photall(k)+sum(scanout.phot);
                    photch(k,s,1:length(scanout.phot))=scanout.phot;
                    timep(k)=timep(k)+scanout.time;
                end
                flpos(k,:)=scanout.flpos;
                timep(k)=timep(k)/numseq;
                xest(k,:)=xest(k,:)+pospattern_beforecenter;
            end
            if bleached
                k=k-1;
            end
            out.raw(1,:,:,:)=photch;
            out.loc.xnm=xest(1:k,1);out.loc.ynm=xest(1:k,2);out.loc.znm=xest(1:k,3);
            out.loc.xfl(1:k,1)=scanout.flpos(1);out.loc.yfl(1:k,1)=scanout.flpos(2);out.loc.zfl(1:k,1)=scanout.flpos(3);
            out.loc.time=timep(1:k,1);
            out.loc.phot=photall(1:k,1);
            out.loc.abortcondition=zeros(size(out.loc.phot));out.loc.abortcondition(end)=1;
        end
        function out=repeatSequence(obj,key,maxloc,repetitions)
            out=[];
            for k=1:repetitions
                obj.fluorophores(1).reset;
                out2=obj.runSequence(key,maxloc);
                out=obj.addout(out,out2);
            end
        end
        function out1=addout(obj,out1,out2)
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

        function xpattern=makeorbitpattern(obj,orbitpoints,usecenter,orbitorder)
            dphi=2*pi/orbitpoints;
            phi=(0:dphi:2*pi-dphi)';
            x=cos(phi);
            y=sin(phi);
            xpattern=horzcat(x,y,0*phi);
            if usecenter
                xpattern(end+1,:)=[0,0,0];
            end
            if nargin>3 && ~isempty(orbitorder)
                xpattern=xpattern(orbitorder,:);
            end
        end
        function lp=locprec(obj,photons,L)
            lp=L./sqrt(8*(photons));
        end

        function locprec=calculateCRB(obj,patternnames,args)
            arguments
                obj
                patternnames
                args.dim=3;
                args.position=obj.fluorophores(1).position(0);
            end
            dim=args.dim;
            flpos=args.position;
            % flpos=obj.fluorophores(1).pos; %the main one
            bg=obj.background/obj.fluorophores(1).brightness;
            
            ih=0;
            % make x,y, z fisher matrix
            eps=1;
            IFisher=zeros(dim);

                pattern=obj.patterns(patternnames);
                for k=length(pattern.zeropos):-1:1
                    for coord=1:dim
                        dposa=[0 0 0]; dposa(coord)=eps/2; dposa2=dposa; dposa2(coord)=-eps/2;
                        dpdc(coord)=(pi(dposa)-pi(dposa2))/eps;
                    end
                    for coord=1:dim
                        for coord2=1:dim
                            IFisher(coord,coord2)=IFisher(coord,coord2)+dpdc(coord)*dpdc(coord2)/(pi([0 0 0])+1e-5);
                        end
                    end
                end
                crlb=(inv(IFisher));
                locprec=diag(sqrt(crlb))';

            function iho=pi(dpos)
                    ih=pattern.psf(k).intensity(flpos-pattern.pos(k,:)-obj.posgalvo,dpos,pattern.psfpar(k),pattern.zeropos(k))+bg;
                    
                    ihm=0;
                    for m=length(pattern.zeropos):-1:1
                         ihm=ihm+pattern.psf(m).intensity(flpos-pattern.pos(m,:)-obj.posgalvo,dpos,pattern.psfpar(m),pattern.zeropos(m))+bg;
                    end
                    iho=ih/ihm;
            end
        end
        function displayresults(obj,out, lp,L)
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
             ff1='%1.1f,';
            disp([ 'photch: ', num2str(photch(:)',ff1),...
                ' mean(phot): ', num2str(mean(phot),ff1)])
             ff='%1.2f,';
            disp(['std: ', num2str(std(xest,'omitnan'),ff),...
                ' rmse: ', num2str(rmse(xest,flpos,'omitnan'),ff),...
                ' pos: ', num2str(mean(xest,'omitnan'),ff),...
                ' bias: ', num2str(mean(xest-flpos,'omitnan'),ff),...
                ' locprec: ', num2str(obj.locprec(mean(phot),L),ff),...
                ' sqrtCRB: ', num2str(lp,ff)])
        end
    end
end




