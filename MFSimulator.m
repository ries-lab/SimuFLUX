classdef MFSimulator<handle
    properties
        patterns=dictionary;
        sequences=dictionary;
        fluorophores
        background=0; %kHz from AF, does not count towards photon budget
        dwelltime=10; % us
        pospattern=[0 0 0]; %nm, position of pattern center
        time=0;
    end
    methods
        function obj=MFSimulator(fl)
            
            obj.fluorophores=fl;
        end
        function definePattern(obj,patternname,psf,args)
            arguments 
                obj
                patternname
                psf
                args.psfpar="";
                args.zeropos=0;
                args.pos=[0, 0, 0];
                args.makepattern=[];
                args.orbitpoints=4;
                args.orbitL=100;
                args.probecenter=true; 
                args.orbitorder=[];
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
            
            for k=1:length(zeropos)
                pattern.psf(k)=psf;
                pattern.psfpar(k)=string(args.psfpar);
                pattern.zeropos(k)=string((zeropos(k)));
                pattern.psf(k).calculatePSFs(args.psfpar,pattern.zeropos(k))

            end
            obj.patterns(patternname)=pattern;
        end
        
        function combinepatterns(obj,newname, allpatterns)
            patternnew=obj.patterns(allpatterns(1));
            for k=2:length(allpatterns)
                pattern2=obj.patterns(allpatterns(k));
                ln=length(pattern2.zeropos);
                patternnew.psf(end+1:end+ln)=pattern2.psf;
                patternnew.psfpar(end+1:end+ln)=pattern2.psfpar;
                patternnew.zeropos(end+1:end+ln)=pattern2.zeropos;
                patternnew.pos(end+1:end+ln,:)=pattern2.pos;
            end
            obj.patterns(newname)=patternnew;
        end

        function out=patternscan(obj,patternname)
            pattern=obj.patterns(patternname);
            fl=obj.fluorophores;
            numpoints=length(pattern.zeropos);
            intall=zeros(numpoints,1);
            flpos=zeros(1,3);
            timep=0;
            for k=1:numpoints
                inten=0;
                timep=timep+obj.time;
                for f=1:length(fl)
                    flposh=fl(f).position(obj.time);
                    ih=pattern.psf(k).intensity(pattern.psfpar(k),pattern.zeropos(k),flposh-pattern.pos(k,:)-obj.pospattern);
                    inten=inten+fl(f).intensity(ih,obj.dwelltime);
                    obj.time=obj.time+obj.dwelltime;
                    flpos(f,:)=flposh+flpos(f,:);
                end
                intall(k)=inten+obj.background;
            end
            out.phot=poissrnd(intall); %later: fl.tophot(intenall): adds bg, multiplies with brightness, does 
            out.photbg=obj.background*numpoints*obj.dwelltime;
            out.flpos=flpos/numpoints;
            out.time=timep/numpoints;
        end

        function out=patternrepeat(obj,patternname,reps)
            out=obj.patternscan(patternname);
            for rep=2:reps
                out2=obj.patternscan(patternname);
                out.phot=out.phot+out2.phot2;
                out.photbg=out.photbg+out2.photbg2;
                out.flpos=out.flpos+out2.flpos;
                out.time=out.time+out2.time;
            end
            out.flpos=out.flpos/reps;
            out.time=out.time/reps;
        end

        function defineSequence(obj,key,sequencelist)
            obj.sequences(key)={sequencelist};
        end

        function out=runSequence(obj,key,repetitions)
            % obj.fluorophores(1).reset;
            seq=obj.sequences(key);
            seq=seq{1};
            for s=size(seq,1):-1:1
                pat=obj.patterns(seq{s,1});
                ls=length(pat.zeropos);
                photch{s}=zeros(repetitions,ls);
            end
            xest=zeros(repetitions,3);
            flpos=zeros(repetitions,3);
            photall=zeros(repetitions,1);
            time=zeros(repetitions,1);
            for k=1:repetitions
                obj.fluorophores(1).reset;
                numseq=size(seq,1);
                pospattern_beforecenter=obj.pospattern;
                for s=1:numseq
                    outh=obj.patternrepeat(seq{s,1},seq{s,2});
                    xh=seq{s,3}(outh.phot);
                    xest(k,seq{s,4})=xh(seq{s,4});
                    flpos(k,seq{s,4})=outh.flpos(seq{s,4});
                    if seq{s,6} %recenter
                        seq{s,5}(xest(k,:));
                    end
                    photall(k)=photall(k)+sum(outh.phot);
                    photch{s}(k,:)=outh.phot;
                    time(k)=time(k)+outh.time;
                end
                time(k)=time(k)/numseq;
                xest(k,:)=xest(k,:)+pospattern_beforecenter;
            end
            out.photall=photall;
            out.photch=photch;
            out.xest=xest;
            out.flpos=flpos;
            out.time=time;
        end
        

        function xpattern=makeorbitpattern(obj,orbitpoints,usecenter,orbitorder)
            dphi=2*pi/orbitpoints;
            phi=(0:dphi:2*pi-dphi)';%+pi/2;
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
                % crlb=(inv(IFisher(1:2,1:2)));
                crlb=(inv(IFisher));
                locprec=diag(sqrt(crlb))';
                % locprec(1)=sqrt(crlb(1,1)); locprec(2)=sqrt(crlb(2,2));
                % locprec(3)=sqrt(crlb(3,3));

            function iho=pi(dpos)
                    ih=pattern.psf(k).intensity(pattern.psfpar(k),pattern.zeropos(k),flpos-pattern.pos(k,:)-obj.pospattern+dpos)+bg;
                    
                    ihm=0;
                    for m=length(pattern.zeropos):-1:1
                         ihm=ihm+pattern.psf(m).intensity(pattern.psfpar(m),pattern.zeropos(m),flpos-pattern.pos(m,:)-obj.pospattern+dpos)+bg;
                    end
                    iho=ih/ihm;
            end
        end
        function displayresults(obj,out, lp,L)
            % fl=obj.fluorophores(1);
            ff='%1.2f,';
            disp(['mean(phot): ', num2str(mean(out.photall),ff),...
                ' std: ', num2str(std(out.xest,'omitnan'),ff),...
                ' rmse: ', num2str(rmse(out.xest,out.flpos,'omitnan'),ff),...
                ' pos: ', num2str(mean(out.xest,'omitnan'),ff),...
                ' bias: ', num2str(mean(out.xest-out.flpos,'omitnan'),ff),...
                ' locprec: ', num2str(obj.locprec(mean(out.photall),L),ff),...
                ' sqrtCRB: ', num2str(lp,ff)])
        end
    end
end
