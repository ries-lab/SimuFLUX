classdef Simulator<handle
    properties
        patterns=dictionary;
        sequences=dictionary;
        fluorophores
        background=0; %kHz from AF, does not count towards photon budget. Should scale with PSF.normfact and laser intensity?
        posgalvo=[0 0 0]; %nm, position of pattern center
        posEOD=[0 0 0]; %nm, not descanned, with respect to posgalvo.
        time=0;
        background_estimated=0;
    end
    methods
        function obj=Simulator(args)
            arguments
                args.fluorophores=[];
                args.background=0;
                args.background_estimated=0;
                args.loadfile={};
            end
            obj.fluorophores=args.fluorophores;
            obj.background=args.background;
            obj.background_estimated=args.background_estimated;
            if ~isempty(args.loadfile)
                obj.loadsequence(args.loadfile{:});
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
                args.pointdwelltime=.01; %us
                args.laserpower=1; %usually we use relative, but can also be absolute.
                args.repetitions=1;
                args.dim=1:2;
            end
            pattern.repetitions=args.repetitions;
            pattern.par=args;
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
                pattern.backgroundfac(k)=1/psf.normfactor(args.phasemask,pattern.zeropos(k));
               
                pattern.laserpower(k)=args.laserpower;
            end
            if length(args.pointdwelltime)==size(pos,1) %one dwelltime per position
                pattern.pointdwelltime=args.pointdwelltime;
            else
                pattern.pointdwelltime=zeros(1,size(pos,1))+args.pointdwelltime(1);
            end
            if length(args.pointdwelltime)==2
                pattern.pointdwelltime(end)=args.pointdwelltime(2);
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
            bgphot=0;
            time=obj.time;
            background=obj.background;
            for r=1:repetitions
                for k=1:numpoints
                    timep=timep+time; %for calculating average time point
                    [flposh,isactive]=fluorophores.position(obj.time,flproperties);
                    flposrel=flposh-posgalvo;
                    [intensityh,pinholehfac]=pattern.psf(k).intensity(flposrel(isactive,:),pattern.pos(k,:)+posEOD,pattern.phasemask(k),pattern.zeropos(k));
                    intensityh=intensityh*pattern.laserpower(k);
                    flint=fluorophores.intensity(intensityh,pattern.pointdwelltime(k),time,pinholehfac,flproperties);
                    intensity=sum(flint);
                    flpos(:,:)=flpos(:,:)+flposh;
                    flintall(isactive,:)=flintall(isactive,:)+flint;
                    time=time+pattern.pointdwelltime(k);
                    bgphoth=pattern.backgroundfac(k)*background*pattern.pointdwelltime(k);
                    bgphot=bgphoth+bgphot;
                    intall(k)=intall(k)+intensity+bgphoth; %sum over repetitions, fluorophores
                end
                fluorophores.updateonoff(time);
            end
            out.phot=poissrnd(intall); %later: fl.tophot(intenall): adds bg, multiplies with brightness, does 
            out.photrate=out.phot./pattern.pointdwelltime';out.photrate=out.photrate/sum(out.photrate)*sum(out.phot); 
            out.pointdwelltimes=pattern.pointdwelltime';
            out.bg_photons_gt=bgphot;
            out.bgphot_est=0;
            out.intensity=intall;
            out.flpos=flpos/numpoints/repetitions;
            out.flint=flintall;
            out.averagetime=timep/numpoints/repetitions;
            out.patterntotaltime=time-timestart;
            out.patterntime=sum(pattern.pointdwelltime);
            out.counter=1;
            out.repetitions=repetitions;
            out.par=pattern.par;
            out.par.L=pattern.L;
            out.par.patternpos=pattern.pos;
            out.par.zeropos=pattern.zeropos;
            out.par.dim=pattern.dim;
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
            photbg=zeros(maxlocalizations*numpat,1);
            loc=initlocs(maxlocalizations*numpat,{'xnm','ynm','znm','xfl1',...
                'yfl1','zfl1','phot','time','abortcondition', 'xgalvo', 'ygalvo', 'zgalvo',...
                'xeod','yeod','zeod'});
   
            bleached=false;
            loccounter=0;
            fl.pos=[];fl.int=[];
            for k=1:maxlocalizations
                if obj.fluorophores.remainingphotons<1
                    bleached=true;
                    break
                end
                posgalvo_beforecenter=obj.posgalvo;
                xest=[0,0,0];
                loccounter_seq=loccounter+1; %XXXX 9.2.25
                for s=1:numseq
                    component=obj.patterns(seq{s});
                    switch component.type
                        case "pattern"
                            scanout=obj.patternscan(seq{s});
                            loccounter=loccounter+1; % every localization gets a new entry, as in abberior
                            loc.phot(loccounter)=loc.phot(loccounter)+sum(scanout.phot);
                            photch(loccounter,1:length(scanout.phot))=scanout.phot;
                            photbg(loccounter)=scanout.bg_photons_gt;
                            loc.time(loccounter)=loc.time(loccounter)+scanout.averagetime;
                            flpos=scanout.flpos(1,:);
                            loc.xfl1(loccounter)=flpos(1);loc.yfl1(loccounter)=flpos(2);loc.zfl1(loccounter)=flpos(3);
                            fl.pos(loccounter,1:size(scanout.flpos,1),:)=scanout.flpos;
                            fl.int(loccounter,1:size(scanout.flpos,1))=scanout.flint;
                        case "estimator"
                            % replace placeholder names by values
                            component_par=replaceinlist(component.parameters,'patternpos',scanout.par.patternpos,'L',scanout.par.L,...
                                'probecenter',scanout.par.probecenter,'bg_photons_gt',scanout.bg_photons_gt, ...
                                'background_estimated',obj.background_estimated,'iteration',s);
                            xesth=component.functionhandle(scanout.photrate,component_par{:});
                            if length(xesth)==3
                                xest(component.dim)=xesth(component.dim);
                            else
                                xest(component.dim)=xesth;
                            end
                            xesttot=xest+posgalvo_beforecenter;
                            
                            loc.xgalvo(loccounter,1)=obj.posgalvo(1);loc.ygalvo(loccounter,1)=obj.posgalvo(2);loc.zgalvo(loccounter,1)=obj.posgalvo(3);
                            loc.xeod(loccounter,1)=obj.posEOD(1);loc.yeod(loccounter,1)=obj.posEOD(2);loc.zeod(loccounter,1)=obj.posEOD(3);
                        
                        case "positionupdater"
                            [posgalvo,posEOD]=component.functionhandle(xest,obj.posgalvo,obj.posEOD,component.parameters{:});
                            obj.posgalvo(component.dim)=posgalvo(component.dim);obj.posEOD(component.dim)=posEOD(component.dim);
                        case "background"
                            component_par=replaceinlist(component.parameters,'patternpos',scanout.par.patternpos,'L',scanout.par.L,...
                                'probecenter',scanout.par.probecenter,'bg_photons_gt',scanout.bg_photons_gt, ...
                                'background_estimated',obj.background_estimated,'iteration',s);
                            scanout=component.functionhandle(scanout,component_par{:});

                        otherwise
                            disp(component.type+ " not defined")
                    end
                end
                % write position estimates for all localizations in sequence together, 9.2.25
                loc.xnm(loccounter_seq:loccounter)=xesttot(1);loc.ynm(loccounter_seq:loccounter)=xesttot(2);loc.znm(loccounter_seq:loccounter)=xesttot(3);

            end
            loc=removeempty(loc,loccounter);
            loc.abortcondition=zeros(size(loc.phot));loc.abortcondition(end)=1+2*bleached;
            out.loc=loc;
            out.raw=photch(1:loccounter,:);
            out.fluorophores=fl;
            out.bg_photons_gt=photbg(1:loccounter);
            % out.bg_photons_est=
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
                out=obj.appendout(out,out2);
            end
            out.sequence=seq;
        end
        function out1=appendout(obj,out1,out2) %helper function
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
        function locprec=calculateCRBdirect(obj,patternnames,args)
            arguments
                obj
                patternnames
                args.dim=1:3;
                args.position=obj.fluorophores.position(0);
            end
            dim=args.dim;
            % flpos=args.position(1,:);
            flpos=args.position;
            % flpos=obj.fluorophores(1).pos; %the main one
            brightness=obj.fluorophores.brightness;
            bg=obj.background/brightness(1); %hack, to not have to calculate with all fluorophores

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
                             % pi([0 0 0])
                        end
                    end
                end
                crlb=(inv(IFisher));
                locprech=diag(sqrt(crlb))';
                locprec=[0 0 0];
                locprec(dim)=locprech;%(dim);

            function iho=pi(dpos)
                    flposh=flpos;
                    flposh(1,:)=flposh(1,:)-dpos;
                    bgh=bg*pattern.backgroundfac(k);
                    ih=sum(pattern.psf(k).intensity(flposh-obj.posgalvo,pattern.pos(k,:),pattern.phasemask(k),pattern.zeropos(k))+bgh);

                    ihm=0;
                    for m=length(pattern.zeropos):-1:1
                        bgh=bg*pattern.backgroundfac(m);
                         ihm=ihm+sum(pattern.psf(m).intensity(flposh-obj.posgalvo,pattern.pos(m,:),pattern.phasemask(m),pattern.zeropos(m))+bgh);
                    end
                    iho=ih/ihm;
            end
        end


        function locprec=calculateCRB(obj,patternnames,args)
            arguments
                obj
                patternnames
                args.dim=1:3;
                args.position=obj.fluorophores.position(0);
            end
            dim=args.dim;
            %patternscan with blinking fluorophores cannot be used for CRB
            if isa(obj.fluorophores,'FlCollectionBlinking') || isa(obj.fluorophores,'FlBlinkBleach')|| isa(obj.fluorophores,'FlMoving')
                locprec=calculateCRBdirect(obj,patternnames,dim=args.dim,position=args.position);
                return
            end
            if isa(obj.fluorophores,'FlCollection')
                fl1=obj.fluorophores.flall{1};
            else
                fl1=obj.fluorophores;
            end
            warning('off','MATLAB:structOnObject')
            oldpar=struct(fl1);
            fl1.remainingphotons=inf;
            posold=fl1.position(obj.time);
            eps=1;
            IFisher=zeros(length(dim));
            out0=obj.patternscan(patternnames);
            pi0=out0.intensity/sum(out0.intensity);
            for coord=1:length(dim)
                dposa=[0 0 0]; dposa(dim(coord))=eps/2; dposa2=dposa; dposa2(dim(coord))=-eps/2;

                fl1.pos=posold+dposa;
                out1=obj.patternscan(patternnames);
                fl1.pos=posold+dposa2;  
                out2=obj.patternscan(patternnames);
                dpdc(:,dim(coord))=((out1.intensity/sum(out1.intensity)-out2.intensity/sum(out2.intensity))/eps);
            end
            for coord=1:length(dim)
                for coord2=1:length(dim)
                     fh=dpdc(:,dim(coord)).*dpdc(:,dim(coord2))./(pi0+1e-12);
                     IFisher((coord),(coord2))=sum(fh);
                end
            end

            fl1.pos=oldpar.pos;
            fl1.remainingphotons=oldpar.remainingphotons;
            crlb=(inv(IFisher));
            locprech=diag(sqrt(crlb))';
            locprec=[0 0 0];
            locprec(dim)=locprech;%(dim);
        end
        
        function st=summarize_results(obj,out,args)%, lpcrb,L)
            arguments
                obj
                out
                args.display=true;
                args.filter=true(size(out.loc.xnm));
            end
            ind=args.filter;
            photraw=out.raw;
            photraw(photraw==-1)=NaN;
            if length(size(photraw))>3
                photraw(photraw<0)=NaN;
                photch=squeeze(mean(mean(photraw(ind,:),1,'omitnan'),2,'omitnan'));
            else
                photch=squeeze(mean(photraw(ind,:),1,'omitnan'));  
            end
            xest=horzcat(out.loc.xnm(ind),out.loc.ynm(ind),out.loc.znm(ind));
            flpos=horzcat(out.loc.xfl1(ind),out.loc.yfl1(ind),out.loc.zfl1(ind));
            phot=out.loc.phot(ind);
            
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

            st.photch=photch;
            st.bg_photons=mean(out.bg_photons_gt(ind));
            st.phot=mean(phot,'omitnan');
            st.phot_signal= st.phot-st.bg_photons;    
            st.pos=mean(xest,1,'omitnan');
            st.std=std(xest,[],1,'omitnan');
            st.rmse=rmse(xest,flpos,1,'omitnan');
            
            st.bias=mean(xest-flpos,1,'omitnan');
            st.locp=lp;
            st.sCRB=sigmaCRB/sqrt(st.phot);
            st.sCRB1=sigmaCRB;
            
            if args.display
                ff1='%1.1f,';
                disp([ 'photch: ', num2str(st.photch,ff1),...
                    ' mean(phot): ', num2str(st.phot,ff1), 'signal phot: ', num2str(st.phot_signal,ff1)])
                 ff='%1.2f,';
                disp(['std:  ', num2str(st.std,ff),...
                    ' rmse: ', num2str(st.rmse,ff),...
                    ' pos: ', num2str(st.pos,ff),...
                    ' bias: ', num2str(st.bias,ff)]);
                disp(['locp: ', num2str(lp,ff),...
                    ' sCRB: ', num2str(st.sCRB,ff),...
                    ' sCRB*sqrt(phot): ', num2str(st.sCRB1,ff1)])
            end
        end
        function so=scan_fov(obj,seq,xcoords,args)
            arguments
                obj
                seq
                xcoords
                args.maxlocs=1000;
                args.repetitions=1;
                args.display=true;
                args.dimplot=1;
                args.dimscan=1;
                args.title="FoV scan";
                args.fluorophorenumber=1;
                args.ax1="std";
                % args.ax2="bias";
                args.clearfigure=true;
                args.tag="";
            end
            args.ax1=string(args.ax1);
            % args.ax2=string(args.ax2);

            so.std=zeros(length(xcoords),3); so.rmse=so.std; so.bias=so.std;
            if isa(obj.fluorophores,'FlCollection')
                fl=obj.fluorophores.flall{args.fluorophorenumber};
            else
                fl=obj.fluorophores;
            end
            coords=["x","y","z"];

            for k=1:length(xcoords)
                fl.pos(args.dimscan)=xcoords(k);
                out=obj.runSequence(seq,maxlocs=args.maxlocs,repetitions=args.repetitions);
                stats=obj.summarize_results(out,"display",false);
                so.std(k,:)=stats.std;
                so.bias(k,:)=stats.bias;
                so.rmse(k,:)=stats.rmse;
                so.sCRB(k,:)=stats.sCRB;
                so.pos(k,:)=stats.pos;
            end
            so.stdrel=so.std./so.sCRB;
            so.biasrel=so.bias./so.pos;
            if args.display
                axh=gca;
                if ~isempty(axh.Legend) && ~args.clearfigure
                    ltxt=string(axh.Legend.String);
                else
                    ltxt={};
                end
                % yyaxis left
                if args.clearfigure
                    hold off
                else
                    hold on
                end
                ylab=coords(args.dimplot)+": ";
                for m=1:length(args.ax1)
                    yval=so.(args.ax1(m));
                    plot(xcoords,yval(:,args.dimplot))
                    hold on
                    ylab=ylab+args.ax1(m);
                    if m<length(args.ax1)
                        ylab=ylab +"/";
                    end
                end
                ylabel(ylab + " (nm)")
                % hold off
    
                % yyaxis right
                % if args.clearfigure
                %     hold off
                % else
                %     hold on
                % end
                % ylab2=coords(args.dimplot)+": ";
                % for m=1:length(args.ax2)
                %     yval=so.(args.ax2(m));
                %     plot(xcoords,yval(:,args.dimplot))
                %     hold on
                %     ylab2=ylab2+args.ax2(m);
                %     if m<length(args.ax2)
                %         ylab2=ylab2 +"/";
                %     end
                % end
                % ylabel(ylab2 + " (nm)")
                % % hold off
                % 
                xlabel(coords(args.dimscan) + " position (nm)")
                title(args.title)
                % legend
                % axh=gca;
                % if ~isempty(axh.Legend)
                %     ltxt=axh.Legend.String;
                % else
                %     ltxt={};
                % end
                ltxt=horzcat(ltxt,args.ax1+": "+args.tag);
                % ltxt=horzcat(ltxt,{args.ax1+": "+args.tag args.ax2+": "+args.tag});
                legend(ltxt)
                if any(contains(args.ax1,"bias"))
                    grid on
                    % plot(xcoords,0*xcoords,'r-')
                end
            end
        end
        function loadsequence(obj,varargin)
            warning('no sequence loader implemented for Simulator class')
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