classdef SimSequencefile<Simulator
    properties
        sequence
        estimators=struct('function',"","par",[],"dim",[]);
        scoutingcoordinates
        psfvec=PsfVectorial; %store only one copy in which different patterns are defined
        % bgcSenseValue=0;
    end
    methods
        function obj=SimSequencefile(varargin)
            obj@Simulator(varargin{:})
            obj.deadtimes=struct('point',0.011,'pattern',0,'estimator',0.015,'positionupdate',0,'localization',0.4);
        end
        function loadsequence(obj,varargin)
            update=true;
            for k=1:length(varargin)
                if strcmp(varargin{k},"noupdate")
                    update=false;
                    continue
                end
                fid=fopen(varargin{k});
                raw=fread(fid,inf);
                str=char(raw');
                fclose(fid);
                obj.sequence=copyfields(obj.sequence,jsondecode(str));
            end
            if update
                obj.makepatterns;
            end
            % obj.sequence_json=jsondecode(str);
        end
        function makepatterns(obj,psfs,phasemasks)
            if nargin<2
                if isfield(obj.sequence,'PSF') && isfield(obj.sequence.PSF,'global')
                    [psfs,phasemasks]=psf_sequence(obj.sequence.PSF,obj.psfvec,obj.sequence);
                else
                    psfs{1}=PsfGauss2D;
                    psfs{2}=PsfDonut2D;
                    phasemasks{1}.phasemask='gauss2D';
                    phasemasks{2}.phasemask='donut2D';
                    obj.loadsequence('defaultestimators.json',"noupdate")
                end
            end
            % if ~isfield(obj.sequence.PSF,'Estimators')
            %     obj.loadsequence('defaultestimators.json')
            % end
            itrs=obj.sequence.Itr;
            for k=1:length(itrs)
                if itrs(k).ccrLimit==-1
                    probecenter=false;
                else
                    probecenter=true;
                end
                L=itrs(k).patGeoFactor*360; %nm
                % betas=[0 0 0 0 0; 1. 4  150 600 -9000; 1 3 20 0 8000; 1 3 60 0 8000];
                % betas=[0 0 0 0 0; 1. 4  0 0 0; 1. 4  0 0 0; 1 3 60 0 8000; ];
                kmin=min(k,length(phasemasks)); % stay for following iterations in last psf
                psf=psfs{kmin};
                phasemask=phasemasks{kmin}.phasemask;

                kmin=min(k,length(phasemasks));
                estimatorh=obj.sequence.PSF.Itr(kmin).estimator;
                % estimatorh.par=replaceinlist(estimatorh.par,"L",L,"probecenter",probecenter); %some default words in parameter to be replaced by values
                obj.estimators(k)=estimatorh;

                switch itrs(k).Mode.pattern
                    case 'hexagon'
                        patternpoints=6;
                    case 'square'
                        patternpoints=4;
                    case 'triangle'
                        patternpoints=3;
                    otherwise
                        disp([itrs(k).Mode.id]+ " not implemented");                        
                end
                
                patterntime=(itrs(k).patDwellTime/itrs(k).patRepeat*(1+probecenter*obj.sequence.ctrDwellFactor))*1e3; %ms
                pointdwelltime=itrs(k).patDwellTime/itrs(k).patRepeat*1e3/patternpoints;
                if probecenter
                    pointdwelltime(2)=pointdwelltime(1)*patternpoints*obj.sequence.ctrDwellFactor;
                end
                % pointdwelltime=patterntime/(patternpoints+probecenter);
                laserpower=itrs(k).pwrFactor;
                obj.definePattern("itr"+k, psf, phasemask=phasemask, makepattern='orbitscan', orbitpoints=patternpoints, orbitL=L,...
                probecenter=probecenter,pointdwelltime=pointdwelltime, laserpower=laserpower,repetitions=itrs(k).patRepeat);
            end
            % obj.estimators=est;
        end
        function out=runSequenceintern(obj)
            %compatible with MFSimulator: repetitions.
            out=[];
            itrs=obj.sequence.Itr;
            maxiter=length(itrs);
            starttime=obj.time;

            deadtimes=obj.deadtimes;
            
            stickiness=obj.sequence.stickiness;
            loclimit=obj.sequence.locLimit;
            if ~isfield(obj.sequence,'maxOffTime') || ~isnumeric(obj.sequence.maxOffTime)&&strcmp(obj.sequence.maxOffTime,'unspecified')
                maxOffTime=3; %ms
            else
                maxOffTime=obj.sequence.maxOffTime*1e6; % to us
            end
            offtimestamp=0;

            if loclimit==-1
                loclimit=10000;%instead of inf. safety to avoid running forever
            end
            loccounter=1;
            
            numitr=0;itr=1;
            abortphot=false;
            stickinesscounter=0;
            xest=[0,0,0];
            par=cell(maxiter,1);

            while numitr<loclimit && stickinesscounter<stickiness && ~abortphot
                %scan iteration
                itrname="itr"+itr;
                stickiness=obj.sequence.stickiness;
                photsum=0;abortccr=0;abortphot=0;
                scanout=[];
                while photsum<itrs(itr).phtLimit && stickinesscounter<stickiness && ~abortphot
                    % repeat scanning until sufficient photons collected,
                    % or abort condition met
                    scanouth=obj.patternscan(itrname);
                    scanouth=backgroundsubtractor(scanouth,obj.background_estimated);
                    scanout=sumstruct(scanout,scanouth);
                    
                    
                    photsum=sum(scanout.phot);
                    photbg=scanout.bg_photons_gt;
                    bgest=scanout.bgphot_est;
                    
                    if itrs(itr).ccrLimit>-1
                        probecenter=true;
                        cfr=scanout.phot(end)/obj.sequence.ctrDwellFactor/sum(scanout.phot(1:end-1));
                        abortccr=cfr>itrs(itr).ccrLimit;
                    else
                        probecenter=false;
                        cfr=-1;
                    end
                    minphot=itrs(itr).bgcThreshold*scanouth.time.patterntotaltime*1e-6;
                    if sum(scanouth.phot)<minphot 
                        if scanouth.time.averagetime>offtimestamp+maxOffTime
                            abortphot=true;
                        end
                    else
                        abortphot=false;
                        offtimestamp=scanouth.time.averagetime;
                    end
                    if abortccr
                        stickinesscounter=stickinesscounter+1;
                    else
                        stickinesscounter=0;
                    end
                end

                if photsum>0 %no photon detected: we don't write localization
                    out.loc.xgalvo(loccounter,1)=obj.posgalvo(1);out.loc.ygalvo(loccounter,1)=obj.posgalvo(2);out.loc.zgalvo(loccounter,1)=obj.posgalvo(3);
                    out.loc.xeod(loccounter,1)=obj.posEOD(1);out.loc.yeod(loccounter,1)=obj.posEOD(2);out.loc.zeod(loccounter,1)=obj.posEOD(3);
                    
                    
                    if ~abortphot && ~abortccr %recenter only for valid
                        vldphotccr=true;
                        %estimate position
                        patternpos=obj.patterns(itrname).pos;
                        L=itrs(itr).patGeoFactor*360; %nm
                        estimator=obj.estimators(itr);
                        estf=str2func(estimator.function);
                        estpar=estimator.par;
                        estpar=replaceinlist(estpar,'patternpos',patternpos,'L',L,'probecenter',probecenter,...
                            'background_est',obj.background_estimated(min(itr,length(obj.background_estimated))),'iteration',itr);
                        xesth=estf(scanout.photrate,estpar{:});
                        xest(estimator.dim)=xesth(estimator.dim);
                        obj.time=obj.time+deadtimes.estimator;

                        % xest=estimators(estimator.estimator,scanout.phot,patternpos,L,estimator.parameters{:});
                        xesttot=xest+obj.posgalvo+obj.posEOD;
                        
                        %recenter
                        if itr==maxiter && ~any(isnan(xesttot))
                            dampf=2^(-obj.sequence.damping);
                            xold=obj.posgalvo;
                            obj.posgalvo=(1-dampf)*obj.posgalvo+dampf*(xesttot);
                            obj.posEOD=obj.posEOD+xold-obj.posgalvo+xest;
                            obj.time=obj.time+deadtimes.positionupdate;
                        elseif ~any(isnan(xesttot))
                            obj.posEOD=obj.posEOD+xest;
                        end
                    else
                        xesttot=[0,0,0]+NaN;
                        vldphotccr=false;
                    end
    
                    if itrs(itr).ccrLimit>-1 %probe center
                        out.loc.eco(loccounter,1)=sum(scanout.phot(1:end-1));
                        out.loc.ecc(loccounter,1)=sum(scanout.phot(end));
                        
                        % orbittime=scanout.patterntotaltime/(1+probecenter*obj.sequence.ctrDwellFactor);
                        % out.loc.efo=out.loc.eco/(orbittime)*1e6;
                        out.loc.efo=out.loc.eco/(sum(scanout.par.pattern.pointdwelltime(1:end-1)))/scanout.repetitions;
                        % out.loc.efc=out.loc.ecc/(scanout.patterntotaltime-orbittime)*1e6;
                        out.loc.efc=out.loc.ecc/(scanout.par.pattern.pointdwelltime(end))/scanout.repetitions;
                    else
                        out.loc.eco(loccounter,1)=sum(scanout.phot);
                        out.loc.efo=out.loc.eco/(sum(scanout.par.pattern.pointdwelltime))/scanout.repetitions;
                        % out.loc.efo=out.loc.eco/(scanout.patterntotaltime);
                        out.loc.ecc(loccounter,1)=-1;
                        out.loc.efc(loccounter,1)=-1;
                        cfr=-1;
                    end
    
                    out.loc.xnm(loccounter,1)=xesttot(1);out.loc.ynm(loccounter,1)=xesttot(2);out.loc.znm(loccounter,1)=xesttot(3);
                    out.loc.xfl1(loccounter,1)=scanout.flpos(1,1)/scanout.counter;out.loc.yfl1(loccounter,1)=scanout.flpos(1,2)/scanout.counter;out.loc.zfl1(loccounter,1)=scanout.flpos(1,3)/scanout.counter;
                    out.loc.time(loccounter,1)=scanout.time.averagetime;
                    out.loc.itr(loccounter,1)=itr;
                    out.loc.numitr(loccounter,1)=numitr;
                    out.loc.loccounter(loccounter,1)=loccounter;
                    out.loc.cfr(loccounter,1)=cfr;
                    out.loc.phot(loccounter,1)=photsum;
                    out.loc.vld(loccounter,1)=stickinesscounter<stickiness & vldphotccr;
                    out.loc.abortcondition(loccounter,1)=1*(abortphot) + 2*(abortccr);
                    out.loc.patternrepeat(loccounter,1)=scanout.counter;
                    out.loc.measuretime(loccounter,1)=scanout.time.patterntotaltime;           
                    out.raw(loccounter,1:length(scanout.phot))=scanout.phot;
                    out.fluorophores.pos(loccounter,1:size(scanout.flpos,1),:)=scanout.flpos/scanout.counter;
                    out.fluorophores.int(loccounter,1:size(scanout.flpos,1))=scanout.flint;
                    out.bg_photons_gt(loccounter,1)=photbg; 
                    out.bg_photons_est(loccounter,1)=sum(bgest);
                    
                    loccounter=loccounter+1;
                    if isempty(par{itr})
                        par{itr}=scanout.par;
                    end
                end

                itr=itr+1;
                if itr>length(itrs)
                    itr=itr+obj.sequence.headstart;
                end
                numitr=numitr+1;
            end
            if ~isempty(out)
                out.sequence={"itr"+max(out.loc.itr)};
                out.duration=obj.time-starttime;
                out.par=par;
            end
            obj.time=obj.time+deadtimes.localization;
            
        end
        function out=runSequence(obj,args)
            arguments
                obj
                % args.maxlocs=10;
                args.repetitions=1;
                args.resetfluorophores=false;
            end
            starttime=obj.time;
            out=[];
            for k=1:args.repetitions
                if args.resetfluorophores
                    obj.fluorophores(1).reset;
                end
                out2=obj.runSequenceintern;
                out=obj.appendout(out,out2);
            end
            if ~isempty(out)
                out.duration=obj.time-starttime;
            end
        end
        function makescoutingpattern(obj,fov,args)
            %fov=[x,y;x2,y2]
            arguments
                obj
                fov
                args.distance=0;
                args.show=false;
            end
            if args.distance==0       
                geofactor=obj.sequence.field.fldGeoFactor;
                args.distance=geofactor*360/640*obj.sequence.Itr(1).wavelength*1e9;
            end
            obj.scoutingcoordinates=makehexgrid(fov,args.distance);
            if args.show
                figure(98);plot(obj.scoutingcoordinates(:,1),obj.scoutingcoordinates(:,2),'o')
            end
        end
        function out=scoutingSequence(obj,args)
            arguments
                obj
                args.maxrep=100000;
                args.maxtime=1e6;
            end
            timestart=obj.time;
            %reset fluorophore?
            % obj.posgalvo(1:2)=obj.scoutingcoordinates(1,:);
            % obj.posEOD=[0 0 0];
            out=[];
            allbleached=false;
            % fprintf(1,'Computation Progress: %3d%%\n',0);
            fprintf(1,"scouting, progress: %3d%%\n",0)
            for reps=1:args.maxrep
                if obj.time>timestart+args.maxtime
                    break
                end
                prog=reps/args.maxrep*100;
                prog=max(prog, (obj.time-timestart)/args.maxtime*100);
                fprintf(1,'\b\b\b\b%3.0f%%',prog)
                % fprintf(num2str(reps,"%2.0f,"))
                for pind=1:size(obj.scoutingcoordinates,1)
                    obj.posgalvo(1:2)=obj.scoutingcoordinates(pind,:);
                    obj.posEOD=[0 0 0];
                    out2=obj.runSequence;
                    out=obj.appendout(out,out2);
                    if obj.fluorophores.allbleached
                        allbleached=true;
                        break
                    end
                end
                if allbleached
                    break
                end
            end
            if ~isempty(out)
                out.duration=obj.time-timestart;
            end
        end
        % function displayresults(obj)
        %     keys=obj.patterns.keys;
        %     displayresults@Sim_Simulator(obj,keys(1),out)
        % end
        function plotpositions(obj, out, args)
            arguments 
                obj
                out
                args.figure=[];
                args.coordinate=1;
                args.axis=[];
                args.xvalues="itr";
            end
            if isempty(out)
                disp('no localizations found')
            end
            xnmn={"xnm","ynm","znm"};
            xfln={"xfl1","yfl1","zfl1"};
            xgn={"xgalvo","ygalvo","zgalvo"};
            xen={"xeod","yeod","zeod"};
            
            if isempty(args.axis)
                if ~isempty(args.figure)
                    f=figure(args.figure);
                end
                ax=gca;
            else
                ax=args.axis;
            end
            switch args.xvalues
                case "itr"
                    xv=out.loc.loccounter;
                    xtxt="time (itr)";
                case "time"
                    xv=out.loc.time;
                    xtxt="time (ms)";
            end
            c=args.coordinate;
            hold(ax,'off')
            plot(ax,xv, out.loc.(xnmn{c}),'k');
            hold(ax,"on")
            plot(ax,xv,out.loc.(xfln{c}),'r');
            
            plot(ax,xv,out.loc.(xgn{c}),"Color",[0 0.7 0])
            plot(ax,xv,out.loc.(xen{c}),'b')
            
            xlabel(ax,xtxt)
            ylabel(ax,'x position(nm)')
            legend(ax,'Estimated','Fluorophore','Galvo','EOD')

        end
        % function so=subtractbackground(obj,si)
        %     bg=obj.bgcSenseValue;
        %     so=si;
        %     if bg>0
        %         pointtime=si.par.pointdwelltime;
        %         so.phot=so.phot-bg*pointtime;
        %         % so.photbg=so.photbg+bg*pointtime; % this is explicitely used in some
        %         % estimators... not consistent
        %     end
        % end
    end
end


function pos=makehexgrid(roi,d)
h=d*cosd(30);
numpx=(roi(2,1)-roi(1,1))/h;
numpy=(roi(2,2)-roi(1,2))/d;
pos=zeros(ceil(numpx*numpy),2);
ind=1;
for l=numpy:-1:-numpx
    for k=0:numpx
        x=k*h+roi(1,1);
        y=l*d+k*d/2+roi(1,2);
        if x<=roi(2,1) && y<=roi(2,2) && x>=roi(1,1) &&y>=roi(1,2)
            pos(ind,:)=[x,y];
            ind=ind+1;
        end
    end
end
pos(ind:end,:)=[];
end
