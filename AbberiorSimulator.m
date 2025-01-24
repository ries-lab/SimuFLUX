classdef AbberiorSimulator<MFSimulator
    properties
        sequence_json
        estimators
        scoutingcoordinates
    end
    methods
        function loadsequence(obj,fname)
            fid=fopen(fname);
            raw=fread(fid,inf);
            str=char(raw');
            fclose(fid);
            obj.sequence_json=jsondecode(str);
        end
        function makepatterns(obj,psfs)
            if nargin<2
                psfscout=PSFMF_gauss2D;
                psfmf=PSFMF_donut2D;
            end
            itrs=obj.sequence_json.Itr;
            for k=1:length(itrs)
                if itrs(k).ccrLimit==-1
                    probecenter=false;
                else
                    probecenter=true;
                end
                L=itrs(k).patGeoFactor*360; %nm
                % L
                betas=[0 0 0 0 0; 1. 4  150 600 -9000; 1 3 20 0 8000; 1 3 60 0 8000];
                % betas=[0 0 0 0 0; 1. 4  0 0 0; 1. 4  0 0 0; 1 3 60 0 8000; ];

                switch itrs(k).Mode.id
                    case 'prbt'
                        psf=psfscout;
                        est(k).estimator="gauss2D";
                        est(k).parameters={psf.sigma, probecenter};
                        % obj.sigmapar(k)=psf.sigma;
                    case 'mflx'
                        psf=psfmf;
                        if probecenter
                            est(k).estimator="donut2DcentermLMS";
                            est(k).parameters={psf.fwhm,betas(k,:)};
                        else
                            est(k).estimator="donut2D";
                            est(k).parameters={psf.fwhm};
                        end
                        % est(k).parameters={psf.fwhm,[1 5]};
                        % obj.sigmapar(k)=psf.fwhm;
                    otherwise
                        disp([itrs(k).Mode.id]+ " not implemented");
                end
                switch itrs(k).Mode.pattern
                    case 'hexagon'
                        patternpoints=6;
                    otherwise
                        disp([itrs(k).Mode.id]+ " not implemented");                        
                end
                
                patterntime=(itrs(k).patDwellTime/itrs(k).patRepeat*(1+probecenter*obj.sequence_json.ctrDwellFactor))*1e6; %us
                pointdwelltime=patterntime/(patternpoints+probecenter);
                laserpower=itrs(k).pwrFactor;
                obj.definePattern("itr"+k, psf, makepattern='orbitscan', orbitpoints=patternpoints, orbitL=L,...
                probecenter=probecenter,pointdwelltime=pointdwelltime, laserpower=laserpower);
            end
            obj.estimators=est;
        end
        function out=runSequence(obj,key,~)
            % obj.fluorophores(1).reset;
            out=[];
            itrs=obj.sequence_json.Itr;
            maxiter=length(itrs);
            
            stickiness=obj.sequence_json.stickiness;
            loclimit=obj.sequence_json.locLimit;
            maxOffTime=obj.sequence_json.maxOfftime;
            if ~isnumeric(maxOffTime)&&strcmp(maxOffTime,'unspecified')
                maxOffTime=3000; %us
            else
                maxOffTime=maxOffTime*1e6; % to us
            end
            offtimestamp=0;

            if loclimit==-1
                loclimit=10000;%instead of inf. safety to avoid running forever
            end
            loccounter=1;
            
            numitr=0;itr=1;
            stickinesscounter=0;

            while numitr<loclimit && stickinesscounter<stickiness 
                %scan iteration
                itrname="itr"+itr;
                stickiness=obj.sequence_json.stickiness;
                photsum=0;abortccr=0;abortphot=0;
                scanout=[];
                while photsum<itrs(itr).phtLimit && stickinesscounter<stickiness 
                    scanouth=obj.patternrepeat(itrname,itrs(itr).patRepeat);
                    scanout=sumstruct(scanout,scanouth);
                    
                    photsum=sum(scanout.phot);
                    
                    if itrs(itr).ccrLimit>-1
                        probecenter=true;
                        cfr=scanout.phot(end)/obj.sequence_json.ctrDwellFactor/sum(scanout.phot(1:end-1));
                        abortccr=cfr>itrs(itr).ccrLimit;
                    else
                        probecenter=false;
                        cfr=-1;
                    end
                    minphot=itrs(itr).bgcThreshold*scanouth.measuretime*1e-6;
                    if sum(scanouth.phot)<minphot 
                        if scanouth.time>offtimestamp+maxOffTime
                            abortphot=true;
                        end
                    else
                        abortphot=false;
                        offtimestamp=scanouth.time;
                    end
                    if abortphot || abortccr
                        stickinesscounter=stickinesscounter+1;
                    else
                        stickinesscounter=0;
                    end
                end

                if photsum>0 %no photon detected: we don't write localization
                    out.loc.xgalvo(loccounter,1)=obj.posgalvo(1);out.loc.ygalvo(loccounter,1)=obj.posgalvo(2);out.loc.zgalvo(loccounter,1)=obj.posgalvo(3);
                    out.loc.xeod(loccounter,1)=obj.posEOD(1);out.loc.yeod(loccounter,1)=obj.posEOD(2);out.loc.zeod(loccounter,1)=obj.posEOD(3);
                    
    
                    if ~abortphot && ~abortccr %recenter only for valid
                        %estimate position
                        patternpos=obj.patterns(itrname).pos;
                        L=itrs(itr).patGeoFactor*360; %nm
                        estimator=obj.estimators(itr);
                        xest=estimators(estimator.estimator,scanout.phot,patternpos,L,estimator.parameters{:});
                        xesttot=xest+obj.posgalvo+obj.posEOD;
                        
                        %recenter
                        if itr==maxiter
                            dampf=2^(-obj.sequence_json.damping);
                            xold=obj.posgalvo;
                           
                            obj.posgalvo=(1-dampf)*obj.posgalvo+dampf*(xesttot);
                            obj.posEOD=obj.posEOD+xold-obj.posgalvo+xest;
                        else
                            obj.posEOD=obj.posEOD+xest;
                        end
                    else
                        xesttot=[0,0,0];
                    end
    
                    
                    
                    if itrs(itr).ccrLimit>-1
                        out.loc.eco(loccounter,1)=sum(scanout.phot(1:end-1));
                        out.loc.ecc(loccounter,1)=sum(scanout.phot(end));
                        
                        orbittime=scanout.measuretime/(1+probecenter*obj.sequence_json.ctrDwellFactor);
                        out.loc.efo=out.loc.eco/(orbittime)*1e6;
                        out.loc.efc=out.loc.ecc/(scanout.measuretime-orbittime)*1e6;
                    else
                        out.loc.eco(loccounter,1)=sum(scanout.phot);
                        out.loc.efo=out.loc.eco/(scanout.measuretime);
                        out.loc.ecc(loccounter,1)=-1;
                        out.loc.efc(loccounter,1)=-1;
                        cfr=-1;
                    end
    
                    out.loc.xnm(loccounter,1)=xesttot(1);out.loc.ynm(loccounter,1)=xesttot(2);out.loc.znm(loccounter,1)=xesttot(3);
                    % out.loc.xfl(loccounter,1)=scanout.flpos(1)/scanout.counter;out.loc.yfl(loccounter,1)=scanout.flpos(2)/scanout.counter;out.loc.zfl(loccounter,1)=scanout.flpos(3)/scanout.counter;
                    out.loc.time(loccounter,1)=scanout.time/scanout.counter;
                    out.loc.itr(loccounter,1)=itr;
                    out.loc.numitr(loccounter,1)=numitr;
                    out.loc.loccounter(loccounter,1)=loccounter;
                    out.loc.cfr(loccounter,1)=cfr;
                    out.loc.phot(loccounter,1)=photsum;
                    out.loc.vld(loccounter,1)=stickinesscounter<stickiness;
                    out.loc.abortcondition(loccounter,1)=1*(abortphot) + 2*(abortccr);
                    out.loc.patternrepeat(loccounter,1)=scanout.counter;
                    out.loc.measuretime(loccounter,1)=scanout.measuretime;           
                    out.raw(loccounter,1:length(scanout.phot))=scanout.phot;
                    out.fluorophores.pos(loccounter,1:size(scanout.flpos,1),:)=scanout.flpos;
                    out.fluorophores.int(loccounter,1:size(scanout.flpos,1))=scanout.flint;
                    
                    loccounter=loccounter+1;
                end

                itr=itr+1;
                if itr>length(itrs)
                    itr=itr+obj.sequence_json.headstart;
                end
                numitr=numitr+1;
            end
        end
        function makescoutingpattern(obj,fov,args)
            %fov=[x,y;x2,y2]
            arguments
                obj
                fov
                args.distance=0;
            end
            if args.distance==0       
                geofactor=obj.sequence_json.field.fldGeoFactor;
                args.distance=geofactor*360/640*obj.sequence_json.Itr(1).wavelength*1e9;
            end
            obj.scoutingcoordinates=makehexgrid(fov,args.distance);
            figure(98);plot(obj.scoutingcoordinates(:,1),obj.scoutingcoordinates(:,2),'o')

        end
        function out=scoutingSequence(obj,args)
            arguments
                obj
                args.maxrep=2;
            end
            %reset fluorophore?
            obj.posgalvo(1:2)=obj.scoutingcoordinates(1,:);
            obj.posEOD=[0 0 0];
            out=[];
            allbleached=false;
            for reps=1:args.maxrep
                for pind=2:size(obj.scoutingcoordinates,1)
                    obj.posgalvo(1:2)=obj.scoutingcoordinates(pind,:);
                    obj.posEOD=[0 0 0];
                    out2=obj.runSequence;
                    out=obj.addout(out,out2);
                    if obj.fluorophores.allbleached
                        allbleached=true;
                        break
                    end
                end
                if allbleached
                    break
                end
            end


        end
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