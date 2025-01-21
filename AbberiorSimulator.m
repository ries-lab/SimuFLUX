classdef AbberiorSimulator<MFSimulator
    properties
        sequence_json
        estimators
        % sigmapar
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
        function out=runSequence(obj,key)
            obj.fluorophores(1).reset;
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

                out.loc.xgalvo(loccounter,1)=obj.posgalvo(1);out.loc.ygalvo(loccounter,1)=obj.posgalvo(2);out.loc.zgalvo(loccounter,1)=obj.posgalvo(3);

                if abortphot<4 && abortccr<4 %recenter only for valid
                    %estimate position
                    patternpos=obj.patterns(itrname).pos;
                    L=itrs(itr).patGeoFactor*360; %nm
                    estimator=obj.estimators(itr);
                    xest=estimators(estimator.estimator,scanout.phot,patternpos,L,estimator.parameters{:});
                    xestg=xest+obj.posgalvo;
                    
                    %recenter
                    if itr==maxiter
                        dampf=2^(-obj.sequence_json.damping);
                    else
                        dampf=1;
                    end
                    obj.posgalvo=(1-dampf)*obj.posgalvo+dampf*(xestg); 
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
                    cfr=-1;
                end

                out.loc.xnm(loccounter,1)=xestg(1);out.loc.ynm(loccounter,1)=xestg(2);out.loc.znm(loccounter,1)=xestg(3);
                out.loc.xfl(loccounter,1)=scanout.flpos(1)/scanout.counter;out.loc.yfl(loccounter,1)=scanout.flpos(2)/scanout.counter;out.loc.zfl(loccounter,1)=scanout.flpos(3)/scanout.counter;
                out.loc.time(loccounter,1)=scanout.time/scanout.counter;
                out.loc.itr(loccounter,1)=itr;
                out.loc.loccounter(loccounter,1)=loccounter;
                out.loc.cfr(loccounter,1)=cfr;
                out.loc.phot(loccounter,1)=photsum;
                out.loc.vld(loccounter,1)=stickinesscounter<stickiness;
                out.loc.abortcondition(loccounter,1)=1*(abortphot) + 2*(abortccr);
                out.loc.patternrepeat(loccounter,1)=scanout.counter;
                out.loc.measuretime(loccounter,1)=scanout.measuretime;
                % out.loc.efo(loccounter,1)
                out.raw(loccounter,1:length(scanout.phot))=scanout.phot;
                

                loccounter=loccounter+1;
                itr=itr+1;
                if itr>length(itrs)
                    itr=itr+obj.sequence_json.headstart;
                end
                numitr=numitr+1;
            end
        end
    end
end