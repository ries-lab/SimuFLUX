classdef AbberiorSimulator<MFSimulator
    properties
        sequence_json
        estimators
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
                psfscout=PSFMF_donut2D;
                psfmf=PSFMF_donut2D;
            end
            itrs=obj.sequence_json.Itr;
            for k=1:length(itrs)
                switch itrs(k).Mode.id
                    case 'prbt'
                        psf=psfscout;
                        est(k)="donut2D";
                    case 'mflx'
                        psf=psfmf;
                        est(k)="donut2D";
                    otherwise
                        disp([itrs(k).Mode.id]+ " not implemented");
                end
                switch itrs(k).Mode.pattern
                    case 'hexagon'
                        patternpoints=6;
                    otherwise
                        disp([itrs(k).Mode.id]+ " not implemented");                        
                end
                if itrs(k).ccrLimit==-1
                    probecenter=false;
                else
                    probecenter=true;
                end
                L=itrs(k).patGeoFactor*360; %nm
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
            
            stickiness=obj.sequence_json.stickiness;
            loclimit=obj.sequence_json.locLimit;
            if loclimit==-1
                loclimit=10000;%instead of inf. safety to avoid running forever
            end

            numitr=0;abortccr=0;abortphot=0;itr=1;
            while numitr<loclimit && abortphot<stickiness && abortccr<stickiness 
                %scan iteration
                itrname="itr"+itr;
                stickiness=obj.sequence_json.stickiness;
                photsum=0;abortccr=0;abortphot=0;
                while photsum<itrs(itr).phtLimit && abortphot<stickiness && abortccr<stickiness 
                    out=obj.patternrepeat(itrname,itrs(itr).patRepeat);
                    photsum=photsum+sum(out.phot);
                    
                    if itrs(itr).ccrLimit>-1
                        cfr=out.phot(end)/obj.sequence_json.ctrDwellFactor/sum(out.phot(1:end-1));
                        if cfr>itrs(itr).ccrLimit
                            abortccr=abortccr+1;
                        end
                    end
                    
                    if sum(out.phot)<itrs(itr).bgcThreshold*out.measuretime*1e-6
                        abortphot=abortphot+1;
                    end
                end
                
                if abortphot<4 && abortccr<4 %recenter only for valid
                    %estimate position
                    patternpos=obj.patterns(itrname).pos;
                    L=itrs(itr).patGeoFactor*360; %nm
                    fwhm=350; %nm arbitrary
                    xest=estimators(obj.estimators(itr),out.phot,patternpos,L,fwhm);
                    
                    %recenter
                    dampf=2^(-obj.sequence_json.damping);
                    obj.pospattern=(1-dampf)*obj.pospattern+dampf*(xest); 
                end
    
                itr=itr+1;
                if itr>length(itrs)
                    itr=itr+obj.sequence_json.headstart;
                end
                numitr=numitr+1;
            end
        end
    end
end