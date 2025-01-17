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
            loccounter=1;
            numitr=0;abortccr=0;abortphot=0;itr=1;
            while numitr<loclimit && abortphot<stickiness && abortccr<stickiness 
                %scan iteration
                itrname="itr"+itr;
                stickiness=obj.sequence_json.stickiness;
                photsum=0;abortccr=0;abortphot=0;
                scanout=[];
                while photsum<itrs(itr).phtLimit && abortphot<stickiness && abortccr<stickiness 
                    scanouth=obj.patternrepeat(itrname,itrs(itr).patRepeat);
                    scanout=sumstruct(scanout,scanouth);
                    
                    photsum=sum(scanout.phot);
                    
                    if itrs(itr).ccrLimit>-1
                        cfr=scanout.phot(end)/obj.sequence_json.ctrDwellFactor/sum(scanout.phot(1:end-1));
                        if cfr>itrs(itr).ccrLimit
                            abortccr=abortccr+1;
                        end
                    else
                        cfr=-1;
                    end
                    minphot=itrs(itr).bgcThreshold*scanouth.measuretime*1e-6;
                    if sum(scanouth.phot)<minphot                    
                        abortphot=abortphot+1;
                    end
                end

                out.loc.xpat(loccounter,1)=obj.pospattern(1);out.loc.ypat(loccounter,1)=obj.pospattern(2);out.loc.zpat(loccounter,1)=obj.pospattern(3);

                if abortphot<4 && abortccr<4 %recenter only for valid
                    %estimate position
                    patternpos=obj.patterns(itrname).pos;
                    L=itrs(itr).patGeoFactor*360; %nm
                    fwhm=350; %nm arbitrary
                    xest=estimators(obj.estimators(itr),scanout.phot,patternpos,L,fwhm);
                    
                    %recenter
                    dampf=2^(-obj.sequence_json.damping);
                    obj.pospattern=(1-dampf)*obj.pospattern+dampf*(xest); 
                end

                out.loc.xnm(loccounter,1)=xest(1);out.loc.ynm(loccounter,1)=xest(2);out.loc.znm(loccounter,1)=xest(3);
                out.loc.xfl(loccounter,1)=scanout.flpos(1)/scanout.counter;out.loc.yfl(loccounter,1)=scanout.flpos(2)/scanout.counter;out.loc.zfl(loccounter,1)=scanout.flpos(3)/scanout.counter;
                out.loc.time(loccounter,1)=scanout.time/scanout.counter;
                out.loc.itr(loccounter,1)=itr;
                out.loc.loccounter(loccounter,1)=loccounter;
                out.loc.cfr(loccounter,1)=cfr;
                out.loc.phot(loccounter,1)=photsum;
                out.loc.vld(loccounter,1)=1*(abortphot<4) && 2*(abortccr<4);
                out.loc.abortcondition(loccounter,1)=1*(abortphot>=4) + 2*(abortccr>=4);
                out.loc.patternrepeat(loccounter,1)=scanout.counter;
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