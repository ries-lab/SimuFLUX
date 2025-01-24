classdef blinkingFluorophorecollection<Fluorophorecollection
    properties
        % flall
        % numberOfFluorophores
        flprop %.moving, isactive, nexton, nextoff, remaining_activations
        switchpar=struct('starton',-1,'tonsmlm',1e3,'toffsmlm',5e4,'photonbudget',1000,'activations',1e6,'brightness',1000) %starton, tonsmlm, toffsmlm, microseconds
        %starton: 0 or 1 (off or on), -1: random (given by ton, toff)
        % tprevious=0;
    end
    methods
        function addmovingfl(obj,fllist)
        end
        function addstatic(obj,poslist)
            %arguments timestart: if obj.time different
            if size(poslist,2)==2 %2D data
                poslist(1,3)=0; %make 3D
            end
            numfl=size(poslist,1);
            switchpar=obj.switchpar;
            for k=1:numfl
                moving(k)=false;
                if switchpar.starton>-1
                    isactive(k)=(switchpar.starton==1);
                else %random
                    isactive(k)=rand<switchpar.tonsmlm/(switchpar.toffsmlm+switchpar.tonsmlm);
                end
                
                addfl=bleachingFluorophore;
                addfl.brightness=switchpar.brightness;
                addfl.photonbudget=switchpar.photonbudget;
                % addfl.activations=switchpar.activations;
                addfl.pos=poslist(k,:);
                addfl.reset;
                nexton(k)=exprnd(switchpar.toffsmlm);
                nextoff(k)=exprnd(switchpar.tonsmlm);
                flall(k)=addfl;
            end
            if isempty(obj.flall) %initialize
                obj.flall=flall;
                obj.flprop.moving=moving;
                obj.flprop.isactive=isactive;
                obj.flprop.nexton=nexton;
                obj.flprop.nextoff=nextoff;
                obj.flprop.activations=double(isactive);
                obj.flprop.remaining_activations=max(1,poissrnd(switchpar.activations,1,numfl));
            else
                obj.flall(end+1:end+numfl)=flall;
                obj.flprop.moving(end+1:end+numfl)=moving;
                obj.flprop.isactive(end+1:end+numfl)=isactive;
                obj.flprop.activations(end+1:end+numfl)=double(isactive);
                obj.flprop.nexton(end+1:end+numfl)=nexton;
                obj.flprop.nextoff(end+1:end+numfl)=nextoff;
                obj.flprop.remaining_activations(end+1:end+numfl)=max(1,poissrnd(switchpar.activations,1,numfl));
            end
            obj.numberOfFluorophores=length(obj.flprop.isactive);
        end
        function updateonoff(obj,time)
            % dt=time-obj.tprevious;
            
            flprop=obj.flprop;
            switchpar=obj.switchpar;
            isactive=flprop.isactive;

            %switch off 
            switchoff=time>=flprop.nextoff & flprop.isactive & ~flprop.moving;
            %add remaining photons=0?
            activef=find(isactive&~flprop.moving);
            for k=1:length(activef)
                if obj.flall(activef(k)).remainingphotons<1
                    switchoff(activef(k))=true;
                end
            end
            ssoff=sum(switchoff);

            if ssoff>0            
                isactive(switchoff)=false;
                obj.flprop.nexton(switchoff)=exprnd(switchpar.toffsmlm,ssoff,1)+time;
            end

            %switch on
            switchon=find(time>=flprop.nexton & ~flprop.isactive & ~flprop.moving & flprop.remaining_activations>0);
            sson=length(switchon);
            if sson>0
                isactive(switchon)=true;
                obj.flprop.activations(switchon)=obj.flprop.activations(switchon)+1;
                obj.flprop.nextoff(switchon)=exprnd(switchpar.tonsmlm,sson,1)+time;
                obj.flprop.remaining_activations(switchon)=obj.flprop.remaining_activations(switchon)-1;
                for k=1:sson
                    obj.flall(switchon(k)).reset;
                end
            end
            if sson>0 || ssoff>0
                obj.flprop.isactive=isactive;  
            end
            % sum(isactive)
        end
        function [pos,isactive]=position(obj,time)
            isactive=obj.flprop.isactive;
            fact=find(isactive);
            pos=zeros(length(fact),3);
            for k=1:length(fact)
                pos(k,:)=obj.flall(fact(k)).position(time);
            end
        end
        function ih=intensity(obj,intin,dt,phfac)
            fact=find(obj.flprop.isactive);
            ih=zeros(length(fact),1);
            for k=1:length(fact)
                ih(k)=obj.flall(fact(k)).intensity(intin(k),dt,phfac(k));
            end

        end
        function setpar(obj,switchpar)
            obj.switchpar=copyfields(obj.switchpar,switchpar);
        end
        % function reset(obj)
        %     %XXX reset fluorophore dynamics if needed
        % end
        function bl=allbleached(obj)
            bl=sum(obj.flprop.remaining_activations)==0;
        end
    end
end
