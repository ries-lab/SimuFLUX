classdef FlCollectionBlinking<FlCollection
    properties
        flprop %.moving, isactive, nexton, nextoff, remaining_activations
        switchpar=struct('starton',-1,'tonsmlm',1e3,'toffsmlm',1e5,'photonbudget',1000,'activations',1e6,'brightness',1000) %starton, tonsmlm, toffsmlm, microseconds
    end
    methods
        function [pos,isactive]=position(obj,time,props)
            flall=obj.flall;
            isactive=obj.flprop.isactive;
            % fact=find(isactive);
            % pos=zeros(length(fact),3);
            pos=zeros(length(flall),3);
            for k=1:length(flall)
                pos(k,:)=flall{k}.position(time);
            end
        end
        function ih=intensity(obj,intin,dt,time,phfac,props)
            fact=find(obj.flprop.isactive);
            ih=zeros(length(fact),1);
            for k=1:length(fact)
                ih(k)=obj.flall{fact(k)}.intensity(intin(k),dt,time,phfac(k));
            end
        end        
        function add(obj,fllist)
            if isa(fllist,"numeric") %call addstatic if poslist
                obj.addstatic(fllist);
                return
            end
            if ~iscell(fllist)
                fllist={fllist};
            end
             obj.flall=cat(2,obj.flall,fllist);
             numfl=length(fllist);
             flproph.isactive=true(1,numfl); flproph.nexton=zeros(1,numfl); flproph.nextoff=zeros(1,numfl);
             flproph.moving=true(1,numfl); flproph.activations=zeros(1,numfl); flproph.remaining_activations=zeros(1,numfl);
             obj.flprop=appendstruct(obj.flprop,flproph);
             obj.numberOfFluorophores=length(obj.flall);
        end
        function addstatic(obj,poslist)
            %arguments timestart: if obj.time different
            if size(poslist,2)==2 %2D data
                poslist(1,3)=0; %make 3D
            end
            numfl=size(poslist,1);
            switchpar=obj.switchpar;
            for k=numfl:-1:1
                flproph.moving(k)=false;
                if switchpar.starton>-1
                    flproph.isactive(k)=(switchpar.starton==1);
                else %random
                    flproph.isactive(k)=rand<switchpar.tonsmlm/(switchpar.toffsmlm+switchpar.tonsmlm);
                end
                
                addfl=FlBleach;
                addfl.brightness=switchpar.brightness;
                addfl.photonbudget=switchpar.photonbudget;
                addfl.pos=poslist(k,:);
                addfl.reset;

                flproph.nexton(k)=exprnd(switchpar.toffsmlm);
                flproph.nextoff(k)=exprnd(switchpar.tonsmlm);
                fllist{k}=addfl;
            end
            flproph.activations=double(flproph.isactive);
            flproph.remaining_activations=max(1,poissrnd(switchpar.activations,1,numfl));
            obj.flprop=appendstruct(obj.flprop,flproph);
            obj.flall=cat(2,obj.flall,fllist);

            obj.numberOfFluorophores=length(obj.flall);
        end
        function updateonoff(obj,time)
            flprop=obj.flprop;
            switchpar=obj.switchpar;
            isactive=flprop.isactive;

            %switch off 
            switchoff=time>=flprop.nextoff & flprop.isactive & ~flprop.moving;
            %add remaining photons=0?
            activef=find(isactive&~flprop.moving);
            for k=1:length(activef)
                if obj.flall{activef(k)}.remainingphotons<1
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
                    obj.flall{switchon(k)}.reset;
                end
            end
            if sson>0 || ssoff>0
                obj.flprop.isactive=isactive;  
            end
        end
        function setpar(obj,varargin)
            if nargin >1 %structure passed on
                spar=struct(varargin{:});
            else
                spar=varargin{1};
            end
            obj.switchpar=copyfields(obj.switchpar,spar);
                
        end
        % function reset(obj)
        %     %XXX reset fluorophore dynamics if needed
        % end
        function bl=allbleached(obj)
            bl=sum(obj.flprop.remaining_activations)==0;
        end
    end
end
