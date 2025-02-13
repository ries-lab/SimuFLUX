 classdef FlCollection<handle
    properties
        flall
        numberOfFluorophores
    end
    methods
        function add(obj,fllist)
            obj.flall=fllist;
            obj.numberOfFluorophores=length(fllist);
        end

        function [pos,isactive]=position(obj,time,props)
            flall=obj.flall;
            % fact=find(isactive);
            pos=zeros(length(flall),3);
            for k=1:length(flall)
                pos(k,:)=obj.flall{k}.position(time);
            end
            isactive=true(length(flall),1);
        end
        function ih=intensity(obj,intin,dt,time,phfac,props)
            flall=obj.flall;
            ih=zeros(length(flall),1);
            for k=1:length(flall)
                ih(k)=flall{k}.intensity(intin(k),dt,time,phfac(k));
            end

        end
        function reset(obj,time)
            %XXX reset fluorophore dynamics if needed
        end
        function updateonoff(obj,time) %dummy function not used here
        end
        function out=allbleached(obj)
            out=false;
        end
        function out=getproperties(obj)
            out=[];
        end
        function out=remainingphotons(obj)
            out=0;
            for k=1:length(obj.flall)
                out=out+(obj.flall{k}.remainingphotons);
            end
        end
        function out=brightness(obj)
            for k=length(obj.flall):-1:1
                out(k)=(obj.flall{k}.brightness);
            end
        end
    end
end
