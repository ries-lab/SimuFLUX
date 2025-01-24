classdef Fluorophorecollection<handle
    properties
        flall
        numberOfFluorophores
        % flprop %.moving, isactive, nexton, nextoff, remaining_activations
        % switchpar=struct('starton',-1,'tonsmlm',1e3,'toffsmlm',5e4,'photonbudget',1000,'activations',1e6) %starton, tonsmlm, toffsmlm, microseconds
        %starton: 0 or 1 (off or on), -1: random (given by ton, toff)
        % tprevious=0;
    end
    methods
        function add(obj,fllist)
            obj.flall=fllist;
            obj.numberOfFluorophores=length(fllist);
        end
        
       
        function [pos,isactive]=position(obj,time)
            
            flall=obj.flall;
            % fact=find(isactive);
            pos=zeros(length(flall),3);
            for k=1:length(flall)
                pos(k,:)=obj.flall(k).position(time);
            end
            isactive=true(length(flall),1);
        end
        function ih=intensity(obj,intin,dt,phfac)
            flall=obj.flall;
            ih=zeros(length(flall),1);
            for k=1:length(flall)
                ih(k)=obj.flall(k).intensity(intin(k),dt,phfac(k));
            end

        end
        function reset(obj)
            %XXX reset fluorophore dynamics if needed
        end
        function updateonoff(obj,time) %dummy function not used here
        end
        function out=allbleached(obj)
            out=false;
        end
        % function bl=allbleached(obj)
        %     bl=sum(obj.flprop.remaining_activations)==0;
        % end
    end
end
