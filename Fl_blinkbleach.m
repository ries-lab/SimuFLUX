classdef Fl_blinkbleach<Fl_Fluorophore
    properties
   
        fast_ton =100;
        fast_toff =0;
        starton=false;
        photonbudget=inf; %1000 photons
        time=0;
        blinkingtrace=[0 0 0];
        tind=1;
        

    end
    methods
        function Io=intensity(obj,I0,dwelltime,time,phfac,varargin)

            fraction=obj.measure(dwelltime,time);
            intensity=(obj.brightness/1000)*I0*dwelltime*fraction;
            remainingphotons=obj.remainingphotons;
            Io=min(intensity,remainingphotons);
            obj.remainingphotons=remainingphotons-Io/phfac;
        end
        function makeblinkingtrace(obj)
            obj.blinkingtrace=calculatetrace(obj.fast_ton,obj.fast_toff);
            obj.tind=1;
            obj.time=0;
        end
        function reset(obj,time) %new fluorophore
            
            obj.makeblinkingtrace;
            if ~obj.starton
                obj.blinkingtrace(:,3)=obj.blinkingtrace(:,3)-rand*obj.blinkingtrace(end,3)/2;
            end
            obj.blinkingtrace(:,3)=obj.blinkingtrace(:,3)+time;
            % obj.blinkingtrace(:,3)=obj.blinkingtrace(:,3)+obj.time;
            % recalculate photon budget
            obj.remainingphotons=exprnd(obj.photonbudget);
        end
        function fraction=measure(obj,dwelltime,t)
            % t=obj.time;
            trace=obj.blinkingtrace;
            
            if t>trace(end,3)-dwelltime*2 %outside range
                obj.extendblinkingtrace;% extend blinkingtrace;
            end
            ind1=obj.tind;
            if trace(ind1,3)>t %time asked that was before
                ind1=find(trace(ind1,3)<=t,1,'last');
                if isempty(ind1)
                    ind1=1;
                end
            end
            % ind1=1;
            while trace(ind1,3)<=t %find last index below
                ind1=ind1+1;
            end
            ind1=max(1,ind1-1);
            ind2=ind1;
            while trace(ind2,3)<t+dwelltime
                ind2=ind2+1;
            end
            ind2=max(1,ind2-1);
            ontime=sum(trace(ind1:ind2-1,1));
            ontime=ontime-min(t-trace(ind1,3),trace(ind1,1)); %first block
            ontime=ontime+min(t+dwelltime-trace(ind2,3),trace(ind2,1)); % last block
            fraction=max(ontime/dwelltime,0);

            obj.tind=ind2;
            % obj.time=t+dwelltime;

        end
        function extendblinkingtrace(obj)
            blt2=calculatetrace(obj.fast_ton,obj.fast_toff);
            blt2(:,3)=blt2(:,3)+sum(obj.blinkingtrace(end,:));
            obj.blinkingtrace=vertcat(obj.blinkingtrace,blt2);
        end
    end
end


function blinkingtrace=calculatetrace(ton,toff)
    maxsimultime=100*(ton+toff);
    numpoints=max(10,ceil(maxsimultime/(ton+toff)*2));
    tont=exprnd(ton,numpoints,1);
    tofft=exprnd(toff,numpoints,1);
    tsum=tont+tofft;
    timestart=cumsum([0;tsum(1:end-1)]);
    blinkingtrace=horzcat(tont,tofft,timestart);
end