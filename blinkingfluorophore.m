classdef blinkingfluorophore<MFfluorophore
    properties
   
        ton =100;
        toff =100;
        starton=false;
        % bleaching = true;
        photonbudget=inf; %1000 photons
        time=0;
        blinkingtrace
        tind=1;
        remainingphotons=1000;

    end
    methods
        % function obj=blinkingfluorophore(varargin)
        %     obj=obj@handle;
        %     obj.reset
        % end
        function Io=intensity(obj,I0,dwelltime)
            fraction=obj.measure(dwelltime);
            intensity=obj.brightness*I0*dwelltime*fraction;
            Io=min(intensity,obj.remainingphotons);
            obj.remainingphotons=obj.remainingphotons-Io;
        end
        function makeblinkingtrace(obj)
            obj.blinkingtrace=calculatetrace(obj.ton,obj.toff);
            obj.tind=1;
            obj.time=0;
        end
        function reset(obj) %new fluorophore
            obj.makeblinkingtrace;
            if ~obj.starton
                obj.blinkingtrace(:,3)=obj.blinkingtrace(:,3)-rand*obj.blinkingtrace(end,3)/2;
            end
            obj.blinkingtrace(:,3)=obj.blinkingtrace(:,3)+obj.time;
            % recalculate photon budget
            obj.remainingphotons=exprnd(obj.photonbudget);
        end
        function fraction=measure(obj,dwelltime)
            t=obj.time;
            trace=obj.blinkingtrace;

            if t>trace(end,3)-dwelltime*2 %outside range
                obj.extendblinkingtrace;% extend blinkingtrace;
            end
            ind1=obj.tind;
            % ind1=1;
            while trace(ind1,3)<=t %find last index below
                ind1=ind1+1;
            end
            ind1=ind1-1;
            ind2=ind1;
            while trace(ind2,3)<t+dwelltime
                ind2=ind2+1;
            end
            ind2=ind2-1;
            ontime=sum(trace(ind1:ind2-1,1));
            ontime=ontime-min(t-trace(ind1,3),trace(ind1,1)); %first block
            ontime=ontime+min(t+dwelltime-trace(ind2,3),trace(ind2,1)); % last block
            fraction=max(ontime/dwelltime,0);

            obj.tind=max(1,ind2-1);
            obj.time=t+dwelltime;

        end
        function extendblinkingtrace(obj)
            blt2=calculatetrace(obj.ton,obj.toff);
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