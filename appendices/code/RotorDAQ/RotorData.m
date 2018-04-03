classdef RotorData < handle
    properties
        X %Horizontal Signal
        Y %Vertical Signal
        ref %Keyphasor Signal
        Fs %Sampling Rate (Samp/Sec)
        Fres %Frequency Resolution
        n %filter synchronous multiplier: 1X, 2X, 3X, ... nX the running speed
    end
    properties (SetAccess = private)
        Zf %Total Amplitude
        Phase %Phase
        Amp %Property that holds X and Y amplitudes
        Speed %RPM
    end
    properties (Access = private)
        refchop %%%%%%%%%%%%
        xchop   %% Organizes the data arrays into matrices based on windows
        ychop   %%%%%%%%%%%%
    end
    properties (Dependent, Access = private)
        NW %Number of windows
        freq %Frequency
        nspw %Number of samples per window
        
    end
    methods
        function obj = RotorData(x, y, r, f, fres)%Fills properties with data
            if nargin == 5
                obj.X = x;
                obj.Y = y;
                obj.ref = r;
                obj.Fs = f;
                obj.Fres = fres;
            end
        end
        function set.X(obj, val)
            obj.X = val;
        end
        function set.Y(obj, val)
            obj.Y = val;
        end
        function set.ref(obj, val)
            obj.ref = val;
        end
        function set.Fs(obj, val)
            obj.Fs = val;
            if isempty(obj.X) || isempty(obj.Y) || isempty(obj.ref) || isempty(obj.Fs) || isempty(obj.Fres)
            else
                obj.update;
            end
        end
        function set.Fres(obj, val)
            obj.Fres = val;
            if isempty(obj.X) || isempty(obj.Y) || isempty(obj.ref) || isempty(obj.Fs) || isempty(obj.Fres)
            else
                obj.update;
            end
        end
        function set.n(obj, val)
            obj.n = val;
            if isempty(obj.X) || isempty(obj.Y) || isempty(obj.ref) || isempty(obj.Fs) || isempty(obj.Fres)
            else
                obj.update;
            end
        end
        function nspw = get.nspw(obj)
            nspw = ceil(obj.Fs/obj.Fres);
        end
        function freq = get.freq(obj)
            df = obj.Fs/obj.nspw;
            f = df*(0:obj.nspw-1);
            Q = ceil((obj.nspw+1)/2); % M/2+1 for M even
            fQ = df*(Q-1);
            freq = f-fQ;
        end
        function NW = get.NW(obj)
            NW = floor(length(obj.X)/obj.nspw); %Rounds number of windows down to closest integer
        end
        function filter(obj)
            disp('filtering')
            for i = 1:1:obj.NW
                Fpass1 = obj.Speed(i)*obj.n/60-5; %Upper frequency of bandpass filter
                Fpass2 = obj.Speed(i)*obj.n/60+5; %Lower frequency of bandpass filter
                if i == 1
                    h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', 10, Fpass1, Fpass2, 50, .1, 50, obj.Fs);
                    Hd(i) = design(h, 'ellip');
                    Hd(i).persistentmemory = true;
                    obj.xchop(:,i) = filter(Hd(i),obj.xchop(:,i)); %Apply filter to data
                    xf = Hd(i).states;
                    reset(Hd(i));
                    obj.ychop(:,i) = filter(Hd(i),obj.ychop(:,i)); %||
                    yf = Hd(i).states;
                else
                    h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', 10, Fpass1, Fpass2, 50, .1, 50, obj.Fs);
                    Hd(i) = design(h, 'ellip');
                    Hd(i).persistentmemory = true;
                    Hd(i).states = xf;
                    obj.xchop(:,i) = filter(Hd(i),obj.xchop(:,i)); %Apply filter to data
                    xf = Hd(i).states;
                    Hd(i).states = yf;
                    obj.ychop(:,i) = filter(Hd(i),obj.ychop(:,i)); %||
                    yf = Hd(i).states;
                end
            end
        end
        function update(obj)
            obj.xchop = zeros(obj.nspw,obj.NW);
            obj.ychop = zeros(obj.nspw,obj.NW);
            obj.refchop = zeros(obj.nspw,obj.NW);
            obj.Speed = zeros(obj.NW,1);
            x = obj.X(1:obj.nspw*obj.NW);
            y = obj.Y(1:obj.nspw*obj.NW);
            r = obj.ref(1:obj.nspw*obj.NW);
            obj.xchop = reshape(x,obj.nspw,obj.NW);
            obj.ychop = reshape(y,obj.nspw,obj.NW);
            obj.refchop = reshape(r,obj.nspw,obj.NW);
            for i = 1:1:obj.NW
                pp = pulseperiod(obj.refchop(:,i),obj.Fs,'Tolerance',7);
                obj.Speed(i) = 1./mean(pp)*60;
            end
            if ~isempty(obj.n) && ~obj.n == 0
                filter(obj)
            end
            
            obj.Amp.XAmp = max(obj.xchop) - min(obj.xchop); %Calculates Amplitude by subtracting min and max from each column
            obj.Amp.YAmp = max(obj.ychop) - min(obj.ychop); %Calculates Amplitude by subtracting min and max from each column
            Zwin = hanning(obj.nspw, 'Periodic');
            Z = obj.xchop + 1i*obj.ychop;
            for i = 1:1:obj.NW
                Z(:,i) = Zwin.*Z(:,i);
            end
            obj.Zf = 2*abs(fft(Z)/obj.nspw);
            obj.Zf = fftshift(obj.Zf',2); %Shifts graph to center on the x axis
            
            obj.Phase.PhaseX = zeros(obj.NW,1);
            obj.Phase.PhaseY = zeros(obj.NW,1);
            distance = .5*obj.Fs/100;
            threshold = .5*(max(obj.refchop(:,ceil(end/2))) - mean(obj.refchop(:,ceil(end/2))));
            for i = 1:1:obj.NW
                [~, locsref] = findpeaks(obj.refchop(:,i),'MinPeakHeight',threshold,'MinPeakDistance',distance);
                [~, locsX] = findpeaks(obj.xchop(:,i),'MinPeakProminence',1,'MinPeakDistance',distance);
                [~, locsY] = findpeaks(obj.ychop(:,i),'MinPeakProminence',1,'MinPeakDistance',distance);
                if length(locsref) < 2 || length(locsX) < 2 || length(locsY) < 2
                    obj.Phase.PhaseX(i) = NaN;
                    obj.Phase.PhaseY(i) = NaN;
                else
                    if length(locsX) ~= length(locsY)
                        if length(locsX) < length(locsY)
                            if abs(locsX(1) - locsY(1)) < abs(locsX(end) - locsY(end))
                                locsY(end) = [];
                            else
                                locsY(1) = [];
                            end
                        else
                            if abs(locsX(1) - locsY(1)) < abs(locsX(end) - locsY(end))
                                locsX(end) = [];
                            else
                                locsX(1) = [];
                            end
                        end
                    end
                    if locsref(1) > locsX(1) || locsref(1) > locsY(1)
                        locsX(1) = [];
                        locsY(1) = [];
                    end
                    if length(locsref) > length(locsX)
                        locsref(length(locsX):end) = [];
                    end
                    phx = zeros(length(locsref)-1,1);
                    phy = zeros(length(locsref)-1,1);
                    for j = 1:1:length(locsref)-1
                        phx(j) = mod((locsX(j) - locsref(j))/(locsref(j+1) - locsref(j))*2*pi,2*pi);
                        phy(j) = mod((locsY(j) - locsref(j))/(locsref(j+1) - locsref(j))*2*pi,2*pi);
                    end
                    obj.Phase.PhaseX(i) = 180/pi*mod(mean(phx(:)),2*pi);
                    obj.Phase.PhaseY(i) = 180/pi*mod(mean(phy(:)),2*pi);
                    clear locsref locsX locsY j
                end
            end
            
        end
        function acqListener(obj, src, event)
            obj.newData = [event.TimeStamps, event.Data]';
            obj.ref = [obj.ref, event.Data(1,:)];
            obj.X = [obj.X, event.Data(2,:)];
            obj.Y = [obj.Y, event.Data(3,:)];
            obj.update;
            
        end
        %Plotting functions below
        bode(obj)
        cascade(obj)
        obj = downsample(obj,r)
        orbit3(obj)
        orbitAnimation(obj)
        mainplot(obj)
    end
end