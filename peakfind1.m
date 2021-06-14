function varargout=peakfind1(signal,Fs,passrange)
%   peakfind1 Finds peaks in a highly repetitive time sequence. This 
%   analysis is based on a fourier transform yielding the dominant
%   frequency in the time series. 
%
%   2013 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute

%^^^%     end
    if ~exist('passrange','var')

%     if ~exist('target','var')
%         target='both';
        passrange=[0 Inf];
    end
    % Remove linear trend from signal
    x=(1:length(signal))';
    [m,b]=linreg(x,signal);
    signal=signal-(m*x+b);
    
    % Remove DC component of signal
    signal=signal-mean(signal);
    
    % Take fourier transform and extract the dominant frequency
    [P,freq]=pfft(signal,Fs);
    P(freq<passrange(1)|freq>passrange(2))=0;
    maxfreq=freq(P==max(P));
    %disp(maxfreq)
    N=length(signal);
    win=2*round(Fs/(2*maxfreq));
    w=ceil(0.9*(win/2));
    
    % Zero pad signal
    signal2=[zeros(w,1);signal;zeros(w,1)];
    N2=length(signal2);
    
    
    maxindex=false(N,1);
    minindex=false(N,1);
    s=0;
    for i=1+w:1:N2-w
        j=i-w;
        if signal2(i)==max(signal2(i-w:i+w)) && s~=1
            maxindex(j)=true;
            s=1;
        end
        if signal2(i)==min(signal2(i-w:i+w)) && s~=-1
            minindex(j)=true;
            s=-1;
        end
    end
    peakindex=maxindex | minindex;
    switch nargout
        case {1,0}
            varargout{1}=peakindex;
        case 2
            varargout{1}=maxindex;
            varargout{2}=minindex;
        case 3
            varargout{1}=peakindex;
            varargout{2}=maxindex;
            varargout{3}=minindex;
    end
end