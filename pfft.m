function [X,freq]=pfft(x,Fs)
%   pfft: Returns the positive real portion of the fft of a 1D signal
%   [X,freq]=pfft(x,Fs)
%   x: initial signal
%   Fs: Signal sampling frequency
%   X: Fourier domain of series
%   Corresponding frequency spectrum of fourier domain series
%
%   2012 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute

%^^^
    % Detrend signal
    x=x-mean(x);
    N=length(x); %get the number of points
    freq=((0:N/2-1)/(N/2)*(Fs/2))';
    X=fft(x,N);
    X=abs(X(1:length(freq))).^2;
    
end

