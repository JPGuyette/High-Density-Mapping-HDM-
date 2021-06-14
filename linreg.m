function varargout=linreg(x,y)
% [m, b] = linreg(u,x)
% [m, b, R2] = linreg(u,x)
%   Outputs the slope of x w.r.t. y using a computationally efficient
%   matrix algebr approach.
%
%   INPUT:
%       x: x values
%       y: y values
%
%   OUTPUT:
%       m: Slope of line
%       b: y-intercept of line
%       R2: (Optional) R-squared value for goodness of fit
%
%

% Check input
    x=squeezecol(x);
    y=squeezecol(y);
    if numel(x)>length(x) || numel(y)>length(y)
        error('Inputs must be arrays')
    end
    
    if length(x) ~= length(y)
        error('Input vectors must be the same length')
    end
    
    
% Calculate necessary values
    n=numel(x);
    m=(n*sum(x.*y)-sum(x).*sum(y))/(n*sum(x.^2)-(sum(x)).^2);
    b=(1/n)*(sum(y)-m*sum(x));
    SSerr=sum((y-(m*x+b)).^2);
    SStot=sum((y-mean(y)).^2);
    R2=1-(SSerr/SStot);
    
    switch nargout
        case {0,1}
            out.m=m;
            out.b=b;
            out.R2=R2;
            varargout={out};
        case 2
            varargout={m,b};
        case 3
            varargout={m,b,R2};
        otherwise
            error('Function has too many outputs');
    end
end

function out=squeezecol(in)
% Forces 1d variable to be a column array
    in=squeeze(in);
    if size(in,2)>1
        in=in';
    end
    if size(in,2)>1
        error('Input must only have data in one dimension')
    end
    out=in;
end