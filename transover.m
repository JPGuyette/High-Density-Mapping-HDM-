function [raw,pix,im]=transover(ol_image,ol_roi,overlay,trans,pos,range)
% transover Creates overlays a matrix on an image given the location of the
% matric
%
%   ol_image: Matlab variable containing the image to overlay
%   ol_roi: xxyy position of the region of interest on the image (pixels)
%   overlay: Matrix to overlay in the ROI
%   trans: Opacity of the overlay
%   pos: The region to actually display
%   range: The strain range to show
%
%   2013 John Favreau
%   The Gaudette Lab at Gateway Park
%   Worcester Polytechnic Institute

%^^^

% Create and resize a pixel-sized interpolation of the matrix image
        raw = overlay;
        pixyx=[1+ol_roi(4)-ol_roi(3),1+ol_roi(2)-ol_roi(1)];
        matyx=size(raw);
        maxratio=max(pixyx/matyx);
        n=nextpow2(maxratio);
        pix = interp2(raw,n);
        pix=imresize(pix,pixyx);
        
    % Create an overlay for just the roi
        ol_mask=poly2mask(pos(:,1)-ol_roi(1)+1,pos(:,2)-ol_roi(3)+1,...
            1+ol_roi(4)-ol_roi(3),1+ol_roi(2)-ol_roi(1));
        
        pix=pix.*ol_mask;
        im=pix;
    
    % Calculate a pixel-sized rgb image for the overlay
        
        if ~exist('range','var');
            minval = min(min(im));
            maxval=max(max(im));
        else
            minval=range(1);
            maxval=range(2);
        end
        im=im-minval;
        im=im/(maxval-minval);
        im=im*255;
        im=im+1;
        im=round(im);
        im=ind2rgb(im,jet(256));
        close gcf

    % Construct strain overlay map and place in the ROI
        ol_image=(double(ol_image)/255);
        ol_image=repmat(ol_image,[1 1 3]);
        im_roi=ol_image(ol_roi(3):ol_roi(4),ol_roi(1):ol_roi(2),:);
        c_a=trans*im;
        C=c_a+(1-trans)*im_roi;
        [imheight,imwidth,~]=size(ol_image);
        fullmask=poly2mask(pos(:,1),pos(:,2),imheight,imwidth);
        mask=repmat(fullmask(ol_roi(3):ol_roi(4),ol_roi(1):ol_roi(2),:),[1 1 3]);
        C(~mask)=im_roi(~mask);
        ol_image(ol_roi(3):ol_roi(4),ol_roi(1):ol_roi(2),:)=C;
        cbar=getcolorbar(range,size(ol_image,1));
        im=zeros(imheight,imwidth+size(cbar,2),3);
        im(:,1:imwidth,:)=ol_image;
        im(:,imwidth+1:end,:)=cbar;
end

function scale=getcbar(min,max)

end