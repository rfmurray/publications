function e = ice( cc )

% ICE  Illuminance contrast energy
% 
%     e = ice( cc )
% 
%     cc are the complex spherical harmonic coefficients of the light probe
%     e is the illuminance contrast energy

% blur by lambertian transfer function; this converts the light probe
% to the illuminance pattern that it generates on a sphere, as explained
% by Basri & Jacobs (2003), IEEE Trans. PAMI, 25, 218-233.
ccb = sphconvc( lambcoef, cc );

% find illuminance contrast energy directly from the spherical harmonic coefficients
ibar = abs(ccb(1))/sqrt(4*pi);   % mean illuminance
i2bar = sum(abs(ccb).^2)/(4*pi); % illuminance energy
e = sqrt( (i2bar-(ibar^2))/(ibar^2) );  % illuminance variance divided by squared mean illuminance

return


function [ c, mk ] = sphconvc( ck, cx )

% SPHCONVC  Spherical convolution on coefficients
%
% c = sphconvc( ck, cx )

% step through orders
c=cx(:);
Lmax=sqrt(size(c,1))-1;
for l=0:Lmax
    
    i=sphlm2i(l,0);		% index   of ith order, 0th degree term in ck
    i1=sphlm2i(l,-l);	% indices of ith order, +-ith degree terms in cx
    i2=sphlm2i(l,l);
    
    % check for nonzonal harmonics
    if any(ck(i1:i-1)) || any(ck(i+1:i2))
        warning('cannot convolve a function with nonzonal harmonics');
    end
    
    % get scale factors
    s=sqrt(4*pi/(2*l+1))*ck(i);
    mk(i1:i2,1)=s;
    
end

% scale coefficients
c=c.*mk;

% scale factors
mk=diag(mk);

return


function c = lambcoef

% LAMBCOEF  Spherical harmonic coefficients of Lambertian kernel
% 
% c = lambcoef

% initialize coefficient vector
L = 2;
c = zeros((L+1)^2,1);

% step through degrees
for i = 0:L

    % coefficients from Basri and Jacobs (2003), IEEE-PAMI, 25(2), equation (15)
    if i==0
        k = sqrt(pi)/2;
    elseif i==1
        k = sqrt(pi/3);
    elseif i/2==round(i/2)
        k = ((-1)^((i/2)+1)) * (factorial(i)/(factorial(i/2)*factorial(i-(i/2)))) * sqrt((2*i+1)*pi) / ((2^i)*(i-1)*(i+2));
    else
        k = 0;
    end

    % assign to zonal harmonic
    c( sphlm2i(i,0) ) = k;

end

return


function i = sphlm2i( l, m )

% SPHLM2I  Map spherical harmonic indices (l,m) to a single linearly ordered index i
% 
% i = sphlm2i( l, m )

if abs(m)>l
	i=NaN;
else
	i=l^2+l+m+1;
end

return
