function [ dprime, hit, fa ] = idealdprime( rmu, rsigma, rclip, lightvec1, lightvec2, atratio, nvec )

% IDEALDPRIME  Simulate ideal performance in lighting direction discrimination task
%
%     [ dprime, hit, fa ] = idealdprime( rmu, rsigma, rclip, lightvec1, lightvec2, atratio, nvec )
% 
%     rmu is the mean of the reflectance distribution
%     rsigma is the nominal standard deviation of the reflectance distribution
%     rclip is the number of nominal standard deviations from the mean at which the reflectance distribution should be clipped
%     lightvec1 is the lighting direction in the first interval
%     lightvec2 is the lighting direction in the second interval
%     atratio is the ambient-to-total light ratio, i.e., eambient/(eambient+edirect)
%     nvec is an n x 3 matrix of the unit normal vectors of the polyhedron's faces

% find direct and ambient illuminances
% - uses fact that atratio = eambient/(edirect+eambient)
% - assumes edirect+eambient = 1 (can assume any value here without
%   changing ideal performance) and solves for edirect and eambient
edirect  = (1-atratio);
eambient = atratio;

% count normal vectors
facen = size(nvec,1);

% set number of trials to simulate
T = 10000;  

% initialize variables
targeti = NaN(T,1);
responsei = targeti;
dvar1 = targeti;
dvar2 = targeti;

% run trials
for t = 1:T
    
    % show progress
    if mod(t,1000)==0
        fprintf('idealdprime %d / %d\n',t,T);
    end
    
    % choose order of lighting conditions
    targeti(t) = 1 + (rand>0.5);  % 1 if lighting order is (1,2), and 2 if lighting order is (2,1)
    if targeti(t) == 1
        lvecA = lightvec1;  % lighting direction in interval A
        lvecB = lightvec2;  % lighting direction in interval B
    else
        lvecA = lightvec2;  % lighting direction in interval A
        lvecB = lightvec1;  % lighting direction in interval B
    end
    
    % assign random reflectances to faces
    refA = cnormrnd(rmu,rsigma,rclip,[ facen 1 ]);  % reflectances in interval A
    refB = cnormrnd(rmu,rsigma,rclip,[ facen 1 ]);  % reflectances in interval B
    
    % render (i.e., convert reflectance and illuminance to luminance)
    lumA = (refA/pi).*( edirect*max(nvec*lvecA',0) + eambient );  % stimulus luminances in interval A
    lumB = (refB/pi).*( edirect*max(nvec*lvecB',0) + eambient );  % stimulus luminances in interval B

    % calculate the illuminances under each of the lighting conditions
    illum1 = edirect*max(nvec*lightvec1',0) + eambient;  % illuminances created by lighting direction 1
    illum2 = edirect*max(nvec*lightvec2',0) + eambient;  % illuminances created by lighting direction 2

    % recover reflectances under the two lighting hypotheses
    refA1 = pi*lumA./illum1;  % reflectances in interval A, if it's illuminated by lighting direction 1
    refA2 = pi*lumA./illum2;  % reflectances in interval A, if it's illuminated by lighting direction 2
    refB1 = pi*lumB./illum1;  % reflectances in interval B, if it's illuminated by lighting direction 1
    refB2 = pi*lumB./illum2;  % reflectances in interval B, if it's illuminated by lighting direction 2
    
    % calculate log probabilities of luminances under the two lighting hypotheses
    lpA1 = sum( log(cnormpdf(refA1,rmu,rsigma,rclip)) - log(illum1) );  % likelihood of luminances in interval A, if it's illuminated from lighting direction 1
    lpA2 = sum( log(cnormpdf(refA2,rmu,rsigma,rclip)) - log(illum2) );  % likelihood of luminances in interval A, if it's illuminated from lighting direction 2
    lpB1 = sum( log(cnormpdf(refB1,rmu,rsigma,rclip)) - log(illum1) );  % likelihood of luminances in interval B, if it's illuminated from lighting direction 1
    lpB2 = sum( log(cnormpdf(refB2,rmu,rsigma,rclip)) - log(illum2) );  % likelihood of luminances in interval B, if it's illuminated from lighting direction 2
    
    % respond
    dvar1(t) = lpA1 + lpB2;  % log likelihood that lighting order is (1,2)
    dvar2(t) = lpA2 + lpB1;  % log likelihood that lighting order is (2,1)
    if dvar1(t)==dvar2(t)
        responsei(t) = 1 + (rand<0.5);
    else
        responsei(t) = 1 + (dvar2(t)>dvar1(t));  % 1 if most likely lighting order is (1,2), and 2 if it's (2,1)
    end
    
end

% calculate sensitivity
hit = mean( responsei(targeti==1)==1 );
fa  = mean( responsei(targeti==2)==1 );
dprime = (z(hit)-z(fa))/sqrt(2);

end


function y = z( x )

% Z  Inverse of standard normal cumulative distribution function
% 
%     y = z( x )

y = sqrt(2)*erfinv(2*(x-0.5));

end


function y = cnormpdf( x, mu, sigma, nclip )

% CNORMPDF  Probability density function of the clipped normal distribution
% 
%     y = cnormpdf( x, mu, sigma, nclip )
% 
%     x is the value of the random variable at which the density is to be evaluated
%     mu is the mean of the clipped normal distribution
%     sigma is the nominal standard deviation of the clipped normal distribution
%     nclip is the number of nominal standard deviations from the mean at which the distribution is clipped

% convert to z-scores
z=(x-mu)./sigma;

% calculate pdf
ptail=0.5*erfc(nclip/sqrt(2));
y=exp(-0.5*z.^2)./((1-2*ptail)*sqrt(2*pi).*sigma);
y(abs(z)>nclip)=0;

end


function x = cnorminv( p, mu, sigma, nclip )

% CNORMINV  Inverse cumulative distribution function of the clipped
%           normal distribution
% 
%     x = cnorminv( p, mu, sigma, nclip )
% 
%     p is the proportion of the distribution that lies below the return value x
%     mu is the mean of the clipped normal distribution
%     sigma is the nominal standard deviation of the clipped normal distribution
%     nclip is the number of nominal standard deviations from the mean at which the distribution is clipped

% calculate inverse cdf
ptail=0.5*erfc(nclip/sqrt(2));
pstar=ptail+p.*(1-2*ptail);
z=-sqrt(2)*erfcinv(2*pstar);
x=mu+sigma.*z;

end


function r = cnormrnd( mu, sigma, nclip, dim )

% CNORMRND  Random number generator for the clipped normal distribution
% 
%     r = cnormrnd( mu, sigma, nclip, dim )
% 
%     mu is the mean of the clipped normal distribution
%     sigma is the nominal standard deviation of the clipped normal distribution
%     nclip is the number of nominal standard deviations from the mean at which the distribution is clipped
%     dim is the dimension of the matrix of random values to be returned

% generate random numbers using the inverse cumulative distribution function
r = cnorminv( rand(dim), mu, sigma, nclip );

end
