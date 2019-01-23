function ie = pcpa2ie( pcorrect, pagree )

% PCPA2IE  Convert proportion correct and proportion of consistent
%          responses to an internal-to-external noise ratio
% 
% ie = pcpa2ie( pcorrect, pagree )

% make the objective function
minfn = @( p ) errfn( p(1), p(2), pcorrect, pagree );

% minimize it
phat = fminsearch( minfn, [ 2*norminv(pcorrect) 1 ] );

% return the estimate of the internal-to-external noise ratio
% ( phat(1) is the estimate of d' )
ie = phat(2);

end


function ss = errfn( dprime, sigmaie, pcorrect, pagree )

% ERRFN  Find the sum-of-squares difference between the proportion
%        correct and proportion of consistent responses implied by
%        dprime and internal-to-external noise ratio ie, and
%        the observer's actual pcorrect and pagree

[ pccalc, pacalc ] = pcpa( dprime, sigmaie );
ss = (pcorrect-pccalc)^2 + (pagree-pacalc)^2;

end


function [ pc, pa, pAA, pBB ] = pcpa( dprime, sigmaie )

% PCPA  Find proportion correct and proportion agreement for
%       a maximum a posteriori observer
%
% [ pc, pa ] = pcpa( dprime, sigmaie )

% proportion correct
pc = normcdf(dprime/2);

% handle case sigmaie=0
if sigmaie==0
    pa = 1;
    return
end

% calculate derived parameters
sigmapi = sqrt( 1 + (1/sigmaie^2) );

% make coordinate vector
u = linspace( -5/sigmaie, 5/sigmaie, 1000 );
du = u(2)-u(1);

% integrals
pAA = trapz( (normcdf( dprime*sigmapi/2 - u ).^2)    .*normpdf(sigmaie*u) )*sigmaie*du;
pBB = trapz( ((1-normcdf( dprime*sigmapi/2 - u )).^2).*normpdf(sigmaie*u) )*sigmaie*du;
pa = pAA + pBB;

end
