% figure3b.m  Generate Figure 3b in Murray, Patel, and Yee (2015)
%             i.e., experiment 2, internal-to-external noise ratios

clear; clc;

% choose subject (subject = 1 to 22)
subject = 1;
fname = sprintf('data/experiment2_subject%02d.txt',subject);
fprintf('processing %s\n',fname);

% check for statistics toolbox
assert(~isempty(which('quantile')),'This program requires the MATLAB Statistics Toolbox');

% load data
trials = load(fname,'-ascii');

% column numbers
matchcol   = 2;    % trial number that this trial is a repeat of (or -1 if not a repeated trial)
sigcstcol  = 6;    % signal contrast
respintcol = 7;    % observer's response (1 = white dot first, 2 = black dot first)
correctcol = 8;    % correct response (1 = correct, 0 = incorrect)

% find two passes
pass2 = trials(trials(:,matchcol)>0,:);  % second pass
pass1 = trials(pass2(:,matchcol),:);     % first pass

% get signal contrast levels
cstlist = unique(pass1(:,sigcstcol))';

% get proportion correct and proportion of consistent responses at each signal contrast level
nrepeat  = [];  % number of pairs of repeated trials
ncorrect = [];  % number of correct responses
nagree   = [];  % number of consistent responses
for cst = cstlist
    i = pass1(:,sigcstcol)==cst;  % find trials at this signal contrast
    nrepeat(end+1)  = sum(i);     % count trials at this signal contrast
    ncorrect(end+1) = sum([ pass1(i,correctcol) ; pass2(i,correctcol) ]);  % count correct responses at this signal contrast
    nagree(end+1)   = sum( pass1(i,respintcol)==pass2(i,respintcol) );     % count consistent responses at this signal contrast
end

% convert counts to proportions
pcorrect = ncorrect./(2*nrepeat);  % proportion correct
pagree =   nagree./nrepeat;        % proportion of consistent responses

% bootstrap internal-to-external noise ratios at each signal level
ienoise = [];
for i = 1:numel(cstlist)
    
    % show progress
    fprintf('  signal contrast %.4f\n',cstlist(i));
    
    % bootstrap
    B = 10000;  % number of bootstrap iterations
    for b = 1:B+1
        
        % show progress
        if mod(b,50)==0, fprintf('    iteration %d of %d\n',b,B); end
        
        if b==1
            % on first iteration use actual pcorrect and pagree
            pcstar = pcorrect(i);
            pastar = pagree(i);
        else
            % on subsequent iterations resample pcorrect and pagree
            pcstar = binornd( 2*nrepeat(i), pcorrect(i) )/(2*nrepeat(i));
            pastar = binornd(   nrepeat(i), pagree(i) )/nrepeat(i);
        end
        
        % find internal-to-external noise ratio from pcstar and pastar
        ienoise(i,b) = pcpa2ie( pcstar, pastar );
        
    end
    
end

% find 95% confidence interval
q025 = quantile( ienoise(:,2:end), 0.025, 2 );
q975 = quantile( ienoise(:,2:end), 0.975, 2 );

% show internal-to-external noise ratios vs proportion correct
figure(1); clf;
errorbar( pcorrect, ienoise(:,1), q025, q975, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k' );
set(gca,'FontName','helvetica','FontWeight','bold','FontSize',18);
xlabel 'proportion correct'
ylabel 'internal-to-external noise ratio'
title(sprintf('subject %d',subject));
axis([ 0.48 1.02 0 2 ]);