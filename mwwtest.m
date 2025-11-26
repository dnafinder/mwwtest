function STATS = mwwtest(x1,x2)
%MWWTEST Mann-Whitney-Wilcoxon test for two unpaired samples.
%
%   Syntax
%   ------
%   STATS = MWWTEST(X1, X2)
%
%   Description
%   -----------
%   MWWTEST performs the non-parametric Mann-Whitney-Wilcoxon test
%   (also known as Wilcoxon rank-sum test) to compare two independent
%   (unpaired) samples X1 and X2.
%
%   The two samples are pooled, ranked (with tie correction), and the
%   Wilcoxon rank-sum statistics W1, W2 and the equivalent Mann-Whitney
%   statistics U1, U2 are computed. Depending on the number of possible
%   rank combinations, the function:
%
%     * uses the exact null distribution of W/U when the number of
%       combinations is < 20,000;
%     * otherwise uses a normal approximation (with tie correction
%       when required).
%
%   The numerical results are equivalent to MATLAB's RANKSUM function,
%   but MWWTEST provides additional detail and formatted tabular output.
%
%   Inputs
%   ------
%   X1 : numeric vector (row or column), real, finite, non-NaN, non-empty.
%   X2 : numeric vector (row or column), real, finite, non-NaN, non-empty.
%
%   Outputs
%   -------
%   STATS : structure with the following fields
%
%     STATS.n   : [n1 n2] sample sizes.
%     STATS.W   : [W1 W2] sums of ranks for group 1 and group 2.
%     STATS.mr  : [mr1 mr2] mean ranks.
%     STATS.U   : [U1 U2] Mann-Whitney U statistics.
%
%     If the normal approximation is used:
%       STATS.method = 'Normal approximation'
%       STATS.mU     : mean of U under H0.
%       STATS.sU     : standard deviation of U under H0.
%       STATS.Z      : Z statistic.
%       STATS.p      : [p_one_tail p_two_tails].
%
%     If the exact distribution is used:
%       STATS.method = 'Exact distribution'
%       STATS.T      : W statistic for the smaller group.
%       STATS.p      : [p_one_tail p_two_tails].
%
%   Notes
%   -----
%   * If the function is called without output argument, results are
%     displayed in the Command Window as formatted tables.
%   * When an output is requested, the same information is stored in
%     the STATS structure.
%
%   Example
%   -------
%   X1 = [181 183 170 173 174 179 172 175 178 176 158 179 180 172 177];
%   X2 = [168 165 163 175 176 166 163 174 175 173 179 180 176 167 176];
%
%   STATS = mwwtest(X1, X2);
%
%   See also
%   --------
%   RANKSUM, TIEDRANK
%
%   References
%   ----------
%   Cardillo G. (2009). MWWTEST: Mann-Whitney-Wilcoxon non parametric
%   test for two unpaired samples. MATLAB Central File Exchange.
%
%   ------------------------------------------------------------------
%   Author : Giuseppe Cardillo
%   Email  : giuseppe.cardillo.75@gmail.com
%   GitHub : https://github.com/dnafinder/mwwtest
%   Created: 2009
%   Updated: 2025-11-26
%   Version: 2.0.0
%   ------------------------------------------------------------------

% Input error handling
ip = inputParser;
ip.FunctionName = 'mwwtest';

addRequired(ip,'x1',@(x) validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','nonnan','nonempty'}));
addRequired(ip,'x2',@(x) validateattributes(x,{'numeric'}, ...
    {'vector','real','finite','nonnan','nonempty'}));

parse(ip,x1,x2);

% Force column vectors internally
x1 = x1(:);
x2 = x2(:);

% Set the basic parameters
n1 = numel(x1);
n2 = numel(x2);
NP = n1*n2;      % total number of pairwise comparisons
N  = n1+n2;      % total sample size
N1 = N+1;
k  = min([n1 n2]);

% Compute the ranks and the tie correction
% A = ranks of pooled data
% B = tie adjustment (scalar, used for variance correction)
[A,B] = tiedrank([x1; x2]);

R1 = A(1:n1);
R2 = A(n1+1:end);

T1 = sum(R1);
T2 = sum(R2);

% U statistics (equivalent Mann-Whitney formulation)
U1 = NP + (n1*(n1+1))/2 - T1;
U2 = NP - U1;

% Display header and basic descriptive table
disp('MANN-WHITNEY-WILCOXON TEST')
disp(' ')

summaryTable = table(...
    [median(x1); n1; T1; T1/n1; U1], ...
    [median(x2); n2; T2; T2/n2; U2], ...
    'VariableNames', {'Group_1','Group_2'}, ...
    'RowNames', {'Median','Numerosity','Sum_of_Rank_W', ...
                 'Mean_Rank','Test_variable_U'});

disp(summaryTable)

% Initialize STATS if requested
if nargout
    STATS.n  = [n1 n2];
    STATS.W  = [T1 T2];
    STATS.mr = [T1/n1 T2/n2];
    STATS.U  = [U1 U2];
end

% Decide between exact distribution and normal approximation
combCount = round(exp(gammaln(N1) - gammaln(k+1) - gammaln(N1-k)));

if combCount > 20000
    % --- Normal approximation branch ---
    mU = NP/2; % mean of U under H0
    
    if B == 0
        % No ties: standard variance of U
        sU = sqrt(NP * N1 / 12);
    else
        % With ties: tie-corrected variance of U
        sU = sqrt((NP/(N^2 - N)) * ((N^3 - N - 2*B)/12));
    end
    
    % Continuity-corrected Z statistic
    Z1 = (abs(U1 - mU) - 0.5) / sU;
    p_one_tail  = 1 - normcdf(Z1);
    p_two_tails = min(1, 2*p_one_tail);
    
    disp('Sample size is large enough to use the normal distribution approximation')
    disp(' ')
    
    resTable = table(mU, sU, Z1, p_one_tail, p_two_tails, ...
        'VariableNames', {'Mean','SD','Z','p_value_one_tail','p_value_two_tails'});
    disp(resTable)
    
    if nargout
        STATS.method = 'Normal approximation';
        STATS.mU     = mU;
        STATS.sU     = sU;
        STATS.Z      = Z1;
        STATS.p      = [p_one_tail p_two_tails];
    end
else
    % --- Exact distribution branch ---
    disp('Sample size is small enough to use the exact Mann-Whitney-Wilcoxon distribution')
    disp(' ')
    
    % W for the smaller group (conventional choice)
    if n1 <= n2
        w = T1;
    else
        w = T2;
    end
    
    % Exact distribution of W by enumerating all combinations
    pdf = sum(nchoosek(A, k), 2); % sums of ranks for all k-combinations
    
    P = [sum(pdf <= w) sum(pdf >= w)] ./ numel(pdf);
    p_one_tail  = min(P);
    p_two_tails = min(1, 2*p_one_tail);
    
    resTable = table(w, p_one_tail, p_two_tails, ...
        'VariableNames', {'W','p_value_one_tail','p_value_two_tails'});
    disp(resTable)
    
    if nargout
        STATS.method = 'Exact distribution';
        STATS.T      = w;
        STATS.p      = [p_one_tail p_two_tails];
    end
end

disp(' ')
end
