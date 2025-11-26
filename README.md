[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/mwwtest&file=mwwtest.m)

# mwwtest

## 📌 Overview

`mwwtest` performs the non-parametric Mann–Whitney–Wilcoxon test (also known as the Wilcoxon rank-sum test) to compare two independent (unpaired) samples.

The function reproduces the numerical results of MATLAB’s built-in `ranksum`, while providing:

- detailed tabular output in the Command Window
- a structured `STATS` output with all key statistics

Depending on the sample sizes, `mwwtest` automatically switches between the exact null distribution and a normal approximation (with tie correction when needed).

---

## 📐 Syntax

STATS = mwwtest(X1, X2)

If no output argument is requested:

mwwtest(X1, X2)

the results are displayed only in the Command Window as formatted tables.

---

## 📥 Inputs

- X1 – numeric vector (row or column), real, finite, non-NaN, non-empty
- X2 – numeric vector (row or column), real, finite, non-NaN, non-empty

Internally, both inputs are reshaped as column vectors; the statistical results are unaffected by the original orientation.

---

## 📤 Outputs

### STATS structure

STATS is a structure with the following fields:

- STATS.n  
  [n1 n2] – sample sizes of X1 and X2.

- STATS.W  
  [W1 W2] – sums of ranks for group 1 and group 2.

- STATS.mr  
  [mr1 mr2] – mean ranks for group 1 and group 2.

- STATS.U  
  [U1 U2] – Mann–Whitney U statistics for group 1 and group 2.

Depending on the branch used:

#### When the normal approximation is used

- STATS.method = 'Normal approximation'
- STATS.mU – mean of U under the null hypothesis
- STATS.sU – standard deviation of U under the null hypothesis (with tie correction if necessary)
- STATS.Z – continuity-corrected Z statistic
- STATS.p – [p_one_tail p_two_tails]

#### When the exact distribution is used

- STATS.method = 'Exact distribution'
- STATS.T – W statistic (sum of ranks for the smaller group)
- STATS.p – [p_one_tail p_two_tails]

---

## 📊 Example

X1 = [181 183 170 173 174 179 172 175 178 176 158 179 180 172 177];
X2 = [168 165 163 175 176 166 163 174 175 173 179 180 176 167 176];

STATS = mwwtest(X1, X2);

Typical Command Window output (abbreviated):

MANN-WHITNEY-WILCOXON TEST

                      Group_1    Group_2
                      _______    _______

    Median              ...        ...
    Numerosity          ...        ...
    Sum_of_Rank_W       ...        ...
    Mean_Rank           ...        ...
    Test_variable_U     ...        ...

Sample size is large enough to use the normal distribution approximation

    Mean       SD        Z       p_value_one_tail    p_value_two_tails
    _____    ______    ______    ________________    _________________
      .         .         .              .                    .

The same information is stored in the STATS structure for further processing.

---

## 🧠 Method

1. Ranking and ties

   - The two samples are pooled and ranked using tiedrank.
   - The function obtains:
     - ranks for each group (W1, W2)
     - the tie-adjustment term required for variance correction.

2. U statistics

   - From the rank sums, the Mann–Whitney statistics U1 and U2 are computed using the standard formulas.

3. Exact vs normal decision

   - Let n1 = length(X1), n2 = length(X2), N = n1 + n2, and k = min(n1, n2).
   - The number of combinations C(N, k) is evaluated (via gammaln).
   - If C(N, k) < 20,000:
     - the exact distribution of the rank sum is enumerated by evaluating all k-combinations of ranks;
   - otherwise:
     - a normal approximation is used with:
       - mean mU = n1*n2/2
       - variance corrected for ties when necessary
       - continuity-corrected Z statistic
       - one- and two-sided p-values.

4. Output

   - A formatted table is printed in the Command Window.
   - If an output argument is requested, all relevant quantities are stored in STATS.

---

## 📦 Requirements

- MATLAB (tested on recent releases; uses inputParser, validateattributes, tiedrank, normcdf, gammaln, nchoosek, table)
- No additional toolboxes are required.

---

## 📚 References

- Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. Annals of Mathematical Statistics, 18(1), 50–60.
- Wilcoxon, F. (1945). Individual comparisons by ranking methods. Biometrics Bulletin, 1(6), 80–83.
- MATLAB documentation for ranksum and tiedrank.

Original File Exchange reference:

Cardillo G
