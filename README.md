# ChenPhelan2021JFS

Replication code for
"International Coordination of Macroprudential Policies with Capital Flows and Financial Asymmetries" by
William Chen and Gregory Phelan, *Journal of Financial Stability*, accepted 2021.

The code in this repository produces the main figures and tables from saved output,
but it does *not* produce all the results (e.g. robustness checks) from scratch.
Several of the paper's results require the calculation of many different different
equilibria and are therefore done in parallel on a high-performance computing cluster.
We provide example scripts for these calculations.
However, the code for producing Appendices E.4 and E.5 is not provided
Furthermore, to generate most of the the paper's robustness results,
users will need to change parameters in the example scripts as needed.

If the user is interested in using this repository's code, then
looking at `numerical_example.m` will be a good starting point.
The key function for solving equilibrium is `eqm_functions/get_eqm.m`.
That function's source code will be helpful in learning how the
rest of the source code works since `get_eqm.m` is intended to be
the user-level interface with the underlying pseudo-transient
relaxation algorithm. Please file an issue if you have any questions.

## Directory Structure

* MATLAB scripts in top-level of directory: produce results for the baseline model,
robustness checks of the baseline model, and the extension to CES production.

* `eqm_functions/`: functions used to compute equilibrium and welfare for the
baseline model and the extension to CES production

* `parameters/`: scripts that initialize economic and numerical parameters

* `data/`: saved output

* `1good/`: functions and scripts to produce results when intermediate goods
are removed from the model

* `linear_invst/`: functions and scripts to produce results when the
investment technology is linear rather than strictly concave

## Tables and Figures Replication Instructions

The following instructions assumes output has already been generated
and the user only wants to produce tables and figures.

1. Run `make_paper_figs.m` in the top-level, `1good/`, and `linear_invst/`;
   and `dominant_strat.m` in `1good/` and `linear_invst/` to produce all the
   main tables and figures reported in the main text. Note that the numbers in
   tables were manually inputted, so these scripts will only generate the
   numbers reported, not the actual tables themselves.

2. The script `make_paper_figs` will sketch out the welfare surface
   for the AE and EME as a function of their leverage constraint policies.
   These welfare surfaces visually show the numerical result that the
   zero constraint is a dominant strategy.
   To further check the Nash result, run `analyze_check_zero.m` and
   `analyze_iternash.m`, which analyze output from the two other
   approaches for checking the Nash result.

3. Run `compare_stat_dist.m` to generate the counterfactual stationary distributions
   in Appendix B.3.

4. Run `bargaining_example.m` to obtain the numbers reported in the bargaining
   example of Appendix D.1.

5. Run `analyze_proddiff_robustness.m` and `analyze_risk_aversion_robustness.m`
   to obtain the numbers reported in Appendices E.1 - E.3.

6. Run `make_paper_ces_tbls_figs.m` to obtain the numbers reported in Appendix E.6
   on the extension to CES production.

7. Run `compare_ce_nash_specific_init_conds.m` to obtain the comparison
   of welfare in the competitive equilibrium to Nash at specific initial conditions.
   These numbers are discussed

## Output Replication Instructions

To re-generate all the saved output from scratch, follow these instructions.

1. Run `numerical_example.m` with `use_saved_guess = 0` to generate the
   converged solution in `data/gamA1p05_gamB2.mat`.

2. Run `compute_eqm_for_benchmarks_etc.m` and `compute_ces_eqms.m` to generate
   some output that is needed by `1good/` and `linear_invst/`

3. Change directory to `hpcc`.

4. Run `hpcc/plot_best_response.m` to sketch out the best responses
   and provide a coarse grid of initial guesses for the other scripts.

5. Run `plot_stat_pdf.m`, `get_best_response.m`, `sparse_get_br.m`, `check_zero.m`, and `iter_nash.m`
   in `hpcc/`.

6. Run `analyze_policy.m` twice after `get_best_response.m` and `sparse_get_br.m` finish, once each
   for the output of the latter two scripts.

7. To check the robustness results for different parameters, re-run steps 4-6 but either adding lines
   to the scripts to change parameters or creating a new parameter script instead of `parameters/baseline_parameters.m`.
