# ChenPhelan2021JFS

Replication code for
"International Coordination of Macroprudential Policies with Capital Flows and Financial Asymmetries" by
William Chen and Gregory Phelan, *Journal of Financial Stability*, accepted 2021.

The code in this repository produces the main figures and tables from saved output,
but it does *not* produce all the results (e.g. robustness checks) from scratch.
Several of the paper's results require the calculation of many different different
equilibria and are therefore done in parallel on a high-performance computing cluster.
We provide example scripts for these calculations. To generate all the paper's results,
users will need to change parameters in the example scripts as needed.
