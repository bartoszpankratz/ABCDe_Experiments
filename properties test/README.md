# ABCDe Speed Test

Codes used to generate the data for [properties test](https://github.com/bartoszpankratz/ABCDe_Experiments/blob/main/ABCD%20properties%20test%20plots.ipynb) of ABCD, LFR  and LFR NetworKit implementation.


### Requirements:

Main body of the experiment was done in Julia 1.6.1 with following packages:
- ABCDGraphGenerator v0.1.0
- ABCDeGraphGenerator v0.2.4
- BenchmarkTools v1.0.0
- Graphs v1.4.1
- StatsBase v0.33.12

Codes used for benchmarking LFR algorithm:
- modified version of [original LFR code](https://github.com/eXascaleInfolab/LFR-Benchmark_UndirWeightOvp) with extension allowing to use external degree and communitie size distribution files (code in <i>LFR-Benchmark_UndirWeightOvp_modified</i> folder).
- NetworKit 8.1 implementation in Python 3.8.8 

