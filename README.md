# ABCDe Experiments

Notebooks and codes designed to test the speed and properties of [ABCD](https://github.com/bkamins/ABCDGraphGenerator.jl),[ABCDe](https://github.com/tolcz/ABCDeGraphGenerator.jl/) and [LFR](https://github.com/eXascaleInfolab/LFR-Benchmark_UndirWeightOvp) algorithms.


### Requirements:

Notebooks are created in Julia 1.7.0 with following packages:
- DataFrames v1.3.0
- Graphs v1.4.1
- IJulia v1.23.2
- PyPlot v2.10.0
- SmoothingSplines v0.3.1
- StatsBase v0.33.13

### Datasets

Three datasets are provided with notebooks:
- [coms_and_degs](https://drive.google.com/file/d/1u1YibNJRwp1LBCL0KEpCY1qqFdM8sP-v/view?usp=sharing) - community and degree distribution files used in the experiments.
- [speed_test](https://drive.google.com/file/d/1MfPL8JOWObJmlvdErOCGsd65zENlJqwe/view?usp=sharing) - data and logfiles for speed benchmark.
- [properties_test](https://drive.google.com/file/d/1IIFcSdkmf4GVAbouW1BPq-fCJ-BJO5DK/view?usp=sharing) - data and logfiles for properties testing.
- [graphs](https://drive.google.com/file/d/1wekpTdvgsKwSuPhNstXl-c9zIM3uLsdK/view?usp=sharing) - dataset contains sample of graphs generated with ABCD, LFR and LFR [NetworKit](https://networkit.github.io/) implementation. Used to generate the $knn(\ell)$ function.