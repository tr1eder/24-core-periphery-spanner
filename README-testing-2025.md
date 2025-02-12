To run a single small graph with one of the four spanners

```bash
bash_convert-to-gbbs.sh -filespecifier Slashdot

bazel run benchmarks/Spanner/FGV_Baseline:Spanner_main -- -s -src 10 ~/spanner/gbbs-com-dblp.txt
```
