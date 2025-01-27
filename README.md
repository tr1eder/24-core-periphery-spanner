# spanner
AAPSP - finding efficient and provably well-approximating spanner using a core-periphery decomposition

### License
This repository is licensed under the MIT License. However, the following subfolders contain code that is independently licensed:

- **GBBS** (Graph Based Benchmark Suite) - [Original Repo URL](https://github.com/paralg/gbbs/)
- **ParlayLib** (A Toolkit for Programming Parallel Algorithms on Shared-Memory Multicore Machines) - [Original Repo URL](https://github.com/cmuparlay/parlaylib)

Please refer to the respective `LICENSE` files for details.

### Install BAZEL - on EULER
```powershell
curl -fsSL https://bazel.build/bazel-release.pub.gpg | gpg --dearmor > ~/bazel-archive-keyring.gpg
curl -LO https://github.com/bazelbuild/bazel/releases/download/7.3.1/bazel-7.3.1-installer-linux-x86_64.sh
chmod +x bazel-7.3.1-installer-linux-x86_64.sh
./bazel-7.3.1-installer-linux-x86_64.sh --user
export PATH="$HOME/bin:$PATH"
bazel --version
```

### Install GCC - on EULER (tested for gcc 11.4. and 12.2.)
```powershell
module load stack/.2024-06-silent
module load gcc/12.2.0
gcc --version
```

### Copy spanner directory without the 4 symlinks - on LAPTOP
```powershell
scp -r ./spanner/ timrieder@euler.ethz.ch:~/
```

### Configure & run bash scripts - on EULER
```powershell
cd ~/spanner
chmod +x bash_convert-to-gbbs.sh
chmod +x bash_run-spanners.sh
chmod +x slurm_convert-to-gbbs.sh
chmod +x slurm_run-spanners.sh
dos2unix slurm_convert-to-gbbs.sh
dos2unix slurm_single-convert-to-gbbs.slurm
dos2unix slurm_run-spanners.sh
dos2unix slurm_single-run-spanners.slurm
./bash_convert-to-gbbs.sh
./bash_run-spanners.sh
./slurm_convert-to-gbbs.sh
./slurm_run-spanners.sh
```

### Run - a test
```powershell
cd ~/spanner/gbbs
bazel run benchmarks/Spanner/FGV_Baseline:Spanner_main -- -s ~/spanner/graphs-sanitized-gbbs/soc-hamsterster.txt
```

### RUN - convert from snap to gbbs format (bash_convert-to-gbbs.sh)
_run this file to convert all graphs directly (synchronous)_
```powershell
#!/bin/bash
cd ~/spanner/gbbs
snapfolder=~/spanner/graphs-sanitized-snap
gbbsfolder=~/spanner/graphs-sanitized-gbbs
for file in "$snapfolder"/*; do     
	filename=$(basename "$file");  
	sed -i '1s/^/# /' "$file";
	bazel run //utils:snap_converter -- -s -i "$file" -o "$gbbsfolder/$filename"; 
done
```

### RUN - run spanner algs (bash_run-spanners.sh)
_run this file to run all graph benchmarks directly (synchronous)_
```powershell
#!/bin/bash
cd ~/spanner/gbbs
snapfolder=~/spanner/graphs-sanitized-snap
gbbsfolder=~/spanner/graphs-sanitized-gbbs
for file in "$gbbsfolder"/*; do
	bazel run benchmarks/Spanner/FGV_Baseline:Spanner_main -- -s "$file"
	bazel run benchmarks/Spanner/MPVX_Baseline:Spanner_main -- -s "$file"
	bazel run benchmarks/Spanner/MPVX_CompactSpanner:Spanner_main -- -s "$file"
done
```

### RUN (slurm_convert-to-gbbs.sh)
_run this file to convert all graphs using slurm (asynchronous & in parallel)_
```powershell
#!/bin/bash
snapfolder=~/spanner/graphs-sanitized-snap
for file in "$snapfolder"/*; do
	sbatch slurm_single-convert-to-gbbs.slurm "$file"
done
```

### RUN (slurm_single-convert-to-gbbs.slurm)
_helper to convert a single graph_
```powershell
#!/bin/bash
#SBATCH --job-name=convert-to-gbbs_%j.slurm                 # Job name
#SBATCH --output=slurmlog/convert-to-gbbs_%j.log            # Output file
#SBATCH --error=slurmerr/convert-to-gbbs_%j.err             # Error file
#SBATCH --ntasks=1                                          # Number of tasks (usually 1 for single commands)
#SBATCH --time=01:00:00                                     # Time limit (e.g., 1 hour)
#SBATCH --mem-per-cpu=8G                                    # Memory requirement
cd ~/spanner/gbbs
snapfolder=~/spanner/graphs-sanitized-snap
gbbsfolder=~/spanner/graphs-sanitized-gbbs
filename=$(basename $1)                                     # Get filename from SLURM_ARRAY_TASK_ID
sed -i '1s/^/# /' "$snapfolder/$filename"
bazel --output_base=~/spanner/tmp/bazel_out_$SLURM_JOB_ID run //utils:snap_converter -- -s -i "$snapfolder/$filename" -o "$gbbsfolder/$filename"
```

### RUN (slurm_run-spanners.sh)
_run this file to run all graph benchmarks using slurm (asynchronous & in parallel)_
```powershell
#!/bin/bash
gbbsfolder=~/spanner/graphs-sanitized-gbbs
for file in "$gbbsfolder"/*; do
	sbatch slurm_single-run-spanners.slurm "$file" "FGV_Baseline"
	sbatch slurm_single-run-spanners.slurm "$file" "MPVX_Baseline"
	sbatch slurm_single-run-spanners.slurm "$file" "MPVX_CompactSpanner"
done
``` 

### RUN (slurm_single-run-spanners.slurm)
_helper to run run a single graph through a single benchmark_
```powershell
#!/bin/bash
#SBATCH --job-name=run-spanners_%j.slurm       # Job name
#SBATCH --output=slurmlog/run-spanners_%j.log  # Output file
#SBATCH --error=slurmerr/run-spanners_%j.err   # Error file
#SBATCH --ntasks=1                             # Number of tasks (usually 1 for single commands)
#SBATCH --time=01:00:00                        # Time limit (e.g., 1 hour)
#SBATCH --mem-per-cpu=8G                       # Memory requirement
cd ~/spanner/gbbs
file=$1
spannerversion=$2
bazel --output_base=~/spanner/tmp/bazel_out_$SLURM_JOB_ID run benchmarks/Spanner/$spannerversion:Spanner_main -- -s "$file"
```

### COPY RESULTS BACK
```
scp -r timrieder@euler.ethz.ch:~/spanner/slurmlog/ ./spanner/graphs-results-timing
scp -r timrieder@euler.ethz.ch:~/spanner/graphs-results/orkut_spanner_MPVX* ./spanner/graphs-results
```

### EXTRACT TIMING INFORMATION (extractTimingInformation.py)
