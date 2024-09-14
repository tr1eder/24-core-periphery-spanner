# spanner - 5
AAPSP - finding efficient and provably well-approximating spanner using a core-periphery decomposition


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
```powershell
#!/bin/bash
snapfolder=~/spanner/graphs-sanitized-snap
for file in "$snapfolder"/*; do
	sbatch slurm_single-convert-to-gbbs.slurm "$file"
done
```

### RUN (slurm_single-convert-to-gbbs.slurm)
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
