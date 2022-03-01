
# CHANGES

1. `-ax map-pb` to `-ax map-hifi`

2. `filter_SV.py` can't find file

3. `combined` missing headers

4. `svnn` call tools using relative path. I removed this so future users need to add all tools to the env vars. This include all tools in the `bin` folder and all python scripts in the `source` folder.

5. `find_SR` has a problem reading sam due to different encodings used by PacBio read quality. To fix it, do `reformat.sh in=test.fq out=fixed.fq qin=33 qout=64 maxcalledquality=41`

# Get Started

```
# install
module load conda3/202011
conda create -n long_reads -c bioconda minimap2
source activate long_reads
conda install -c bioconda ngmlr
conda install -c bioconda samtools=1.9
conda install -c bioconda bamtools
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
conda install -c bioconda vulcan
conda install --channel bioconda svim
conda install sniffles
pip install -e .
conda install pandas
conda install -c intel scikit-learn
conda install seaborn

```

```
# test run
svnn -r ~/Data/Human/hg38/fasta/hg38.main.fa -q fixed.fq -s1 1 -s2 1
```







