# Troubleshooting Conda Environment Issues

## The Problem You Encountered

The error occurs because `fastp=0.23.4` has very specific requirements for `libdeflate` that conflict with other packages in the environment. This is a common issue with bioinformatics tools that have strict dependency requirements.

## Quick Fixes (Try in Order)

### Fix 1: Use the Updated Script (Recommended)
The original script has been updated to use flexible versioning:
```bash
./scripts/01_setup_environment.sh
```

### Fix 2: Use Mamba Instead of Conda
Mamba is much better at resolving complex dependencies:
```bash
./scripts/01_setup_environment_mamba.sh
```

### Fix 3: Gradual Setup (Most Reliable)
If the above fail, use the safe setup script that creates separate environments:
```bash
./scripts/01_setup_environment_safe.sh
```

## Manual Troubleshooting Steps

### Step 1: Clean Up Failed Environment
```bash
# Remove any partially created environment
conda env remove -n metagenome_qc -y
```

### Step 2: Install Mamba (Better Dependency Resolver)
```bash
conda install mamba -n base -c conda-forge -y
```

### Step 3: Create Environment with Mamba
```bash
mamba create -n metagenome_qc -c bioconda -c conda-forge -y \
    fastqc multiqc fastp trimmomatic
```

### Step 4: Alternative - Use Latest Versions
```bash
# Create environment without version pinning
mamba create -n metagenome_qc -c bioconda -c conda-forge -y \
    fastqc multiqc fastp trimmomatic
```

### Step 5: Minimal Environment Approach
If all else fails, create minimal environments:
```bash
# QC environment with just essential tools
mamba create -n metagenome_qc -c bioconda -c conda-forge -y fastqc multiqc

# Separate trimming environment
mamba create -n trimming -c bioconda -c conda-forge -y fastp trimmomatic
```

## Why This Happens

1. **Strict Dependencies**: Bioinformatics tools often have very specific version requirements
2. **Library Conflicts**: Different tools may require incompatible versions of the same library
3. **Package Repository Issues**: Sometimes the conda package metadata has conflicts
4. **Platform Differences**: Linux/macOS/Windows may have different available packages

## Prevention Strategies

### Use Mamba Instead of Conda
```bash
# Install mamba once
conda install mamba -n base -c conda-forge -y

# Then use mamba for all environment creation
mamba create -n myenv -c bioconda -c conda-forge package1 package2
```

### Use Flexible Versioning
```bash
# Instead of: fastp=0.23.4
# Use: fastp (latest version)
# Or: fastp>=0.23 (minimum version)
```

### Separate Complex Tools
Create separate environments for tools with many dependencies:
- GTDB-Tk (many Python dependencies)
- CheckM (database requirements)
- VAMB (PyTorch dependencies)

## Testing Your Environment

After creating the QC environment, test it:
```bash
conda activate metagenome_qc
fastqc --version
fastp --version
multiqc --version
trimmomatic -version
```

## Alternative Tools

If fastp continues to cause issues, you can use alternatives:

### Replace fastp with Trimmomatic
```bash
# In your QC script, replace fastp with:
trimmomatic PE -threads ${THREADS} \
    ${DATA_DIR}/${sample}_R1.fastq.gz \
    ${DATA_DIR}/${sample}_R2.fastq.gz \
    ${RESULTS_DIR}/trimmed/${sample}_R1_trimmed.fastq.gz \
    ${RESULTS_DIR}/trimmed/${sample}_R1_unpaired.fastq.gz \
    ${RESULTS_DIR}/trimmed/${sample}_R2_trimmed.fastq.gz \
    ${RESULTS_DIR}/trimmed/${sample}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```

### Use cutadapt instead
```bash
mamba install -c bioconda cutadapt
```

## Docker Alternative

If conda continues to cause issues, consider using Docker:
```bash
# Pull pre-built image with all tools
docker pull quay.io/biocontainers/fastp:0.23.4--h5f740d0_0
```

## Next Steps

1. Try the updated setup script first
2. If that fails, use the mamba version
3. If still having issues, use the safe gradual setup
4. For persistent problems, consider the minimal environment approach

The pipeline is designed to be robust - even if some tools fail to install, you can run individual steps with available tools.
