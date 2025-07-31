# Host Removal for Tree Microbiome Analysis

## Why Host Removal is Important

When studying tree microbiomes, your sequencing data will contain:
1. **Microbial DNA** (bacteria, archaea, fungi) - what we want to study
2. **Tree/plant host DNA** - what we want to remove
3. **Other contaminants** - environmental DNA, human DNA, etc.

Host removal is crucial because:
- **Tree DNA can be 10-90%** of total sequencing reads
- **Reduces computational burden** for assembly and binning
- **Improves assembly quality** by removing host contamination
- **Focuses analysis** on the microbial community

## Host Removal Strategy

### Step 1: Reference Genome Selection
We use multiple plant reference genomes to catch different types of host contamination:

1. **Arabidopsis thaliana** - Model plant, catches general plant sequences
2. **Tree-specific genomes** - Based on your study system:
   - **Populus trichocarpa** - Poplar/aspen trees
   - **Quercus robur** - Oak trees  
   - **Pinus taeda** - Pine trees
   - **Picea abies** - Spruce trees

### Step 2: Mapping and Filtering
1. **Map reads** to combined host genomes using BWA
2. **Extract unmapped reads** (non-host sequences)
3. **Generate statistics** on host contamination levels

### Step 3: Quality Control
- Monitor **host removal efficiency**
- Check for **over-aggressive filtering**
- Ensure **microbial reads are retained**

## Configuration

Edit `config/host_removal_config.sh` to specify which tree genomes to use:

```bash
# Enable/disable specific tree genomes
DOWNLOAD_ARABIDOPSIS=true    # Always recommended
DOWNLOAD_POPULUS=true        # If studying poplar/aspen
DOWNLOAD_OAK=true           # If studying oak
DOWNLOAD_PINE=false         # If studying pine
DOWNLOAD_SPRUCE=false       # If studying spruce
```

## Expected Results

### Typical Host Contamination Levels:
- **Leaf samples**: 30-70% host DNA
- **Bark samples**: 10-40% host DNA  
- **Root samples**: 20-50% host DNA
- **Soil samples**: 5-20% host DNA

### Output Files:
- `clean_reads/` - Host-removed FASTQ files
- `host_removal_summary.txt` - Summary statistics
- `stats/` - Per-sample detailed statistics

## Quality Assessment

Good host removal should show:
- **Reasonable host percentage** (see ranges above)
- **Majority of reads retained** (>50% for most samples)
- **Consistent results** across similar sample types

### Red Flags:
- **>95% reads removed** - may be over-aggressive
- **<5% host DNA** in leaf samples - may indicate poor host genome match
- **Highly variable** removal rates between similar samples

## Troubleshooting

### High Host Contamination (>80% removed)
- Check if you're using the **correct host genome**
- Consider using **multiple tree genomes**
- Verify **sample collection methods**

### Low Host Contamination (<5% removed)  
- May indicate **good sample preparation**
- Or **poor host genome match**
- Check if using **appropriate reference genome**

### Assembly Issues After Host Removal
- If assembly quality decreases, check if removing too aggressively
- Consider adjusting **mapping stringency**
- Verify **microbial content** in samples

## Integration with Pipeline

The host removal step is integrated between QC and assembly:

```
Raw reads → QC → Host Removal → Assembly → Binning
```

The assembly script automatically detects and uses host-removed reads if available.

## Advanced Options

### Custom Tree Genomes
If studying a specific tree species not included:

1. Find genome at NCBI: https://www.ncbi.nlm.nih.gov/genome/
2. Add URL to `config/host_removal_config.sh`
3. Enable download in configuration

### Alternative Approaches
- **Kraken2/Bracken** - Taxonomic classification-based removal
- **BMTagger** - Alternative host screening tool
- **Deconseq** - Contamination detection tool

## References

- **BWA**: Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997
- **Tree genomes**: Various genome consortiums (JGI, NCBI)
- **Best practices**: Eisenhofer et al. (2019) "Contamination in Low Microbial Biomass Microbiome Studies"
