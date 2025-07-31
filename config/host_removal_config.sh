# Host Removal Configuration for Tree Microbiomes
# Edit this file to customize host removal for your specific tree species

# Common tree and plant genomes for host removal
# Uncomment and modify based on your study system

# Model plant (always included for general plant contamination)
ARABIDOPSIS_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz"

# Deciduous trees
POPULUS_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/775/GCF_000002775.5_Pop_tri_v4.0/GCF_000002775.5_Pop_tri_v4.0_genomic.fna.gz"  # Poplar
OAK_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/294/415/GCF_932294415.1_dhQue1.1/GCF_932294415.1_dhQue1.1_genomic.fna.gz"  # Oak

# Evergreen/coniferous trees  
PINE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/404/065/GCF_000404065.1_Ptaeda2.0/GCF_000404065.1_Ptaeda2.0_genomic.fna.gz"  # Loblolly Pine
SPRUCE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/411/955/GCF_000411955.1_Pabies1.0/GCF_000411955.1_Pabies1.0_genomic.fna.gz"  # Norway Spruce

# Fruit trees
APPLE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/114/115/GCF_002114115.1_ASM211411v1/GCF_002114115.1_ASM211411v1_genomic.fna.gz"

# Set which genomes to download
# Change these to true/false based on your tree species
DOWNLOAD_ARABIDOPSIS=true
DOWNLOAD_POPULUS=false
DOWNLOAD_OAK=false  
DOWNLOAD_PINE=false
DOWNLOAD_SPRUCE=false
DOWNLOAD_APPLE=false

# Host removal parameters
BWA_THREADS=8
MIN_MAPPING_QUALITY=20

# If you know your specific tree species, you can add custom URLs here:
# CUSTOM_TREE_1_URL=""
# CUSTOM_TREE_1_NAME=""
