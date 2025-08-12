#!/bin/bash

# Environment Health Check Script
# Comprehensive check of all metagenomic analysis environments
# Author: Generated for farting_trees project

set -e

echo "=============================================="
echo "🔍 Environment Health Check"
echo "=============================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to check if conda environment exists
check_env_exists() {
    local env_name=$1
    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        return 1
    fi
}

# Function to check if a command exists in an environment
check_command_in_env() {
    local env_name=$1
    local command=$2
    local description=$3
    
    # Activate environment and check command
    if conda run -n "$env_name" which "$command" >/dev/null 2>&1; then
        printf "    ✅ %-20s %s\n" "$command" "$description"
        return 0
    else
        printf "    ❌ %-20s %s (NOT FOUND)\n" "$command" "$description"
        return 1
    fi
}

# Function to check environment health
check_environment() {
    local env_name=$1
    local description=$2
    shift 2
    local commands=("$@")
    
    echo -e "${BLUE}📦 Environment: ${env_name}${NC} (${description})"
    
    if ! check_env_exists "$env_name"; then
        echo -e "  ${RED}❌ Environment not found${NC}"
        echo ""
        return 1
    fi
    
    local failed_commands=0
    local total_commands=${#commands[@]}
    
    for ((i=0; i<${#commands[@]}; i+=2)); do
        local command="${commands[i]}"
        local desc="${commands[i+1]}"
        
        if ! check_command_in_env "$env_name" "$command" "$desc"; then
            ((failed_commands++))
        fi
    done
    
    local success_rate=$(( (total_commands/2 - failed_commands) * 100 / (total_commands/2) ))
    
    if [ $failed_commands -eq 0 ]; then
        echo -e "  ${GREEN}🟢 Status: PERFECT (${success_rate}%)${NC}"
    elif [ $failed_commands -le 1 ]; then
        echo -e "  ${YELLOW}🟡 Status: GOOD (${success_rate}%)${NC}"
    else
        echo -e "  ${RED}🔴 Status: NEEDS ATTENTION (${success_rate}%)${NC}"
    fi
    
    echo ""
    return $failed_commands
}

# Start health check
echo "Checking conda environments for metagenomic analysis pipeline..."
echo ""

total_issues=0

# Check each environment
check_environment "metagenome_assembly" "Assembly & Host Removal" \
    "bwa" "BWA aligner" \
    "samtools" "SAMtools"
total_issues=$((total_issues + $?))

check_environment "metagenome_binning" "Genome Binning" \
    "metabat2" "MetaBAT2 binner" \
    "vamb" "VAMB binner" \
    "jgi_summarize_bam_contig_depths" "JGI depth calculator"
total_issues=$((total_issues + $?))

check_environment "maxbin2_env" "MaxBin2 Binning" \
    "MaxBin" "MaxBin2 binary" \
    "run_MaxBin.pl" "MaxBin2 Perl wrapper"
total_issues=$((total_issues + $?))

check_environment "dastool_env" "DAS Tool Integration" \
    "DAS_Tool" "DAS Tool" \
    "R" "R statistical software"
total_issues=$((total_issues + $?))

check_environment "checkm_env" "Genome Quality Assessment" \
    "checkm" "CheckM quality tool" \
    "prodigal" "Prodigal gene caller"
total_issues=$((total_issues + $?))

# Summary
echo "=============================================="
echo "📊 Health Check Summary"
echo "=============================================="

if [ $total_issues -eq 0 ]; then
    echo -e "${GREEN}🎉 All environments are HEALTHY!${NC}"
    echo "✅ Your metagenomic binning pipeline is ready to use"
    echo ""
    echo "🚀 Next steps:"
    echo "  1. Run quality control: ./code/01_quality_control.sh"
    echo "  2. Run assembly: ./code/03_assembly.sh"
    echo "  3. Run binning: ./code/04_binning.sh"
elif [ $total_issues -le 2 ]; then
    echo -e "${YELLOW}⚠️  Minor issues detected (${total_issues} problems)${NC}"
    echo "🔧 Most tools are working, some manual fixes may be needed"
    echo ""
    echo "💡 Suggested actions:"
    echo "  - Check individual package installations"
    echo "  - Consider reinstalling problematic packages"
elif [ $total_issues -le 5 ]; then
    echo -e "${YELLOW}🟡 Several issues detected (${total_issues} problems)${NC}"
    echo "🔨 Multiple environments need attention"
    echo ""
    echo "💡 Suggested actions:"
    echo "  - Run environment setup again: ./code/00_environment_setup.sh"
    echo "  - Check conda channel priorities"
else
    echo -e "${RED}🚨 Major issues detected (${total_issues} problems)${NC}"
    echo "🛠️  Significant environment problems found"
    echo ""
    echo "💡 Suggested actions:"
    echo "  - Recreate environments: ./code/00_environment_setup.sh"
    echo "  - Check conda installation"
    echo "  - Verify internet connectivity for package downloads"
fi

echo ""
echo "📋 Environment Status Legend:"
echo "  🟢 PERFECT (100%) - All tools working"
echo "  🟡 GOOD (80-99%) - Minor issues, mostly functional"  
echo "  🔴 NEEDS ATTENTION (<80%) - Significant problems"
echo ""

# Provide environment-specific guidance
if [ $total_issues -gt 0 ]; then
    echo "🔧 Quick fixes available:"
    echo "  - MaxBin2 issues: Run ./code/fix_maxbin2.sh"
    echo "  - Package missing: conda install -n ENV_NAME -c bioconda PACKAGE"
    echo "  - Environment recreation: conda env remove -n ENV_NAME && ./code/00_environment_setup.sh"
    echo ""
fi

echo "For detailed troubleshooting, see: TROUBLESHOOTING.md"
