#!/usr/bin/env python3
"""
Regenerate file_mapping.json after filtering datasets
This ensures the R pipeline can find all the correct files
"""

import json
from pathlib import Path
import os

# Configuration
FILTERED_DIR = Path("../../02_filtered_references")

def check_dataset_files(dataset_name):
    """Check which files exist for a dataset and return their paths"""
    dataset_dir = FILTERED_DIR / dataset_name
    
    if not dataset_dir.exists():
        print(f"❌ Directory does not exist: {dataset_dir}")
        return None
    
    print(f"\n✓ Checking {dataset_name}")
    
    # Check for each sex type
    sex_types = ["female", "male", "mixed"]
    files_found = {}
    
    for sex in sex_types:
        files_found[sex] = {}
        
        # Check metadata file - try multiple patterns
        metadata_patterns = [
            f"metadata_{sex}.csv",
            f"metadata_{sex}_harmonized.csv",
            f"metadata_HIP_sex{sex}.csv"
        ]
        
        metadata_found = False
        for mf in metadata_patterns:
            if (dataset_dir / mf).exists():
                files_found[sex]['metadata'] = str(dataset_dir / mf)
                metadata_found = True
                print(f"  ✓ {sex}: Found metadata - {mf}")
                break
        
        if not metadata_found:
            print(f"  ❌ {sex}: No metadata file found")
            
        # Check expression matrix files - try multiple patterns
        expr_patterns = [
            f"expression_matrix_{sex}.h5",
            f"matrix_HIP_sex{sex}_filtered.csv",
            f"expression_{sex}.h5"
        ]
        
        expr_found = False
        for ef in expr_patterns:
            if (dataset_dir / ef).exists():
                # Check file size to ensure it's not empty
                size_mb = os.path.getsize(dataset_dir / ef) / (1024 * 1024)
                if size_mb < 0.1:  # Less than 100KB is suspicious
                    print(f"  ⚠️  {sex}: Found {ef} but it's very small ({size_mb:.2f} MB)")
                else:
                    files_found[sex]['expression'] = str(dataset_dir / ef)
                    expr_found = True
                    print(f"  ✓ {sex}: Found expression - {ef} ({size_mb:.1f} MB)")
                break
        
        if not expr_found:
            print(f"  ❌ {sex}: No expression matrix found")
    
    return files_found

def main():
    print("="*80)
    print("REGENERATING FILE MAPPING FOR R PIPELINE")
    print("="*80)
    
    # Define datasets
    datasets = ["mouse10x_2020", "mouse_smartseq_2019", "yao_hippo_10x"]
    
    # Create mapping
    file_mapping = {}
    
    for dataset in datasets:
        result = check_dataset_files(dataset)
        if result:
            file_mapping[dataset] = result
    
    # Save mapping
    mapping_file = FILTERED_DIR / "file_mapping.json"
    with open(mapping_file, "w") as f:
        json.dump(file_mapping, f, indent=2)
    
    print(f"\n✓ File mapping saved to: {mapping_file}")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    complete_count = 0
    for dataset in datasets:
        if dataset in file_mapping:
            # Check if all sex types have both metadata and expression
            all_complete = True
            for sex in ["female", "male", "mixed"]:
                if sex not in file_mapping[dataset] or \
                   'metadata' not in file_mapping[dataset][sex] or \
                   'expression' not in file_mapping[dataset][sex]:
                    all_complete = False
                    break
            
            if all_complete:
                print(f"✓ {dataset}: Complete")
                complete_count += 1
            else:
                print(f"⚠ {dataset}: Partially complete")
                # Show what's missing
                for sex in ["female", "male", "mixed"]:
                    if sex not in file_mapping[dataset]:
                        print(f"    Missing: {sex}")
                    elif 'metadata' not in file_mapping[dataset][sex]:
                        print(f"    Missing: {sex} metadata")
                    elif 'expression' not in file_mapping[dataset][sex]:
                        print(f"    Missing: {sex} expression")
        else:
            print(f"❌ {dataset}: Not found")
    
    print(f"\nTotal complete: {complete_count}/{len(datasets)}")
    
    # Show the mapping
    print("\nFile mapping content:")
    print(json.dumps(file_mapping, indent=2))

if __name__ == "__main__":
    main()