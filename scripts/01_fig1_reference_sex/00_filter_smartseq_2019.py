#!/usr/bin/env python3
"""
CORRECTED filtering script for Mouse Whole Cortex & Hippocampus — SMART-seq (2019)
Properly identifies hippocampal cells using cluster_label patterns
"""

import pandas as pd
import numpy as np
import h5py
from scipy.sparse import csc_matrix
from pathlib import Path
import time
import json

# === CONFIGURATION ===
DATASET_NAME = "mouse_smartseq_2019"
METADATA_PATH = r"D:/1Reference Datasets/Mouse Whole Cortex & Hippocampus — SMART-seq (2019)/metadata.csv"
EXPRESSION_PATH = r"D:/1Reference Datasets/Mouse Whole Cortex & Hippocampus — SMART-seq (2019)/expression_matrix.hdf5"
OUTPUT_DIR = Path(f"../../02_filtered_references/{DATASET_NAME}")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === HIPPOCAMPUS IDENTIFICATION ===
# Since region_label != "HIP" for this dataset, use cluster_label patterns
HIPPO_CLUSTER_PATTERNS = [
    # Main hippocampal regions
    r"\d+_CA1", r"\d+_CA2", r"\d+_CA3", r"\d+_DG",
    # Specific patterns from Allen Brain Atlas
    r"CA1-ProS", r"CA2-IG-FC", r"CA3",
    # Dentate gyrus patterns
    r"DG-", r"\d+_DG",
    # Subiculum patterns
    r"\d+_SUB", r"SUB-ProS", r"CT SUB", r"NP SUB",
    # Additional hippocampal interneurons
    r"Sst Chodl"  # Hippocampus-enriched interneuron type
]

# Expected hippocampal subclass labels
HIPPO_SUBCLASS_LABELS = [
    "CA1-ProS", "CA2-IG-FC", "CA3", "DG", 
    "SUB-ProS", "CT SUB", "NP SUB"
]

SEX_LABELS = {"female": "F", "male": "M", "mixed": None}

def identify_hippocampal_cells(metadata):
    """Identify hippocampal cells using cluster patterns since region_label != HIP"""
    
    print("\nIdentifying hippocampal cells...")
    
    # Initialize mask
    hippo_mask = pd.Series(False, index=metadata.index)
    
    # Method 1: Check cluster_label for hippocampal patterns
    if "cluster_label" in metadata.columns:
        print("  Checking cluster_label patterns...")
        for pattern in HIPPO_CLUSTER_PATTERNS:
            pattern_mask = metadata["cluster_label"].astype(str).str.contains(
                pattern, case=False, na=False, regex=True
            )
            n_found = pattern_mask.sum()
            if n_found > 0:
                print(f"    Pattern '{pattern}': {n_found} cells")
                hippo_mask |= pattern_mask
    
    # Method 2: Check subclass_label
    if "subclass_label" in metadata.columns:
        print("  Checking subclass_label...")
        subclass_mask = metadata["subclass_label"].isin(HIPPO_SUBCLASS_LABELS)
        n_subclass = subclass_mask.sum()
        if n_subclass > 0:
            print(f"    Subclass labels: {n_subclass} cells")
            hippo_mask |= subclass_mask
    
    # Method 3: For cells in SUB-ProS region, check if they're truly hippocampal
    if "region_label" in metadata.columns:
        # SUB-ProS is hippocampal, but ACA is not
        subiculum_mask = metadata["region_label"] == "SUB-ProS"
        n_sub = subiculum_mask.sum()
        if n_sub > 0:
            print(f"    SUB-ProS region: {n_sub} cells")
            hippo_mask |= subiculum_mask
    
    total_hippo = hippo_mask.sum()
    print(f"\nTotal hippocampal cells identified: {total_hippo}")
    
    # Validate by showing cell type distribution
    if total_hippo > 0 and "cluster_label" in metadata.columns:
        hippo_clusters = metadata.loc[hippo_mask, "cluster_label"].value_counts().head(10)
        print("\nTop 10 hippocampal cluster types:")
        for cluster, count in hippo_clusters.items():
            print(f"  {cluster}: {count}")
    
    return hippo_mask

def load_smartseq_sparse_matrix(h5_path, sample_indices):
    """Load SMART-seq sparse matrix data"""
    
    print("\nLoading SMART-seq expression matrix...")
    
    with h5py.File(h5_path, "r") as f:
        # SMART-seq uses sparse format with exon + intron counts
        print("  Reading sparse matrix components...")
        
        # Get matrix components
        i = f["data/exon/i"][:]  # Row indices (genes)
        p = f["data/exon/p"][:]  # Column pointers (cells)
        x = f["data/exon/x"][:]  # Expression values
        
        # Get dimensions
        dims = f["data/exon/dims"][:]
        n_cells, n_genes = int(dims[0]), int(dims[1])
        
        print(f"  Full matrix: {n_cells} cells x {n_genes} genes")
        
        # Create full sparse matrix
        full_matrix = csc_matrix((x, i, p), shape=(n_cells, n_genes))
        
        # Get gene names
        gene_names = [g.decode('utf-8') for g in f["gene_names"][:n_genes]]
        
        # Get all sample names
        all_sample_names = [s.decode('utf-8') for s in f["sample_names"][:n_cells]]
        
    # Create sample lookup
    sample_lookup = {name: idx for idx, name in enumerate(all_sample_names)}
    
    # Map our samples to indices
    valid_indices = []
    valid_samples = []
    missing_samples = []
    
    for idx, sample in enumerate(sample_indices):
        if sample in sample_lookup:
            valid_indices.append(sample_lookup[sample])
            valid_samples.append(sample)
        else:
            missing_samples.append(sample)
    
    if missing_samples:
        print(f"  WARNING: {len(missing_samples)} samples not found in expression matrix")
    
    print(f"  Extracting data for {len(valid_indices)} cells...")
    
    # Extract subset (cells are rows, genes are columns)
    subset_matrix = full_matrix[valid_indices, :].toarray()
    
    # Transpose to genes x cells
    subset_matrix = subset_matrix.T
    
    print(f"  Subset matrix: {subset_matrix.shape[0]} genes x {subset_matrix.shape[1]} cells")
    
    return subset_matrix, gene_names, valid_samples

def save_filtered_data(sex_key, sex_label):
    """Filter and save data for specific sex"""
    
    print(f"\n{'='*60}")
    print(f"Processing {sex_key.upper()} reference")
    print('='*60)
    
    # Load metadata
    print("Loading metadata...")
    metadata = pd.read_csv(METADATA_PATH, low_memory=False)
    print(f"Total cells in dataset: {len(metadata)}")
    
    # Show region distribution
    if "region_label" in metadata.columns:
        print("\nRegion distribution:")
        print(metadata["region_label"].value_counts())
    
    # Identify hippocampal cells
    hippo_mask = identify_hippocampal_cells(metadata)
    
    if hippo_mask.sum() == 0:
        print("ERROR: No hippocampal cells found!")
        return None, None
    
    # Apply sex filter
    if sex_label:
        sex_mask = metadata["donor_sex_label"] == sex_label
        final_mask = hippo_mask & sex_mask
        print(f"\nFiltering for sex='{sex_label}'")
        print(f"  Sex-specific cells: {sex_mask.sum()}")
        print(f"  Hippocampal + sex filter: {final_mask.sum()}")
    else:
        final_mask = hippo_mask
        print(f"\nKeeping all sexes (mixed reference)")
    
    if final_mask.sum() == 0:
        print(f"WARNING: No cells remaining after filtering for {sex_key}")
        return None, None
    
    # Get filtered metadata
    filtered_metadata = metadata[final_mask].copy()
    print(f"\nFinal cells after filtering: {len(filtered_metadata)}")
    
    # Get sample names for expression matrix
    sample_names = filtered_metadata["sample_name"].tolist()
    
    # Load expression data
    expression_matrix, gene_names, valid_samples = load_smartseq_sparse_matrix(
        EXPRESSION_PATH, sample_names
    )
    
    # Update metadata to match valid samples
    filtered_metadata = filtered_metadata[filtered_metadata["sample_name"].isin(valid_samples)]
    
    # Save outputs
    print("\nSaving filtered data...")
    
    # Save metadata
    metadata_output = OUTPUT_DIR / f"metadata_{sex_key}.csv"
    filtered_metadata.to_csv(metadata_output, index=False)
    print(f"  Metadata saved: {metadata_output}")
    
    # Create expression DataFrame
    expr_df = pd.DataFrame(
        expression_matrix,
        index=gene_names,
        columns=valid_samples
    )
    
    # Save as HDF5 for consistency with mouse10x format
    h5_output = OUTPUT_DIR / f"expression_matrix_{sex_key}.h5"
    with h5py.File(h5_output, 'w') as f:
        # Save expression data
        f.create_dataset('data/expression', data=expression_matrix, compression='gzip')
        
        # Save gene names
        gene_names_encoded = [g.encode('utf-8') for g in gene_names]
        f.create_dataset('data/gene_names', data=gene_names_encoded)
        
        # Save cell names
        cell_names_encoded = [c.encode('utf-8') for c in valid_samples]
        f.create_dataset('data/cell_names', data=cell_names_encoded)
        
        # Add metadata attributes
        f.attrs['n_genes'] = len(gene_names)
        f.attrs['n_cells'] = len(valid_samples)
        f.attrs['dataset'] = DATASET_NAME
        f.attrs['sex_filter'] = sex_key
    
    print(f"  Expression matrix saved: {h5_output}")
    
    # Also save as CSV for compatibility
    csv_output = OUTPUT_DIR / f"expression_matrix_{sex_key}.csv"
    expr_df.to_csv(csv_output)
    print(f"  Expression matrix (CSV) saved: {csv_output}")
    
    # Save summary
    summary = {
        "dataset": DATASET_NAME,
        "sex_filter": sex_key,
        "total_cells_original": len(metadata),
        "hippocampal_cells_identified": int(hippo_mask.sum()),
        "cells_after_sex_filter": len(filtered_metadata),
        "cells_with_expression": len(valid_samples),
        "genes": len(gene_names),
        "hippocampal_identification_method": "cluster_label patterns + subclass_label + SUB-ProS region",
        "validation": {
            "region_distribution": filtered_metadata["region_label"].value_counts().to_dict() if "region_label" in filtered_metadata else {},
            "top_clusters": filtered_metadata["cluster_label"].value_counts().head(20).to_dict() if "cluster_label" in filtered_metadata else {},
            "sex_distribution": filtered_metadata["donor_sex_label"].value_counts().to_dict() if "donor_sex_label" in filtered_metadata else {}
        }
    }
    
    summary_output = OUTPUT_DIR / f"filtering_summary_{sex_key}.json"
    with open(summary_output, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"  Summary saved: {summary_output}")
    
    return filtered_metadata, expression_matrix

def main():
    """Main execution"""
    print("="*80)
    print("CORRECTED HIPPOCAMPAL FILTERING FOR MOUSE SMART-SEQ 2019")
    print("="*80)
    
    overall_start = time.time()
    
    # Process each sex group
    results = {}
    for sex_key, sex_label in SEX_LABELS.items():
        try:
            metadata, expression = save_filtered_data(sex_key, sex_label)
            if metadata is not None:
                results[sex_key] = {
                    "success": True,
                    "n_cells": len(metadata),
                    "n_genes": expression.shape[0] if expression is not None else 0
                }
            else:
                results[sex_key] = {
                    "success": False,
                    "error": "No cells found after filtering"
                }
        except Exception as e:
            print(f"\nERROR processing {sex_key}: {str(e)}")
            import traceback
            traceback.print_exc()
            results[sex_key] = {
                "success": False,
                "error": str(e)
            }
    
    # Summary
    print("\n" + "="*80)
    print("PROCESSING COMPLETE")
    print("="*80)
    print(f"Total time: {time.time() - overall_start:.1f} seconds")
    print("\nResults:")
    for sex_key, result in results.items():
        if result["success"]:
            print(f"  {sex_key}: ✓ {result['n_cells']} cells, {result['n_genes']} genes")
        else:
            print(f"  {sex_key}: ✗ Failed - {result['error']}")

if __name__ == "__main__":
    main()