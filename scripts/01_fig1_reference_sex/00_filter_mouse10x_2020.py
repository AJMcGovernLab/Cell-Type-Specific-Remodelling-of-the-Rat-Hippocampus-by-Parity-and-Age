#!/usr/bin/env python3
"""
CORRECTED filtering script for Mouse Whole Cortex & Hippocampus — 10x Genomics (2020)
Properly filters for true hippocampal cells based on region_label and cluster patterns
"""

import pandas as pd
import numpy as np
import h5py
from pathlib import Path
import time
import json

# === CONFIGURATION ===
DATASET_NAME = "mouse10x_2020"
METADATA_PATH = r"D:/1Reference Datasets/Mouse Whole Cortex & Hippocampus — 10x Genomics (2020)/metadata.csv"
EXPRESSION_PATH = r"D:/1Reference Datasets/Mouse Whole Cortex & Hippocampus — 10x Genomics (2020)/expression_matrix.hdf5"
OUTPUT_DIR = Path(f"../../02_filtered_references/{DATASET_NAME}")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === HIPPOCAMPUS IDENTIFICATION CRITERIA ===
# Primary method: region_label
HIPPO_REGION = "HIP"

# Secondary validation: cluster_label patterns for hippocampal subregions
HIPPO_CLUSTER_PATTERNS = [
    "CA1", "CA2", "CA3", "DG",  # Main hippocampal subfields
    "SUB", "ProS",               # Subiculum regions
    "CT SUB", "NP SUB"           # Specific subiculum types
]

# Subclass labels known to be hippocampal
HIPPO_SUBCLASS_LABELS = [
    "CA1-ProS", "CA2-IG-FC", "CA3", "DG", 
    "SUB-ProS", "CT SUB", "NP SUB"
]

SEX_LABELS = {"female": "F", "male": "M", "mixed": None}

def identify_hippocampal_cells(metadata):
    """
    Identify true hippocampal cells using multiple criteria
    """
    # Primary filter: region_label == "HIP"
    region_mask = metadata["region_label"] == HIPPO_REGION
    
    # If no HIP region found, try alternative approaches
    if region_mask.sum() == 0:
        print("WARNING: No cells found with region_label='HIP'")
        print("Attempting alternative identification methods...")
        
        # Method 2: Check cluster_label for hippocampal patterns
        cluster_mask = metadata["cluster_label"].astype(str).str.contains(
            '|'.join(HIPPO_CLUSTER_PATTERNS), 
            case=False, 
            na=False
        )
        
        # Method 3: Check subclass_label
        subclass_mask = metadata["subclass_label"].isin(HIPPO_SUBCLASS_LABELS)
        
        # Combine methods
        hippo_mask = cluster_mask | subclass_mask
        
        if hippo_mask.sum() > 0:
            print(f"Found {hippo_mask.sum()} hippocampal cells using cluster/subclass patterns")
    else:
        hippo_mask = region_mask
        print(f"Found {hippo_mask.sum()} cells in HIP region")
    
    return hippo_mask

def validate_hippocampal_cells(metadata, hippo_mask):
    """
    Validate and report on identified hippocampal cells
    """
    hippo_cells = metadata[hippo_mask]
    
    print("\n=== Hippocampal Cell Validation ===")
    print(f"Total hippocampal cells: {len(hippo_cells)}")
    
    # Show region distribution
    print("\nRegion distribution:")
    print(hippo_cells["region_label"].value_counts())
    
    # Show subclass distribution
    print("\nTop 20 subclass labels:")
    print(hippo_cells["subclass_label"].value_counts().head(20))
    
    # Show cluster patterns
    print("\nCluster labels containing hippocampal markers:")
    clusters = hippo_cells["cluster_label"].dropna().unique()
    hippo_clusters = [c for c in clusters if any(
        pattern in str(c) for pattern in HIPPO_CLUSTER_PATTERNS
    )]
    print(f"Found {len(hippo_clusters)} hippocampal cluster types")
    for cluster in sorted(hippo_clusters)[:20]:
        count = (hippo_cells["cluster_label"] == cluster).sum()
        print(f"  {cluster}: {count} cells")
    
    return hippo_cells

def read_expression_matrix_chunked(h5_file, cell_indices, chunk_size=5000):
    """
    Read expression matrix in chunks for memory efficiency
    """
    n_genes = h5_file["data/counts"].shape[0]
    n_cells = len(cell_indices)
    
    # Sort indices for efficient reading
    sorted_idx = np.argsort(cell_indices)
    sorted_indices = np.array(cell_indices)[sorted_idx]
    
    # Initialize output matrix
    expression_matrix = np.zeros((n_genes, n_cells), dtype=np.float32)
    
    print(f"Reading expression matrix ({n_genes} genes x {n_cells} cells)...")
    
    # Read in chunks
    for i in range(0, n_cells, chunk_size):
        end_idx = min(i + chunk_size, n_cells)
        chunk_indices = sorted_indices[i:end_idx]
        
        # Read chunk
        expression_matrix[:, i:end_idx] = h5_file["data/counts"][:, chunk_indices]
        
        if (i + chunk_size) % 10000 == 0:
            print(f"  Processed {i + chunk_size}/{n_cells} cells...")
    
    # Restore original order
    unsorted_matrix = np.zeros_like(expression_matrix)
    for i, orig_pos in enumerate(sorted_idx):
        unsorted_matrix[:, orig_pos] = expression_matrix[:, i]
    
    return unsorted_matrix

def save_filtered_data(sex_key, sex_label):
    """
    Filter and save data for specific sex
    """
    print(f"\n{'='*60}")
    print(f"Processing {sex_key.upper()} reference")
    print('='*60)
    
    # Load metadata
    print("Loading metadata...")
    metadata = pd.read_csv(METADATA_PATH, low_memory=False)
    print(f"Total cells in dataset: {len(metadata)}")
    
    # Identify hippocampal cells
    hippo_mask = identify_hippocampal_cells(metadata)
    
    # Apply sex filter if specified
    if sex_label:
        sex_mask = metadata["donor_sex_label"] == sex_label
        final_mask = hippo_mask & sex_mask
        print(f"\nFiltering for sex='{sex_label}'")
    else:
        final_mask = hippo_mask
        print("\nKeeping all sexes (mixed reference)")
    
    # Get filtered metadata
    filtered_metadata = metadata[final_mask]
    print(f"Cells after filtering: {len(filtered_metadata)}")
    
    # Validate cells
    validate_hippocampal_cells(metadata, final_mask)
    
    # Get cell IDs and indices
    cell_ids = filtered_metadata["sample_name"].tolist()
    
    # Load expression matrix
    print("\n=== Loading Expression Matrix ===")
    with h5py.File(EXPRESSION_PATH, "r") as f:
        # Get gene and cell names
        all_genes = [g.decode("utf-8") for g in f["data/gene"][:]]
        all_cells = [c.decode("utf-8") for c in f["data/samples"][:]]
        
        # Create lookup
        cell_lookup = {cid: idx for idx, cid in enumerate(all_cells)}
        
        # Get indices for our cells
        cell_indices = []
        missing_cells = []
        for cid in cell_ids:
            if cid in cell_lookup:
                cell_indices.append(cell_lookup[cid])
            else:
                missing_cells.append(cid)
        
        if missing_cells:
            print(f"WARNING: {len(missing_cells)} cells not found in expression matrix")
        
        print(f"Loading expression data for {len(cell_indices)} cells...")
        
        # Read expression data
        start_time = time.time()
        expression_data = read_expression_matrix_chunked(f, cell_indices)
        print(f"Expression matrix loaded in {time.time() - start_time:.1f} seconds")
    
    # Create output paths
    metadata_output = OUTPUT_DIR / f"metadata_{sex_key}.csv"
    expression_output = OUTPUT_DIR / f"expression_matrix_{sex_key}.h5"
    summary_output = OUTPUT_DIR / f"filtering_summary_{sex_key}.json"
    
    # Save filtered metadata
    filtered_metadata.to_csv(metadata_output, index=False)
    print(f"\nMetadata saved to: {metadata_output}")
    
    # Save expression matrix in HDF5 format
    print("Saving expression matrix...")
    with h5py.File(expression_output, "w") as f:
        # Create datasets
        f.create_dataset("data/expression", data=expression_data, compression="gzip")
        f.create_dataset("data/gene_names", data=[g.encode() for g in all_genes])
        f.create_dataset("data/cell_names", data=[all_cells[i].encode() for i in cell_indices])
        
        # Add metadata
        f.attrs["n_genes"] = len(all_genes)
        f.attrs["n_cells"] = len(cell_indices)
        f.attrs["dataset"] = DATASET_NAME
        f.attrs["sex_filter"] = sex_key
    
    print(f"Expression matrix saved to: {expression_output}")
    
    # Save summary
    summary = {
        "dataset": DATASET_NAME,
        "sex_filter": sex_key,
        "total_cells_original": len(metadata),
        "hippocampal_cells": int(hippo_mask.sum()),
        "cells_after_sex_filter": len(filtered_metadata),
        "genes": len(all_genes),
        "missing_cells": len(missing_cells),
        "filtering_criteria": {
            "primary": f"region_label == '{HIPPO_REGION}'",
            "secondary": "cluster_label patterns",
            "tertiary": "subclass_label matching"
        },
        "validation": {
            "region_distribution": filtered_metadata["region_label"].value_counts().to_dict(),
            "sex_distribution": filtered_metadata["donor_sex_label"].value_counts().to_dict() if "donor_sex_label" in filtered_metadata else {},
            "top_subclasses": filtered_metadata["subclass_label"].value_counts().head(10).to_dict()
        }
    }
    
    with open(summary_output, "w") as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary saved to: {summary_output}")
    
    return filtered_metadata, expression_data

def main():
    """
    Main execution function
    """
    print("="*80)
    print("CORRECTED HIPPOCAMPAL CELL FILTERING FOR MOUSE 10X 2020 DATASET")
    print("="*80)
    
    overall_start = time.time()
    
    # Process each sex group
    results = {}
    for sex_key, sex_label in SEX_LABELS.items():
        try:
            metadata, expression = save_filtered_data(sex_key, sex_label)
            results[sex_key] = {
                "success": True,
                "n_cells": len(metadata),
                "n_genes": expression.shape[0]
            }
        except Exception as e:
            print(f"\nERROR processing {sex_key}: {str(e)}")
            results[sex_key] = {
                "success": False,
                "error": str(e)
            }
    
    # Final summary
    print("\n" + "="*80)
    print("PROCESSING COMPLETE")
    print("="*80)
    print(f"Total time: {time.time() - overall_start:.1f} seconds")
    print("\nResults summary:")
    for sex_key, result in results.items():
        if result["success"]:
            print(f"  {sex_key}: ✓ {result['n_cells']} cells, {result['n_genes']} genes")
        else:
            print(f"  {sex_key}: ✗ Failed - {result['error']}")
    
    # Save overall summary
    with open(OUTPUT_DIR / "overall_summary.json", "w") as f:
        json.dump({
            "dataset": DATASET_NAME,
            "processing_time": time.time() - overall_start,
            "results": results,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }, f, indent=2)

if __name__ == "__main__":
    main()