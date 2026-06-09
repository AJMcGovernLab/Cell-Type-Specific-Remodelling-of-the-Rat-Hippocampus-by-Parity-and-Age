#!/usr/bin/env python3
"""
CORRECTED filtering script for Yao Hippocampus 10x dataset
Properly identifies hippocampal cells and handles large expression matrix
"""

import pandas as pd
import numpy as np
import h5py
from pathlib import Path
import time
import json

# === CONFIGURATION ===
DATASET_NAME = "yao_hippo_10x"
METADATA_PATH = r"D:/1Reference Datasets/A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain 2023/CTX_Hip_anno_10x.csv"
EXPRESSION_PATH = r"D:/1Reference Datasets/A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain 2023/CTX_Hip_counts_10x.h5"
OUTPUT_DIR = Path(f"../../02_filtered_references/{DATASET_NAME}")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# === HIPPOCAMPUS IDENTIFICATION ===
# Use cluster patterns since this dataset has different region labels
HIPPO_CLUSTER_PATTERNS = [
    # Numbered patterns from Allen Brain Atlas
    r"333_CA1-ProS", r"342_CA1", r"343_CA1", r"344_CA1",
    r"359_CA2-IG-FC", r"360_CA2-IG-FC",
    r"354_CA3",
    r"362_DG", r"363_DG", r"364_DG",
    r"294_CT SUB", r"273_NP SUB", r"320_SUB",
    # General patterns
    r"\d+_CA[123]", r"\d+_DG", r"\d+_SUB",
    r"CA[123]-", r"DG-", r"SUB-",
    # Interneuron patterns specific to hippocampus
    r"\d+_Sst Chodl"  # Hippocampus-enriched
]

# Subclass labels for hippocampus
HIPPO_SUBCLASS_LABELS = [
    "CA1-ProS", "CA2-IG-FC", "CA3", "DG", 
    "SUB-ProS", "CT SUB", "NP SUB"
]

SEX_LABELS = {"female": "F", "male": "M", "mixed": None}

def identify_hippocampal_cells(metadata):
    """Identify hippocampal cells using cluster and subclass patterns"""
    
    print("\nIdentifying hippocampal cells...")
    
    # Initialize mask
    hippo_mask = pd.Series(False, index=metadata.index)
    
    # Method 1: Check cluster_label patterns
    if "cluster_label" in metadata.columns:
        print("  Checking cluster_label patterns...")
        cluster_labels = metadata["cluster_label"].astype(str)
        
        # First try exact matches for known hippocampal clusters
        exact_matches = [
            "333_CA1-ProS", "342_CA1", "343_CA1", "344_CA1",
            "359_CA2-IG-FC", "360_CA2-IG-FC",
            "354_CA3", "362_DG", "363_DG", "364_DG",
            "294_CT SUB", "273_NP SUB", "320_SUB"
        ]
        
        for cluster in exact_matches:
            exact_mask = cluster_labels == cluster
            n_found = exact_mask.sum()
            if n_found > 0:
                print(f"    Cluster {cluster}: {n_found} cells")
                hippo_mask |= exact_mask
        
        # Then check patterns
        for pattern in HIPPO_CLUSTER_PATTERNS:
            pattern_mask = cluster_labels.str.contains(
                pattern, case=False, na=False, regex=True
            )
            n_found = pattern_mask.sum() - hippo_mask.sum()  # New cells found
            if n_found > 0:
                print(f"    Pattern '{pattern}': +{n_found} new cells")
                hippo_mask |= pattern_mask
    
    # Method 2: Check subclass_label
    if "subclass_label" in metadata.columns:
        print("  Checking subclass_label...")
        subclass_mask = metadata["subclass_label"].isin(HIPPO_SUBCLASS_LABELS)
        n_new = (subclass_mask & ~hippo_mask).sum()
        if n_new > 0:
            print(f"    Subclass labels: +{n_new} new cells")
            hippo_mask |= subclass_mask
    
    # Method 3: Check if any cells from hippocampal regions got missed
    if "region_label" in metadata.columns:
        # Even though regions are labeled differently, check for subiculum
        if "SUB" in metadata["region_label"].unique() or "ProS" in metadata["region_label"].unique():
            sub_pattern = metadata["region_label"].astype(str).str.contains("SUB|ProS", na=False)
            # Only add if cluster also suggests hippocampus
            if "cluster_label" in metadata.columns:
                sub_pattern = sub_pattern & metadata["cluster_label"].astype(str).str.contains(
                    "SUB|ProS|CA|DG", na=False
                )
            n_new = (sub_pattern & ~hippo_mask).sum()
            if n_new > 0:
                print(f"    Subiculum regions: +{n_new} new cells")
                hippo_mask |= sub_pattern
    
    total_hippo = hippo_mask.sum()
    print(f"\nTotal hippocampal cells identified: {total_hippo}")
    
    # Validate
    if total_hippo > 0 and "cluster_label" in metadata.columns:
        hippo_clusters = metadata.loc[hippo_mask, "cluster_label"].value_counts().head(20)
        print("\nTop 20 hippocampal clusters:")
        for cluster, count in hippo_clusters.items():
            print(f"  {cluster}: {count}")
    
    return hippo_mask

def load_yao_expression_h5(h5_path, cell_names):
    """Load expression data from Yao HDF5 file"""
    
    print("\nLoading expression matrix from HDF5...")
    
    with h5py.File(h5_path, 'r') as f:
        # Explore structure
        print("  HDF5 structure:")
        def print_structure(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"    {name}: {obj.shape} (Dataset)")
            elif isinstance(obj, h5py.Group):
                print(f"    {name}: (Group)")
        f.visititems(print_structure)
        def print_structure(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"    {name}: {obj.shape} (Dataset)")
            elif isinstance(obj, h5py.Group):
                print(f"    {name}: (Group)")
        f.visititems(print_structure)
        
        # Common HDF5 structures for 10x data
        # Try different possible paths
        matrix_paths = [
            'matrix',
            'data/matrix',
            'data/counts',  # This is the actual path in this file
            'X',
            'counts'
        ]
        
        matrix_data = None
        matrix_path = None
        for path in matrix_paths:
            if path in f:
                # Check if it's actually a dataset (not a group)
                if isinstance(f[path], h5py.Dataset):
                    matrix_data = f[path]
                    matrix_path = path
                    print(f"  Found expression matrix at: {path}")
                    break
        
        if matrix_data is None:
            raise ValueError("Could not find expression matrix in HDF5 file")
        
        # Get cell barcodes
        barcode_paths = [
            'barcodes',
            'data/barcodes',
            'data/samples',  # Added this path
            'obs_names',
            'cell_names',
            'data/cell_names'
        ]
        
        all_barcodes = None
        for path in barcode_paths:
            try:
                if path in f:
                    barcode_data = f[path][:]
                    if len(barcode_data) > 0:
                        # Handle fixed-length byte strings (|S46 dtype)
                        if barcode_data.dtype.kind == 'S':  # Fixed-length bytes
                            all_barcodes = [b.decode('utf-8').strip() for b in barcode_data]
                        elif isinstance(barcode_data[0], bytes):
                            all_barcodes = [b.decode('utf-8') for b in barcode_data]
                        else:
                            all_barcodes = list(barcode_data)
                        print(f"  Found cell barcodes at: {path}")
                        print(f"  Total cells in matrix: {len(all_barcodes)}")
                        break
            except Exception as e:
                print(f"  Error reading {path}: {e}")
                continue
        
        if all_barcodes is None:
            raise ValueError("Could not find cell barcodes in HDF5 file")
        
        # Get gene names
        gene_paths = [
            'genes',
            'data/genes',
            'data/gene',  # This is the actual path in this file
            'var_names',
            'features',
            'data/features',
            'gene_names'
        ]
        
        gene_names = None
        for path in gene_paths:
            try:
                if path in f:
                    gene_data = f[path][:]
                    if len(gene_data) > 0:
                        # Handle fixed-length byte strings (|S30 dtype)
                        if gene_data.dtype.kind == 'S':  # Fixed-length bytes
                            gene_names = [g.decode('utf-8').strip() for g in gene_data]
                        elif isinstance(gene_data[0], bytes):
                            gene_names = [g.decode('utf-8') for g in gene_data]
                        else:
                            gene_names = list(gene_data)
                        print(f"  Found gene names at: {path}")
                        print(f"  Total genes: {len(gene_names)}")
                        break
            except:
                continue
        
        if gene_names is None:
            # If no gene names, create generic ones
            n_genes = matrix_data.shape[0] if matrix_data.shape[0] < matrix_data.shape[1] else matrix_data.shape[1]
            gene_names = [f"Gene_{i}" for i in range(n_genes)]
            print(f"  No gene names found, using generic names")
        
        # Create barcode lookup
        barcode_to_idx = {bc: idx for idx, bc in enumerate(all_barcodes)}
        
        # Find matching cells
        matching_indices = []
        matching_cells = []
        unmatched = []
        
        # Try exact matches first
        for cell in cell_names:
            if cell in barcode_to_idx:
                matching_indices.append(barcode_to_idx[cell])
                matching_cells.append(cell)
            else:
                unmatched.append(cell)
        
        # If no exact matches, try to match by suffix (common with 10x barcodes)
        if len(matching_indices) == 0 and len(unmatched) > 0:
            print("  No exact matches found, trying barcode suffix matching...")
            
            # Show examples to help debug
            print("\n  Sample name examples from metadata:")
            for i, cell in enumerate(cell_names[:3]):
                print(f"    {cell}")
            print("\n  Barcode examples from expression matrix:")
            for i, bc in enumerate(all_barcodes[:3]):
                print(f"    {bc}")
            
            # Extract barcode part from sample names
            # Format: 10X_cells.ACGCAGCAGACCGGAT-L8TX_180221_01_C11
            for cell in unmatched:
                if '.' in cell and '-' in cell:
                    # Extract barcode between . and -
                    try:
                        barcode = cell.split('.')[1].split('-')[0]
                        # Look for this barcode in the matrix
                        for idx, bc in enumerate(all_barcodes):
                            if barcode in bc or bc in barcode:
                                matching_indices.append(idx)
                                matching_cells.append(cell)
                                break
                    except:
                        continue
                        
            # If still no matches, try matching the full sample name
            if len(matching_indices) == 0:
                print("  Still no matches, checking if sample names are substrings...")
                for i, cell in enumerate(cell_names[:100]):  # Check first 100
                    for idx, bc in enumerate(all_barcodes):
                        if cell in bc or bc in cell:
                            matching_indices.append(idx)
                            matching_cells.append(cell)
                            if len(matching_indices) >= 10:  # Stop after finding some matches
                                print(f"    Found match: '{cell}' <-> '{bc}'")
                            break
        
        print(f"  Matched {len(matching_indices)} cells out of {len(cell_names)}")
        
        if len(matching_indices) == 0:
            # Show examples to help debug
            print("\n  Sample name examples from metadata:")
            for i, cell in enumerate(cell_names[:3]):
                print(f"    {cell}")
            print("\n  Barcode examples from expression matrix:")
            for i, bc in enumerate(all_barcodes[:3]):
                print(f"    {bc}")
            raise ValueError("No matching cells found between metadata and expression matrix")
        
        # Load expression data for matched cells
        print(f"  Loading expression data for {len(matching_indices)} cells...")
        
        # The matrix is genes x cells (31053 x 1169213)
        # We need to extract specific columns (cells)
        
        # For HDF5, we can use fancy indexing to load only the cells we need
        # This is much more memory efficient than loading the entire matrix
        print(f"  Matrix shape: {matrix_data.shape}")
        
        # Sort indices for more efficient HDF5 access
        sorted_idx = np.argsort(matching_indices)
        sorted_indices = np.array(matching_indices)[sorted_idx]
        
        # Load data in chunks if there are many cells
        if len(sorted_indices) > 10000:
            print(f"  Loading in chunks due to large number of cells...")
            chunk_size = 5000
            expression_chunks = []
            
            for i in range(0, len(sorted_indices), chunk_size):
                end_idx = min(i + chunk_size, len(sorted_indices))
                chunk_indices = sorted_indices[i:end_idx]
                chunk_data = matrix_data[:, chunk_indices]
                expression_chunks.append(chunk_data)
                if (i + chunk_size) % 10000 == 0:
                    print(f"    Loaded {i + chunk_size} cells...")
            
            # Combine chunks
            expression_matrix = np.hstack(expression_chunks)
            
            # Restore original order
            unsort_idx = np.empty_like(sorted_idx)
            unsort_idx[sorted_idx] = np.arange(len(sorted_idx))
            expression_matrix = expression_matrix[:, unsort_idx]
        else:
            # Load all at once for smaller datasets
            expression_matrix = matrix_data[:, matching_indices]
        
        # Convert to numpy array if needed
        expression_matrix = np.array(expression_matrix, dtype=np.float32)
        
        print(f"  Final expression matrix: {expression_matrix.shape[0]} genes x {expression_matrix.shape[1]} cells")
    
    return expression_matrix, gene_names, matching_cells

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
        region_counts = metadata["region_label"].value_counts()
        for region, count in region_counts.head(10).items():
            print(f"  {region}: {count}")
    
    # Identify hippocampal cells
    hippo_mask = identify_hippocampal_cells(metadata)
    
    if hippo_mask.sum() == 0:
        print("ERROR: No hippocampal cells found!")
        return None, None
    sex_col = None
    # Apply sex filter
    if sex_label:
        if "sex_label" in metadata.columns:
            sex_col = "sex_label"
        elif "donor_sex_label" in metadata.columns:
            sex_col = "donor_sex_label"
        else:
            print("WARNING: No sex column found, creating mixed dataset")
            sex_mask = pd.Series(True, index=metadata.index)
            sex_col = None
        
        if sex_col:
            sex_mask = metadata[sex_col] == sex_label
            print(f"\nFiltering for sex='{sex_label}'")
            print(f"  Sex-specific cells: {sex_mask.sum()}")
        
        final_mask = hippo_mask & sex_mask
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
    
    # Get sample names
    sample_names = filtered_metadata["sample_name"].tolist()
    
    # Process expression data
    try:
        expression_matrix, gene_names, valid_samples = load_yao_expression_h5(
            EXPRESSION_PATH, sample_names
        )
        
        # Update metadata to match valid samples
        filtered_metadata = filtered_metadata[filtered_metadata["sample_name"].isin(valid_samples)]
        
    except Exception as e:
        print(f"ERROR processing expression matrix: {e}")
        # Save metadata anyway
        metadata_output = OUTPUT_DIR / f"metadata_{sex_key}.csv"
        filtered_metadata.to_csv(metadata_output, index=False)
        return filtered_metadata, None
    
    # Save outputs
    print("\nSaving filtered data...")
    
    # Save metadata
    metadata_output = OUTPUT_DIR / f"metadata_{sex_key}.csv"
    filtered_metadata.to_csv(metadata_output, index=False)
    print(f"  Metadata saved: {metadata_output}")
    
    # Save expression as HDF5 (much more efficient than CSV)
    h5_output = OUTPUT_DIR / f"expression_matrix_{sex_key}.h5"
    print(f"  Saving expression matrix to HDF5...")
    
    with h5py.File(h5_output, 'w') as f:
        # Save expression data with compression
        f.create_dataset('data/expression', 
                        data=expression_matrix.astype(np.float32), 
                        compression='gzip',
                        compression_opts=4)  # Medium compression
        
        # Save gene names
        gene_names_encoded = [str(g).encode('utf-8') for g in gene_names]
        f.create_dataset('data/gene_names', data=gene_names_encoded)
        
        # Save cell names
        cell_names_encoded = [c.encode('utf-8') for c in valid_samples]
        f.create_dataset('data/cell_names', data=cell_names_encoded)
        
        # Add metadata attributes
        f.attrs['n_genes'] = len(gene_names)
        f.attrs['n_cells'] = len(valid_samples)
        f.attrs['dataset'] = DATASET_NAME
        f.attrs['sex_filter'] = sex_key
        f.attrs['created'] = time.strftime("%Y-%m-%d %H:%M:%S")
    
    print(f"  Expression matrix saved: {h5_output}")
    
    # Save summary
    summary = {
        "dataset": DATASET_NAME,
        "sex_filter": sex_key,
        "total_cells_original": len(metadata),
        "hippocampal_cells_identified": int(hippo_mask.sum()),
        "cells_after_sex_filter": len(filtered_metadata),
        "cells_with_expression": len(valid_samples),
        "genes": len(gene_names),
        "file_size_mb": h5_output.stat().st_size / (1024 * 1024),
        "hippocampal_identification_method": "cluster_label exact matches + patterns",
        "validation": {
            "region_distribution": filtered_metadata["region_label"].value_counts().head(10).to_dict() 
                if "region_label" in filtered_metadata else {},
            "top_clusters": filtered_metadata["cluster_label"].value_counts().head(20).to_dict() 
                if "cluster_label" in filtered_metadata else {},
            "sex_distribution": filtered_metadata[sex_col].value_counts().to_dict() 
                if sex_col and sex_col in filtered_metadata else {}
        }
    }
    
    summary_output = OUTPUT_DIR / f"filtering_summary_{sex_key}.json"
    with open(summary_output, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"  Summary saved: {summary_output}")
    
    return filtered_metadata, True

def main():
    """Main execution"""
    print("="*80)
    print("CORRECTED HIPPOCAMPAL FILTERING FOR YAO HIPPOCAMPUS 10X")
    print("="*80)
    print("\nNOTE: This will create HDF5 files instead of CSV to handle the large matrices")
    
    overall_start = time.time()
    
    # Process each sex group
    results = {}
    for sex_key, sex_label in SEX_LABELS.items():
        try:
            metadata, success = save_filtered_data(sex_key, sex_label)
            if metadata is not None:
                results[sex_key] = {
                    "success": success is not None,
                    "n_cells": len(metadata),
                    "expression_saved": success is not None
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
        if result.get("success", False):
            status = "✓ (with expression)" if result.get("expression_saved") else "⚠ (metadata only)"
            print(f"  {sex_key}: {status} {result['n_cells']} cells")
        else:
            error_msg = result.get('error', 'Unknown error')
            print(f"  {sex_key}: ✗ Failed - {error_msg}")
    
    print("\nNOTE: The old CSV files can be deleted after verifying the new HDF5 files work correctly")

if __name__ == "__main__":
    main()