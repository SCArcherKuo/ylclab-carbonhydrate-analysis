"""
Utility to retry failed identifiers and CIDs

This module provides functions to load and retry previously failed identifiers
that were saved during processing.
"""

import json
from pathlib import Path
from typing import Dict, List, Any, Optional
from loguru import logger
from . import config
from .pubchem_api import PubChemClient


def list_failed_files(failed_dir: Optional[str] = None) -> Dict[str, List[Path]]:
    """
    List all failed identifier files in the cache directory.
    
    Parameters:
    -----------
    failed_dir : str, optional
        Directory containing failed identifier files (defaults to config.FAILED_IDENTIFIERS_DIR)
    
    Returns:
    --------
    dict
        Dictionary with keys 'inchikey', 'smiles', 'cids' containing lists of file paths
    """
    cache_dir = Path(failed_dir or config.FAILED_IDENTIFIERS_DIR)
    
    if not cache_dir.exists():
        logger.warning(f"Cache directory does not exist: {cache_dir}")
        return {'inchikey': [], 'smiles': [], 'cids': []}
    
    failed_files = {
        'inchikey': list(cache_dir.glob('failed_inchikeys_*.json')),
        'smiles': list(cache_dir.glob('failed_smiless_*.json')),
        'cids': list(cache_dir.glob('failed_cids_*.json'))
    }
    
    return failed_files


def load_failed_identifiers(file_path: Path) -> Dict[str, Any]:
    """
    Load failed identifiers from a JSON file.
    
    Parameters:
    -----------
    file_path : Path
        Path to the failed identifiers JSON file
    
    Returns:
    --------
    dict
        Dictionary containing the failed data
    """
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    logger.info(f"Loaded {data.get('count', 0)} failed items from {file_path.name}")
    return data


def retry_failed_identifiers(
    file_path: Path,
    client: Optional[PubChemClient] = None
) -> Dict[str, Any]:
    """
    Retry failed identifiers from a saved file.
    
    Parameters:
    -----------
    file_path : Path
        Path to the failed identifiers JSON file
    client : PubChemClient, optional
        PubChem client to use (creates new one if not provided)
    
    Returns:
    --------
    dict
        Dictionary mapping identifiers to CIDs (or None if still failed)
    """
    data = load_failed_identifiers(file_path)
    
    if client is None:
        client = PubChemClient()
    
    identifier_type = data.get('identifier_type')
    identifiers = data.get('identifiers', [])
    
    if not identifiers:
        logger.warning("No identifiers to retry")
        return {}
    
    logger.info(f"Retrying {len(identifiers)} {identifier_type}s...")
    print(f"\nRetrying {len(identifiers)} failed {identifier_type}s...")
    
    # Clear the failed list for this retry attempt
    if identifier_type in client.failed_identifiers:
        client.failed_identifiers[identifier_type] = []
    
    results = client.resolve_identifiers_to_cids(identifiers, identifier_type)
    
    # Report results
    successful = sum(1 for cid in results.values() if cid is not None)
    still_failed = len(results) - successful
    
    logger.info(f"Retry complete: {successful} succeeded, {still_failed} still failed")
    print(f"\n✓ Retry complete: {successful} succeeded, {still_failed} still failed")
    
    return results


def retry_failed_cids(
    file_path: Path,
    client: Optional[PubChemClient] = None
) -> Dict[int, Dict[str, Any]]:
    """
    Retry fetching information for failed CIDs.
    
    Parameters:
    -----------
    file_path : Path
        Path to the failed CIDs JSON file
    client : PubChemClient, optional
        PubChem client to use (creates new one if not provided)
    
    Returns:
    --------
    dict
        Dictionary mapping CIDs to compound information
    """
    data = load_failed_identifiers(file_path)
    
    if client is None:
        client = PubChemClient()
    
    cids = data.get('cids', [])
    
    if not cids:
        logger.warning("No CIDs to retry")
        return {}
    
    logger.info(f"Retrying {len(cids)} CIDs...")
    print(f"\nRetrying {len(cids)} failed CIDs...")
    
    # Clear the failed list for this retry attempt
    client.failed_cids = []
    
    results = client.get_compound_info_batch(cids)
    
    # Report results
    successful = len(results)
    still_failed = len(cids) - successful
    
    logger.info(f"Retry complete: {successful} succeeded, {still_failed} still failed")
    print(f"\n✓ Retry complete: {successful} succeeded, {still_failed} still failed")
    
    return results


def main():
    """
    Interactive CLI for retrying failed identifiers and CIDs.
    """
    print("Failed Identifiers Retry Utility")
    print("=" * 50)
    
    # List all failed files
    failed_files = list_failed_files()
    
    all_files = []
    for file_type, files in failed_files.items():
        if files:
            print(f"\n{file_type.upper()} files:")
            for i, file_path in enumerate(files, 1):
                print(f"  [{len(all_files) + 1}] {file_path.name}")
                all_files.append((file_type, file_path))
    
    if not all_files:
        print("\nNo failed identifier files found.")
        return
    
    print(f"\nTotal: {len(all_files)} files")
    print("\nOptions:")
    print("  Enter file number to retry that file")
    print("  Enter 'all' to retry all files")
    print("  Enter 'q' to quit")
    
    choice = input("\nYour choice: ").strip().lower()
    
    if choice == 'q':
        return
    
    client = PubChemClient()
    
    if choice == 'all':
        for file_type, file_path in all_files:
            print(f"\n{'=' * 50}")
            print(f"Processing: {file_path.name}")
            print('=' * 50)
            
            if file_type == 'cids':
                retry_failed_cids(file_path, client)
            else:
                retry_failed_identifiers(file_path, client)
    else:
        try:
            idx = int(choice) - 1
            if 0 <= idx < len(all_files):
                file_type, file_path = all_files[idx]
                
                if file_type == 'cids':
                    retry_failed_cids(file_path, client)
                else:
                    retry_failed_identifiers(file_path, client)
            else:
                print("Invalid file number")
        except ValueError:
            print("Invalid input")


if __name__ == "__main__":
    main()
