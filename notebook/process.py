# %%
"""
Process metabolomics data and classify carbohydrates

This script provides functions to parse metabolomics Excel files and classify
compounds as carbohydrates using ChEBI ontology.
"""

# %%
import pandas as pd
import sys
from typing import List, Optional
from tqdm import tqdm
from loguru import logger
from datetime import datetime
from carbonhydrate_analysis import get_compound_info_pubchem
from carbonhydrate_analysis.utils import extract_ontology_terms_from_node

# Configure loguru logger
log_file = f"{datetime.now().strftime('%Y%m%d')}.log"

# Remove default handler
logger.remove()

# Add file handler
logger.add(
    log_file,
    rotation="100 MB",
    retention="10 days",
    level="DEBUG",
    format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {name}:{function}:{line} - {message}"
)

# Add console handler that uses tqdm.write() to not interfere with progress bars
logger.add(
    lambda msg: tqdm.write(msg, end=""),
    format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
    level="INFO",
    colorize=True
)

logger.info("="*80)
logger.info("Starting carbohydrate analysis processing")
logger.info("="*80)

# %%
def parse_metabolomics_xlsx(
    file_path: str,
    header_row: int = 4
) -> pd.DataFrame:
    """
    Parse metabolomics Excel file.
    
    Parameters:
    -----------
    file_path : str
        Path to the Excel file containing metabolomics data
    header_row : int, optional
        Row index where the header is located (default: 4)
    
    Returns:
    --------
    pd.DataFrame
        Complete dataframe with all columns
    
    Examples:
    ---------
    >>> df = parse_metabolomics_xlsx(
    ...     'data/raw/20251126_CH_Algae_T3_metabolomes_NEG_Area.xlsx'
    ... )
    >>> print(f"Total columns: {len(df.columns)}")
    """
    # Read the Excel file with proper header
    df = pd.read_excel(file_path, header=header_row)
    
    return df

# %%
def classify_carbohydrates(
    metadata_df: pd.DataFrame,
    inchikey_column: str = 'INCHIKEY'
) -> pd.DataFrame:
    """
    Classify compounds as carbohydrates using InChIKey and PubChem API.
    
    This function takes a metadata dataframe and adds three new columns:
    - 'Is Carbohydrate': Boolean indicating if the compound is a carbohydrate
    - 'Main Class': Main classification category
    - 'Subclass': Specific subclass name
    
    Parameters:
    -----------
    metadata_df : pd.DataFrame
        Dataframe containing metabolite metadata with an INCHIKEY column
    inchikey_column : str, optional
        Name of the column containing InChIKey information (default: 'INCHIKEY')
    
    Returns:
    --------
    pd.DataFrame
        Dataframe with three additional columns: 'Is Carbohydrate', 'Main Class', 'Subclass'
    
    Examples:
    ---------
    >>> metadata_df = parse_metabolomics_xlsx('data.xlsx')[1]
    >>> classified_df = classify_carbohydrates(metadata_df)
    >>> print(f"Carbohydrates found: {classified_df['Is Carbohydrate'].sum()}")
    """
    logger.info("Starting carbohydrate classification")
    logger.info(f"Input dataframe shape: {metadata_df.shape}")
    
    # Create a copy to avoid modifying the original
    result_df = metadata_df.copy()
    
    # Initialize new columns
    result_df['Is Carbohydrate'] = False
    result_df['Main Class'] = None
    result_df['Subclass'] = None
    
    # Collect all valid InChIKeys
    valid_inchikeys = []
    inchikey_to_idx = {}
    
    logger.info(f"Collecting valid InChIKeys from column '{inchikey_column}'")
    for idx, row in result_df.iterrows():
        inchikey_value = row.get(inchikey_column)
        
        # Skip if InChIKey is missing or NaN
        if pd.notna(inchikey_value) and isinstance(inchikey_value, str) and inchikey_value.strip():
            inchikey = inchikey_value.strip()
            valid_inchikeys.append(inchikey)
            if inchikey not in inchikey_to_idx:
                inchikey_to_idx[inchikey] = []
            inchikey_to_idx[inchikey].append(idx)
    
    if not valid_inchikeys:
        logger.warning("No valid InChIKeys found in the dataset")
        print("No valid InChIKeys found in the dataset")
        return result_df
    
    logger.info(f"Found {len(valid_inchikeys)} valid InChIKeys ({len(inchikey_to_idx)} unique)")
    print(f"Found {len(valid_inchikeys)} valid InChIKeys ({len(inchikey_to_idx)} unique)")
    
    # Get unique InChIKeys for batch query
    unique_inchikeys = list(inchikey_to_idx.keys())
    
    # Query PubChem in batch
    logger.info(f"Querying PubChem for {len(unique_inchikeys)} unique compounds")
    print(f"\nQuerying PubChem for {len(unique_inchikeys)} unique compounds...")
    
    try:
        pubchem_results = get_compound_info_pubchem(unique_inchikeys, identifier_type='inchikey')
        if pubchem_results is None:
            logger.error("PubChem query returned None")
            print("Error: PubChem query failed")
            return result_df
        if not isinstance(pubchem_results, list):
            logger.error(f"PubChem query returned unexpected type: {type(pubchem_results)}")
            print("Error: Unexpected response from PubChem")
            return result_df
        logger.info(f"PubChem query completed, received {len([r for r in pubchem_results if r is not None])} results")
    except Exception as e:
        logger.exception(f"Fatal error during PubChem query: {str(e)}")
        raise
    
    # Process results with progress bar
    logger.info("Processing classification results")
    print("\nProcessing results...")
    
    carb_count = 0
    processed_count = 0
    failed_count = 0
    
    for inchikey, compound_info in tqdm(
        zip(unique_inchikeys, pubchem_results),
        total=len(unique_inchikeys),
        desc="Classifying compounds"
    ):
        try:
            # Get all rows with this InChIKey
            row_indices = inchikey_to_idx[inchikey]
            
            if compound_info is not None:
                is_carb = compound_info.get('is_carbohydrate', False)
                main_class = compound_info.get('carbohydrate_main_class')
                subclass = compound_info.get('carbohydrate_subclass')
                
                # Update all rows with this InChIKey
                for idx in row_indices:
                    result_df.loc[idx, 'Is Carbohydrate'] = is_carb  # type: ignore
                    if is_carb:
                        result_df.loc[idx, 'Main Class'] = main_class  # type: ignore
                        result_df.loc[idx, 'Subclass'] = subclass  # type: ignore
                        carb_count += 1
                
                processed_count += 1
                logger.debug(f"Processed {inchikey[:20]}... -> is_carb={is_carb}, class={main_class}")
            else:
                failed_count += 1
                logger.warning(f"No compound info for {inchikey[:20]}...")
        
        except Exception as e:
            logger.exception(f"Error processing InChIKey {inchikey[:20]}...: {str(e)}")
            failed_count += 1
            continue
    
    logger.info(f"Classification complete: {processed_count} processed, {failed_count} failed, {carb_count} carbohydrates found")
    return result_df

# %%
def process_metabolomics_file(
    file_path: str,
    header_row: int = 4,
    inchikey_column: str = 'INCHIKEY'
) -> pd.DataFrame:
    """
    Complete processing pipeline: parse file and classify carbohydrates.
    
    This convenience function combines parsing and classification into one step.
    
    Parameters:
    -----------
    file_path : str
        Path to the Excel file containing metabolomics data
    header_row : int, optional
        Row index where the header is located (default: 4)
    inchikey_column : str, optional
        Name of the column containing InChIKey information (default: 'INCHIKEY')
    
    Returns:
    --------
    pd.DataFrame
        Complete dataframe with carbohydrate classification columns added
    
    Examples:
    ---------
    >>> df = process_metabolomics_file(
    ...     'data/raw/20251126_CH_Algae_T3_metabolomes_NEG_Area.xlsx'
    ... )
    >>> carb_count = df['Is Carbohydrate'].sum()
    >>> print(f"Found {carb_count} carbohydrates")
    """
    # Parse the file
    df = parse_metabolomics_xlsx(file_path, header_row)
    
    # Classify carbohydrates using InChIKey
    classified_df = classify_carbohydrates(df, inchikey_column)
    
    return classified_df

# %%
# Example usage
if __name__ == '__main__':
    import os
    import glob
    
    logger.info("Main script execution started")
    
    # Determine the correct path based on where the script is run from
    if os.path.exists('../data/raw'):
        raw_dir = '../data/raw'
        output_dir = '../data/processed'
    else:
        raw_dir = 'data/raw'
        output_dir = 'data/processed'
    
    logger.info(f"Raw data directory: {raw_dir}")
    logger.info(f"Output directory: {output_dir}")
    
    # Find all xlsx files in data/raw
    xlsx_files = glob.glob(os.path.join(raw_dir, '*.xlsx'))
    
    if not xlsx_files:
        logger.error(f"No Excel files found in {raw_dir}")
        print(f"No Excel files found in {raw_dir}")
        exit(1)
    
    logger.info(f"Found {len(xlsx_files)} Excel file(s) to process")
    print(f"Found {len(xlsx_files)} Excel file(s) to process:\n")
    for i, file in enumerate(xlsx_files, 1):
        print(f"  {i}. {os.path.basename(file)}")
        logger.info(f"  File {i}: {file}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each file
    total_carbohydrates = 0
    for file_idx, file_path in enumerate(xlsx_files, 1):
        print(f"\n{'='*80}")
        print(f"Processing file {file_idx}/{len(xlsx_files)}: {os.path.basename(file_path)}")
        print('='*80)
        
        logger.info(f"{'='*80}")
        logger.info(f"Processing file {file_idx}/{len(xlsx_files)}: {file_path}")
        logger.info(f"{'='*80}")
        
        try:
            logger.info(f"Calling process_metabolomics_file for {file_path}")
            df = process_metabolomics_file(file_path)
            
            logger.info(f"File structure - Total columns: {len(df.columns)}, Total rows: {len(df)}")
            print(f"\nFile structure:")
            print(f"  Total columns: {len(df.columns)}")
            print(f"  Total rows: {len(df)}")
            
            carb_count = df['Is Carbohydrate'].sum()
            logger.info(f"Metabolite classification - Total: {len(df)}, Carbohydrates: {carb_count}")
            print(f"\nMetabolite classification:")
            print(f"  Total metabolites: {len(df)}")
            print(f"  Carbohydrates found: {carb_count}")
            total_carbohydrates += carb_count
            
            if carb_count > 0:
                print(f"\nCarbohydrate breakdown:")
                class_counts = df[df['Is Carbohydrate']]['Main Class'].value_counts()
                for class_name, count in class_counts.items():
                    print(f"  {class_name}: {count}")
                    logger.info(f"  {class_name}: {count}")
            
            # Extract base filename without extension
            base_filename = os.path.splitext(os.path.basename(file_path))[0]
            
            # Save classified data
            classified_output = os.path.join(output_dir, f"{base_filename}_classified.csv")
            df.to_csv(classified_output, index=False)
            logger.info(f"Saved classified data to: {classified_output}")
            print(f"\n✓ Saved classified data to: {classified_output}")
            
            # Save carbohydrates only (if any found)
            if carb_count > 0:
                carb_only = df[df['Is Carbohydrate']].copy()
                carb_output = os.path.join(output_dir, f"{base_filename}_carbohydrates_only.csv")
                carb_only.to_csv(carb_output, index=False)
                logger.info(f"Saved carbohydrates only to: {carb_output}")
                print(f"✓ Saved carbohydrates only to: {carb_output}")
            
            logger.info(f"Successfully completed processing of {file_path}")
            
        except KeyboardInterrupt:
            logger.warning(f"Keyboard interrupt received while processing {file_path}")
            print(f"\n✗ Processing interrupted by user")
            raise
        except Exception as e:
            logger.exception(f"Error processing {file_path}: {str(e)}")
            print(f"✗ Error processing {os.path.basename(file_path)}: {str(e)}")
            continue
    
    print(f"\n{'='*80}")
    print(f"SUMMARY")
    print('='*80)
    print(f"Files processed: {len(xlsx_files)}")
    print(f"Total carbohydrates found: {total_carbohydrates}")
    print(f"All results saved to {output_dir}/")
    
    logger.info(f"{'='*80}")
    logger.info(f"FINAL SUMMARY")
    logger.info(f"Files processed: {len(xlsx_files)}")
    logger.info(f"Total carbohydrates found: {total_carbohydrates}")
    logger.info(f"All results saved to {output_dir}/")
    logger.info(f"{'='*80}")
    logger.info("Script execution completed successfully")

# %%
