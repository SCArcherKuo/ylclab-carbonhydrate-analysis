"""
Carbohydrate Analysis Package

A package for identifying and classifying carbohydrates in metabolomics data
using PubChem and ChEBI databases.
"""

from .pubchem_api import get_compound_info_pubchem, PubChemClient
from .chebi_api import ChEBIClient, get_chebi_children, get_main_groups
from .classification import CarbohydrateClassifier, classify_carbohydrate
from .utils import extract_term_string, extract_ontology_terms_from_node
from .retry_failed import (
    list_failed_files,
    load_failed_identifiers,
    retry_failed_identifiers,
    retry_failed_cids
)

__version__ = "0.1.0"

__all__ = [
    # Main API function
    'get_compound_info_pubchem',
    
    # Clients
    'PubChemClient',
    'ChEBIClient',
    'CarbohydrateClassifier',
    
    # ChEBI functions
    'get_chebi_children',
    'get_main_groups',
    
    # Classification
    'classify_carbohydrate',
    
    # Utilities
    'extract_term_string',
    'extract_ontology_terms_from_node',
    
    # Retry utilities
    'list_failed_files',
    'load_failed_identifiers',
    'retry_failed_identifiers',
    'retry_failed_cids',
]
