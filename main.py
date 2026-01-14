"""
Carbohydrate Filtering for Metabolomics Data

This module provides functions to query compound information from PubChem API
and identify carbohydrates based on chemical classifications and synonyms.
"""

import requests
import json
import time
from typing import Optional, Dict, Any, List, Tuple
from urllib.parse import quote
import config

# Cache for ChEBI API responses to avoid repeated calls
_chebi_children_cache = {}


def get_chebi_children(chebi_id: int) -> List[Dict[str, Any]]:
    """
    Get direct children of a ChEBI ID from ChEBI backend API.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID (without 'CHEBI:' prefix)
    
    Returns:
    --------
    list of dict
        List of direct children with 'is a' relationship.
        Each dict contains: 'init_id', 'init_name', 'relation_type'
    """
    # Check cache first
    if chebi_id in _chebi_children_cache:
        return _chebi_children_cache[chebi_id]
    
    try:
        url = f'{config.CHEBI_API_BASE_URL}/ontology/children/{chebi_id}/'
        response = requests.get(url, headers={'accept': '*/*'}, timeout=config.API_TIMEOUT)
        
        if response.status_code != 200:
            print(f"Warning: ChEBI API failed for {chebi_id}: HTTP {response.status_code}")
            _chebi_children_cache[chebi_id] = []
            return []
        
        data = response.json()
        incoming = data.get('ontology_relations', {}).get('incoming_relations', [])
        
        # Filter for "is a" relations (direct children)
        children = [rel for rel in incoming if rel.get('relation_type') == 'is a']
        
        _chebi_children_cache[chebi_id] = children
        return children
        
    except Exception as e:
        print(f"Warning: Error fetching ChEBI children for {chebi_id}: {str(e)}")
        _chebi_children_cache[chebi_id] = []
        return []


def get_main_groups(chebi_id: int) -> List[str]:
    """
    Get main groups (children with >1 children) of a ChEBI ID.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID to get main groups for
    
    Returns:
    --------
    list of str
        List of main group names (children that have >1 children themselves)
    """
    children = get_chebi_children(chebi_id)
    main_groups = []
    
    for child in children:
        child_id = child.get('init_id')
        child_name = child.get('init_name', '')
        
        if child_id:
            # Check if this child has more than 1 children
            grandchildren = get_chebi_children(child_id)
            if len(grandchildren) > 1:
                main_groups.append(child_name)
    
    return main_groups


def extract_term_string(term: Any) -> str:
    """
    Extract string from StringWithMarkup format or return as-is.
    
    Parameters:
    -----------
    term : any
        Term that might be in StringWithMarkup format or plain string
    
    Returns:
    --------
    str
        Extracted string value
    """
    if isinstance(term, dict):
        if 'StringWithMarkup' in term:
            string_data = term['StringWithMarkup']
            if isinstance(string_data, dict) and 'String' in string_data:
                return string_data['String']
        # Try direct String key
        if 'String' in term:
            return term['String']
    return str(term)


def classify_carbohydrate(chebi_ontology: List[Any]) -> Tuple[Optional[str], Optional[str]]:
    """
    Classify a compound based on its ChEBI ontology terms.
    
    Parameters:
    -----------
    chebi_ontology : list
        List of ChEBI ontology terms (may be strings or StringWithMarkup dicts)
    
    Returns:
    --------
    tuple of (main_class, subclass)
        main_class: One of ['main carbohydrate group', 'other carbohydrate', 
                           'main carbohydrate derivative group', 'other carbohydrate derivative', 
                           'other', None]
        subclass: Specific subclass name or None
    """
    if not chebi_ontology:
        return None, None
    
    # Extract clean strings from all terms
    clean_terms = [extract_term_string(term) for term in chebi_ontology]
    clean_terms_lower = [t.lower() for t in clean_terms]
    
    # Check if compound belongs to carbohydrates and carbohydrate derivatives (CHEBI:78616)
    if 'carbohydrates and carbohydrate derivatives' not in clean_terms_lower:
        return None, None
    
    # Get main groups for carbohydrate and carbohydrate derivative
    carb_main_groups = get_main_groups(config.CHEBI_CARBOHYDRATE_ID)
    carb_deriv_main_groups = get_main_groups(config.CHEBI_CARBOHYDRATE_DERIVATIVE_ID)
    
    carb_main_groups_lower = [g.lower() for g in carb_main_groups]
    carb_deriv_main_groups_lower = [g.lower() for g in carb_deriv_main_groups]
    
    # Check for carbohydrate (CHEBI:16646)
    if 'carbohydrate' in clean_terms_lower:
        # Check if matches main carbohydrate groups
        for i, term_lower in enumerate(clean_terms_lower):
            if term_lower in carb_main_groups_lower:
                idx = carb_main_groups_lower.index(term_lower)
                return 'main carbohydrate group', carb_main_groups[idx]
        
        # If no main group match, it's "other carbohydrate"
        # Get all direct children and find which one matches
        carb_children = get_chebi_children(config.CHEBI_CARBOHYDRATE_ID)
        carb_children_names = [child['init_name'].lower() for child in carb_children]
        
        # Find the direct child that appears in the ontology terms
        for term in clean_terms:
            term_lower = term.lower()
            if term_lower in carb_children_names:
                # Return the original capitalization from the ontology
                return 'other carbohydrate', term
        
        return 'other carbohydrate', 'carbohydrate'
    
    # Check for carbohydrate derivative (CHEBI:63299)
    if 'carbohydrate derivative' in clean_terms_lower:
        # Check if matches main carbohydrate derivative groups
        for i, term_lower in enumerate(clean_terms_lower):
            if term_lower in carb_deriv_main_groups_lower:
                idx = carb_deriv_main_groups_lower.index(term_lower)
                return 'main carbohydrate derivative group', carb_deriv_main_groups[idx]
        
        # If no main group match, it's "other carbohydrate derivative"
        # Get all direct children and find which one matches
        carb_deriv_children = get_chebi_children(config.CHEBI_CARBOHYDRATE_DERIVATIVE_ID)
        carb_deriv_children_names = [child['init_name'].lower() for child in carb_deriv_children]
        
        # Find the direct child that appears in the ontology terms
        for term in clean_terms:
            term_lower = term.lower()
            if term_lower in carb_deriv_children_names:
                # Return the original capitalization from the ontology
                return 'other carbohydrate derivative', term
        
        return 'other carbohydrate derivative', 'carbohydrate derivative'
    
    # If under root but not in carbohydrate or carbohydrate derivative
    # This is the "other" category
    # Get all direct children and find which one matches
    root_children = get_chebi_children(config.CHEBI_ROOT_ID)
    root_children_names = [child['init_name'].lower() for child in root_children]
    
    # Find the direct child that appears in the ontology terms
    for term in clean_terms:
        term_lower = term.lower()
        if term_lower in root_children_names and term_lower != 'carbohydrates and carbohydrate derivatives':
            # Return the original capitalization from the ontology
            return 'other', term
    
    return 'other', 'carbohydrates and carbohydrate derivatives'


def get_compound_info_pubchem(identifier: str, identifier_type: str = 'auto') -> Optional[Dict[str, Any]]:
    """
    Fetch compound information from PubChem using InChIKey or SMILES.
    
    Parameters:
    -----------
    identifier : str
        The InChIKey or SMILES string of the compound
    identifier_type : str, optional
        Type of identifier: 'inchikey', 'smiles', or 'auto' (default: 'auto')
        'auto' will detect based on format
    
    Returns:
    --------
    dict or None
        Dictionary containing compound information:
        - pubchem_cid: PubChem Compound ID
        - name: IUPAC name
        - formula: Molecular formula
        - molecular_weight: Molecular weight
        - inchi: InChI string
        - inchikey: InChIKey
        - smiles: Canonical SMILES
        - classifications: List of chemical classification hierarchies
        - chebi_ontology: ChEBI ontology terms extracted from classifications
        - synonyms: List of compound synonyms (max 20)
        - is_carbohydrate: Boolean flag (True if compound belongs to CHEBI:78616)
        - carbohydrate_main_class: Main classification category (5 categories or None)
        - carbohydrate_subclass: Subclass name or None
    
    Example:
    --------
    >>> result = get_compound_info_pubchem("RBNPOMFGQQGHHO-UHFFFAOYSA-N")
    >>> print(f"Compound: {result['name']}, Is carbohydrate: {result['is_carbohydrate']}")
    >>> if result['chebi_ontology']:
    >>>     print(f"ChEBI Terms: {result['chebi_ontology']}")
    """
    # Auto-detect identifier type
    if identifier_type == 'auto':
        if '-' in identifier and len(identifier.split('-')) == 3:
            identifier_type = 'inchikey'
        else:
            identifier_type = 'smiles'
    
    try:
        base_url = config.PUBCHEM_BASE_URL
        
        # Search for compound
        if identifier_type == 'inchikey':
            search_url = f"{base_url}/compound/inchikey/{identifier}/JSON"
        else:  # SMILES
            # URL encode SMILES to handle special characters
            encoded_smiles = quote(identifier, safe='')
            search_url = f"{base_url}/compound/smiles/{encoded_smiles}/JSON"
        
        response = requests.get(search_url, timeout=config.API_TIMEOUT)
        
        if response.status_code != 200:
            print(f"PubChem search failed: HTTP {response.status_code}")
            return None
        
        data = response.json()
        
        if 'PC_Compounds' not in data or len(data['PC_Compounds']) == 0:
            print(f"No compound found in PubChem for {identifier}")
            return None
        
        compound = data['PC_Compounds'][0]
        cid = compound['id']['id']['cid']
        
        print(f"Found PubChem CID: {cid}")
        
        # Get properties
        props_url = f"{base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,InChI,InChIKey,CanonicalSMILES,IUPACName/JSON"
        props_response = requests.get(props_url, timeout=config.API_TIMEOUT)
        
        properties = {}
        if props_response.status_code == 200:
            props_data = props_response.json()
            if 'PropertyTable' in props_data and 'Properties' in props_data['PropertyTable']:
                properties = props_data['PropertyTable']['Properties'][0]
        
        # Get classification (with error handling for malformed JSON)
        time.sleep(config.RATE_LIMIT_DELAY)  # Rate limiting
        class_url = f"{base_url}/compound/cid/{cid}/classification/JSON"
        class_response = requests.get(class_url, timeout=config.API_TIMEOUT)
        
        classifications = []
        if class_response.status_code == 200:
            try:
                class_data = class_response.json()
                if 'Hierarchies' in class_data and 'Hierarchy' in class_data['Hierarchies']:
                    for hierarchy in class_data['Hierarchies']['Hierarchy']:
                        if 'Node' in hierarchy:
                            classifications.append(hierarchy)
            except json.JSONDecodeError:
                print(f"Warning: Classification data has JSON errors (PubChem issue), skipping...")
        
        # Get synonyms
        time.sleep(config.RATE_LIMIT_DELAY)  # Rate limiting
        syn_url = f"{base_url}/compound/cid/{cid}/synonyms/JSON"
        syn_response = requests.get(syn_url, timeout=config.API_TIMEOUT)
        
        synonyms = []
        if syn_response.status_code == 200:
            syn_data = syn_response.json()
            if 'InformationList' in syn_data and 'Information' in syn_data['InformationList']:
                if len(syn_data['InformationList']['Information']) > 0:
                    synonyms = syn_data['InformationList']['Information'][0].get('Synonym', [])
        
        # Extract ChEBI ontology terms from classifications
        chebi_ontology = []
        
        def extract_ontology_terms(node, terms_list):
            """Recursively extract ontology terms from classification nodes."""
            if isinstance(node, dict):
                if 'Information' in node:
                    info = node['Information']
                    if 'Name' in info:
                        terms_list.append(info['Name'])
                
                # Recursively process children
                if 'Children' in node and 'Node' in node['Children']:
                    children = node['Children']['Node']
                    if isinstance(children, list):
                        for child in children:
                            extract_ontology_terms(child, terms_list)
                    else:
                        extract_ontology_terms(children, terms_list)
        
        for classification in classifications:
            # Check if this is ChEBI classification
            source = classification.get('SourceName', '')
            if 'chebi' in source.lower() or 'ChEBI' in source:
                # Extract all terms from this hierarchy
                if 'Node' in classification:
                    nodes = classification['Node'] if isinstance(classification['Node'], list) else [classification['Node']]
                    for node in nodes:
                        extract_ontology_terms(node, chebi_ontology)
        
        # Classify carbohydrate based on ChEBI ontology
        main_class, subclass = classify_carbohydrate(chebi_ontology)
        
        # Determine is_carbohydrate flag based on CHEBI:78616
        is_carbohydrate = main_class is not None
        
        # Compile results
        result = {
            'pubchem_cid': cid,
            'name': properties.get('IUPACName', 'Unknown'),
            'formula': properties.get('MolecularFormula'),
            'molecular_weight': properties.get('MolecularWeight'),
            'inchi': properties.get('InChI'),
            'inchikey': properties.get('InChIKey'),
            'smiles': properties.get('CanonicalSMILES'),
            'classifications': classifications,
            'chebi_ontology': chebi_ontology,
            'synonyms': synonyms[:config.MAX_SYNONYMS],
            'is_carbohydrate': is_carbohydrate,
            'carbohydrate_main_class': main_class,
            'carbohydrate_subclass': subclass
        }
        
        return result
        
    except requests.exceptions.Timeout:
        print("PubChem request timeout")
        return None
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """
    Main function for testing the compound lookup functionality.
    """
    print("Carbohydrate Filtering Tool")
    print("=" * 50)
    
    # Example usage
    test_inchikey = "RBNPOMFGQQGHHO-UHFFFAOYSA-N"
    print(f"\nTesting with InChIKey: {test_inchikey}")
    
    result = get_compound_info_pubchem(test_inchikey)
    
    if result:
        print(f"\n✓ Compound found:")
        print(f"  Name: {result['name']}")
        print(f"  Formula: {result['formula']}")
        print(f"  PubChem CID: {result['pubchem_cid']}")
        print(f"  Is Carbohydrate: {result['is_carbohydrate']}")
        
        if result['is_carbohydrate']:
            print(f"  Main Class: {result['carbohydrate_main_class']}")
            print(f"  Subclass: {result['carbohydrate_subclass']}")
        
        if result['chebi_ontology']:
            print(f"\n  ChEBI Ontology Terms ({len(result['chebi_ontology'])}):")
            for term in result['chebi_ontology'][:10]:
                print(f"    - {term}")
            if len(result['chebi_ontology']) > 10:
                print(f"    ... and {len(result['chebi_ontology']) - 10} more")
    else:
        print("\n✗ Failed to retrieve compound information")


if __name__ == "__main__":
    main()
