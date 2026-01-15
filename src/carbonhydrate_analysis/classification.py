"""
Carbohydrate Classification

This module provides functions to classify compounds as carbohydrates based on
their ChEBI ontology terms and hierarchical relationships.
"""

from typing import List, Any, Tuple, Optional
from loguru import logger
from .chebi_api import get_chebi_children, get_main_groups, get_all_ancestors
from .utils import extract_term_string
from . import config


class CarbohydrateClassifier:
    """
    Classifier for identifying and categorizing carbohydrates.
    
    This class implements the classification logic for determining whether
    a compound is a carbohydrate and categorizing it into one of five main
    categories based on ChEBI ontology.
    """
    
    def __init__(self):
        """Initialize the carbohydrate classifier."""
        self.root_id = config.CHEBI_ROOT_ID
        self.carbohydrate_id = config.CHEBI_CARBOHYDRATE_ID
        self.carbohydrate_derivative_id = config.CHEBI_CARBOHYDRATE_DERIVATIVE_ID
    
    def classify(self, chebi_ontology: List[Any]) -> Tuple[Optional[str], Optional[str]]:
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
        
        Examples:
        ---------
        >>> classifier = CarbohydrateClassifier()
        >>> ontology = ['carbohydrates and carbohydrate derivatives', 'carbohydrate', 'monosaccharide']
        >>> main_class, subclass = classifier.classify(ontology)
        >>> print(main_class, subclass)
        'main carbohydrate group' 'monosaccharide'
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
        carb_main_groups = get_main_groups(self.carbohydrate_id)
        carb_deriv_main_groups = get_main_groups(self.carbohydrate_derivative_id)
        
        carb_main_groups_lower = [g.lower() for g in carb_main_groups]
        carb_deriv_main_groups_lower = [g.lower() for g in carb_deriv_main_groups]
        
        # Check for carbohydrate (CHEBI:16646)
        if 'carbohydrate' in clean_terms_lower:
            return self._classify_carbohydrate(
                clean_terms, clean_terms_lower,
                carb_main_groups, carb_main_groups_lower
            )
        
        # Check for carbohydrate derivative (CHEBI:63299)
        if 'carbohydrate derivative' in clean_terms_lower:
            return self._classify_carbohydrate_derivative(
                clean_terms, clean_terms_lower,
                carb_deriv_main_groups, carb_deriv_main_groups_lower
            )
        
        # If under root but not in carbohydrate or carbohydrate derivative
        return self._classify_other(clean_terms, clean_terms_lower)
    
    def _classify_carbohydrate(
        self,
        clean_terms: List[str],
        clean_terms_lower: List[str],
        main_groups: List[str],
        main_groups_lower: List[str]
    ) -> Tuple[str, str]:
        """Classify carbohydrate into main or other group."""
        # Check if matches main carbohydrate groups
        for i, term_lower in enumerate(clean_terms_lower):
            if term_lower in main_groups_lower:
                idx = main_groups_lower.index(term_lower)
                return 'main carbohydrate group', main_groups[idx]
        
        # If no main group match, it's "other carbohydrate"
        # Get all direct children and find which one matches
        carb_children = get_chebi_children(self.carbohydrate_id)
        carb_children_names = [child['init_name'].lower() for child in carb_children]
        
        # Find the direct child that appears in the ontology terms
        for term in clean_terms:
            term_lower = term.lower()
            if term_lower in carb_children_names:
                # Return the original capitalization from the ontology
                return 'other carbohydrate', term
        
        return 'other carbohydrate', 'carbohydrate'
    
    def _classify_carbohydrate_derivative(
        self,
        clean_terms: List[str],
        clean_terms_lower: List[str],
        main_groups: List[str],
        main_groups_lower: List[str]
    ) -> Tuple[str, str]:
        """Classify carbohydrate derivative into main or other group."""
        # Check if matches main carbohydrate derivative groups
        for i, term_lower in enumerate(clean_terms_lower):
            if term_lower in main_groups_lower:
                idx = main_groups_lower.index(term_lower)
                return 'main carbohydrate derivative group', main_groups[idx]
        
        # If no main group match, it's "other carbohydrate derivative"
        # Get all direct children and find which one matches
        carb_deriv_children = get_chebi_children(self.carbohydrate_derivative_id)
        carb_deriv_children_names = [child['init_name'].lower() for child in carb_deriv_children]
        
        # Find the direct child that appears in the ontology terms
        for term in clean_terms:
            term_lower = term.lower()
            if term_lower in carb_deriv_children_names:
                # Return the original capitalization from the ontology
                return 'other carbohydrate derivative', term
        
        return 'other carbohydrate derivative', 'carbohydrate derivative'
    
    def _classify_other(
        self,
        clean_terms: List[str],
        clean_terms_lower: List[str]
    ) -> Tuple[str, str]:
        """Classify as 'other' category."""
        # Get all direct children and find which one matches
        root_children = get_chebi_children(self.root_id)
        root_children_names = [child['init_name'].lower() for child in root_children]
        
        # Find the direct child that appears in the ontology terms
        for term in clean_terms:
            term_lower = term.lower()
            if term_lower in root_children_names and term_lower != 'carbohydrates and carbohydrate derivatives':
                # Return the original capitalization from the ontology
                return 'other', term
        
        return 'other', 'carbohydrates and carbohydrate derivatives'


# =============================================================================
# Module-level convenience function for backward compatibility
# =============================================================================

# Global classifier instance
_default_classifier = CarbohydrateClassifier()


def classify_carbohydrate(chebi_ontology: List[Any]) -> Tuple[Optional[str], Optional[str]]:
    """
    Classify a compound based on its ChEBI ontology terms.
    
    Convenience function that uses the default CarbohydrateClassifier.
    
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
    return _default_classifier.classify(chebi_ontology)


def classify_by_chebi_ancestry(chebi_id: int) -> Tuple[Optional[str], Optional[str]]:
    """
    Classify a compound based on its ChEBI ID ancestry using direct ChEBI API calls.
    
    This is more efficient than the PubChem classification method as it directly
    queries ChEBI and walks up the ontology tree to determine classification.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID (without 'CHEBI:' prefix)
    
    Returns:
    --------
    tuple of (main_class, subclass)
        main_class: One of ['main carbohydrate group', 'other carbohydrate', 
                           'main carbohydrate derivative group', 'other carbohydrate derivative', 
                           'other', None]
        subclass: Specific subclass name or None
    
    Examples:
    ---------
    >>> main_class, subclass = classify_by_chebi_ancestry(15365)  # D-glucose
    >>> print(main_class, subclass)
    'main carbohydrate group' 'monosaccharide'
    """
    logger.debug(f"Classifying ChEBI {chebi_id} by ancestry")
    
    try:
        # Get all ancestors
        ancestors = get_all_ancestors(chebi_id)
        logger.debug(f"ChEBI {chebi_id} has {len(ancestors)} ancestors")
        
        # Check if compound belongs to carbohydrates root (CHEBI:78616)
        if config.CHEBI_ROOT_ID not in ancestors:
            logger.debug(f"ChEBI {chebi_id} not under carbohydrate root")
            return None, None
        
        logger.debug(f"ChEBI {chebi_id} is under carbohydrate root")
        
        # Get main groups for carbohydrate and carbohydrate derivative
        carb_main_groups = get_main_groups(config.CHEBI_CARBOHYDRATE_ID)
        carb_deriv_main_groups = get_main_groups(config.CHEBI_CARBOHYDRATE_DERIVATIVE_ID)
        
        # Check if it's a carbohydrate (CHEBI:16646)
        if config.CHEBI_CARBOHYDRATE_ID in ancestors:
            logger.debug(f"ChEBI {chebi_id} is a carbohydrate")
            return _classify_by_ancestry_path(
                chebi_id, 
                ancestors, 
                config.CHEBI_CARBOHYDRATE_ID,
                carb_main_groups,
                'main carbohydrate group',
                'other carbohydrate'
            )
        
        # Check if it's a carbohydrate derivative (CHEBI:63299)
        if config.CHEBI_CARBOHYDRATE_DERIVATIVE_ID in ancestors:
            logger.debug(f"ChEBI {chebi_id} is a carbohydrate derivative")
            return _classify_by_ancestry_path(
                chebi_id,
                ancestors,
                config.CHEBI_CARBOHYDRATE_DERIVATIVE_ID,
                carb_deriv_main_groups,
                'main carbohydrate derivative group',
                'other carbohydrate derivative'
            )
        
        # If under root but not in carbohydrate or carbohydrate derivative
        logger.debug(f"ChEBI {chebi_id} under root but not carbohydrate or derivative")
        root_children = get_chebi_children(config.CHEBI_ROOT_ID)
        for child in root_children:
            child_id = child.get('init_id')
            if child_id and child_id in ancestors:
                result = ('other', child.get('init_name', 'carbohydrates and carbohydrate derivatives'))
                logger.debug(f"ChEBI {chebi_id} classified as: {result}")
                return result
        
        result = ('other', 'carbohydrates and carbohydrate derivatives')
        logger.debug(f"ChEBI {chebi_id} classified as: {result}")
        return result
        
    except Exception as e:
        logger.exception(f"Error classifying ChEBI {chebi_id} by ancestry: {str(e)}")
        raise
    
    return 'other', 'carbohydrates and carbohydrate derivatives'


def _classify_by_ancestry_path(
    chebi_id: int,
    ancestors: List[int],
    parent_category_id: int,
    main_groups: List[str],
    main_class_label: str,
    other_class_label: str
) -> Tuple[str, str]:
    """
    Helper function to classify based on ancestry path.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID being classified
    ancestors : list of int
        List of all ancestor IDs
    parent_category_id : int
        ID of the parent category (carbohydrate or carbohydrate derivative)
    main_groups : list of str
        List of main group names
    main_class_label : str
        Label for main class matches
    other_class_label : str
        Label for other category matches
    
    Returns:
    --------
    tuple of (main_class, subclass)
    """
    # Get children of the parent category
    parent_children = get_chebi_children(parent_category_id)
    
    # Build mapping of child_id to child_name
    child_id_to_name = {child['init_id']: child['init_name'] for child in parent_children}
    
    # Find which direct child of parent_category is in the ancestry
    for ancestor_id in ancestors:
        if ancestor_id in child_id_to_name:
            child_name = child_id_to_name[ancestor_id]
            
            # Check if this is a main group
            if child_name in main_groups:
                return main_class_label, child_name
            else:
                return other_class_label, child_name
    
    # Fallback to parent category name
    parent_name = 'carbohydrate' if parent_category_id == config.CHEBI_CARBOHYDRATE_ID else 'carbohydrate derivative'
    return other_class_label, parent_name

