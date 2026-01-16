"""
ChEBI API Client

This module provides functions for interacting with the ChEBI (Chemical Entities of
Biological Interest) API to retrieve ontology information about chemical compounds.
"""

import requests
import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from collections import OrderedDict
from loguru import logger
from . import config


class LRUCache:
    """
    Least Recently Used (LRU) cache with maximum size limit.
    
    Automatically evicts least recently used items when cache is full.
    This prevents unbounded memory growth during long-running processes.
    """
    
    def __init__(self, max_size: int = 200):
        """
        Initialize LRU cache.
        
        Parameters:
        -----------
        max_size : int
            Maximum number of items to keep in cache (default: 200)
        """
        self.cache = OrderedDict()
        self.max_size = max_size
    
    def get(self, key, default=None):
        """Get value from cache, moving it to end (most recently used)."""
        if key in self.cache:
            # Move to end (most recent)
            self.cache.move_to_end(key)
            return self.cache[key]
        return default
    
    def __getitem__(self, key):
        """Dictionary-style access."""
        value = self.get(key)
        if value is None:
            raise KeyError(key)
        return value
    
    def __contains__(self, key):
        """Check if key exists in cache."""
        return key in self.cache
    
    def set(self, key, value):
        """Set value in cache, evicting oldest item if cache is full."""
        if key in self.cache:
            # Update existing and move to end
            self.cache.move_to_end(key)
        self.cache[key] = value
        
        # Evict oldest item if over limit
        if len(self.cache) > self.max_size:
            self.cache.popitem(last=False)  # Remove first (oldest) item
    
    def __setitem__(self, key, value):
        """Dictionary-style assignment."""
        self.set(key, value)
    
    def clear(self):
        """Clear all items from cache."""
        self.cache.clear()
    
    def __len__(self):
        """Return number of items in cache."""
        return len(self.cache)
    
    def items(self):
        """Return cache items."""
        return self.cache.items()
    
    def keys(self):
        """Return cache keys."""
        return self.cache.keys()
    
    def values(self):
        """Return cache values."""
        return self.cache.values()


# Cache for ChEBI API responses to avoid repeated calls
# Using LRU cache to prevent unbounded memory growth
_chebi_children_cache: LRUCache = LRUCache(max_size=200)
_chebi_parents_cache: LRUCache = LRUCache(max_size=200)
_chebi_ancestors_cache: LRUCache = LRUCache(max_size=200)


class ChEBIClient:
    """
    Client for interacting with ChEBI API.
    
    This class provides methods to query the ChEBI ontology database
    and retrieve hierarchical relationships between chemical entities.
    Supports persistent caching to reduce API calls.
    """
    
    def __init__(self, base_url: Optional[str] = None, timeout: Optional[int] = None, 
                 cache_dir: Optional[Path] = None, use_cache: bool = True):
        """
        Initialize ChEBI API client.
        
        Parameters:
        -----------
        base_url : str, optional
            Base URL for ChEBI API (defaults to config.CHEBI_API_BASE_URL)
        timeout : int, optional
            Request timeout in seconds (defaults to config.API_TIMEOUT)
        cache_dir : Path, optional
            Directory for cache storage (defaults to data/cache/)
        use_cache : bool, optional
            Whether to use persistent caching (default: True)
        """
        self.base_url = base_url or config.CHEBI_API_BASE_URL
        self.timeout = timeout or config.API_TIMEOUT
        self.use_cache = use_cache
        
        # Set up cache directory
        if cache_dir:
            self.cache_dir = cache_dir
        else:
            # Default to data/cache relative to project root
            project_root = Path(__file__).parent.parent.parent
            self.cache_dir = project_root / 'data' / 'cache'
        
        self.cache_file = self.cache_dir / 'chebi_children_cache.json'
        self.parents_cache_file = self.cache_dir / 'chebi_parents_cache.json'
        self.ancestors_cache_file = self.cache_dir / 'chebi_ancestors_cache.json'
        
        # Use global LRU caches for backward compatibility
        self.cache = _chebi_children_cache
        self.parents_cache = _chebi_parents_cache
        self.ancestors_cache = _chebi_ancestors_cache
        
        # Track if cache has been modified
        self._cache_dirty = False
        self._cache_save_counter = 0
        self._save_batch_size = 10  # Save after every N new entries
        
        # Load persistent cache if enabled
        if self.use_cache:
            self._load_cache()
            self._load_parents_cache()
            self._load_ancestors_cache()
    
    def get_children(self, chebi_id: int) -> List[Dict[str, Any]]:
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
        
        Examples:
        ---------
        >>> client = ChEBIClient()
        >>> children = client.get_children(16646)  # carbohydrate
        >>> print(len(children))
        """
        # Check cache first
        if chebi_id in self.cache:
            logger.debug(f"Cache hit for ChEBI children: {chebi_id}")
            return self.cache[chebi_id]
        
        logger.debug(f"Fetching ChEBI children for {chebi_id}")
        try:
            url = f'{self.base_url}/ontology/children/{chebi_id}/'
            response = requests.get(url, headers={'accept': '*/*'}, timeout=self.timeout)
            
            if response.status_code != 200:
                logger.warning(f"ChEBI API failed for {chebi_id}: HTTP {response.status_code}")
                # Don't cache failures - allow retry on next call
                return []
            
            data = response.json()
            incoming = data.get('ontology_relations', {}).get('incoming_relations', [])
            
            # Filter for "is a" relations (direct children)
            children = [rel for rel in incoming if rel.get('relation_type') == 'is a']
            logger.debug(f"Found {len(children)} children for ChEBI {chebi_id}")
            
            self.cache[chebi_id] = children
            self._cache_dirty = True
            self._cache_save_counter += 1
            
            # Save cache periodically (batch saves)
            if self.use_cache and self._cache_save_counter >= self._save_batch_size:
                self._save_cache()
                self._cache_save_counter = 0
            
            return children
            
        except requests.exceptions.Timeout as e:
            logger.error(f"Timeout fetching ChEBI children for {chebi_id}: {str(e)}")
            # Don't cache errors - allow retry
            return []
        except requests.exceptions.RequestException as e:
            logger.error(f"Request error fetching ChEBI children for {chebi_id}: {str(e)}")
            # Don't cache errors - allow retry
            return []
        except Exception as e:
            logger.exception(f"Unexpected error fetching ChEBI children for {chebi_id}: {str(e)}")
            # Don't cache errors - allow retry
            return []
    
    def get_parents(self, chebi_id: int) -> List[Dict[str, Any]]:
        """
        Get direct parents of a ChEBI ID from ChEBI backend API.
        
        Parameters:
        -----------
        chebi_id : int
            ChEBI ID (without 'CHEBI:' prefix)
        
        Returns:
        --------
        list of dict
            List of direct parents with 'is a' relationship.
            Each dict contains: 'init_id', 'init_name', 'relation_type', 'final_id', 'final_name'
        
        Examples:
        ---------
        >>> client = ChEBIClient()
        >>> parents = client.get_parents(4167)  # D-glucopyranose
        >>> print(parents)
        """
        # Check cache first
        if chebi_id in self.parents_cache:
            logger.debug(f"Cache hit for ChEBI parents: {chebi_id}")
            return self.parents_cache[chebi_id]
        
        logger.debug(f"Fetching ChEBI parents for {chebi_id}")
        try:
            url = f'{self.base_url}/ontology/parents/{chebi_id}/'
            response = requests.get(url, headers={'accept': '*/*'}, timeout=self.timeout)
            
            if response.status_code != 200:
                logger.warning(f"ChEBI API failed for {chebi_id}: HTTP {response.status_code}")
                # Don't cache failures - allow retry on next call
                return []
            
            data = response.json()
            outgoing = data.get('ontology_relations', {}).get('outgoing_relations', [])
            
            # Filter for "is a" relations (direct parents)
            parents = [rel for rel in outgoing if rel.get('relation_type') == 'is a']
            logger.debug(f"Found {len(parents)} parents for ChEBI {chebi_id}")
            
            self.parents_cache[chebi_id] = parents
            self._cache_dirty = True
            self._cache_save_counter += 1
            
            # Save cache periodically (batch saves)
            if self.use_cache and self._cache_save_counter >= self._save_batch_size:
                self._save_cache()
                self._save_parents_cache()
                self._cache_save_counter = 0
            
            return parents
            
        except requests.exceptions.Timeout as e:
            logger.error(f"Timeout fetching ChEBI parents for {chebi_id}: {str(e)}")
            # Don't cache errors - allow retry
            return []
        except requests.exceptions.RequestException as e:
            logger.error(f"Request error fetching ChEBI parents for {chebi_id}: {str(e)}")
            # Don't cache errors - allow retry
            return []
        except Exception as e:
            logger.exception(f"Unexpected error fetching ChEBI parents for {chebi_id}: {str(e)}")
            # Don't cache errors - allow retry
            return []
    
    def get_all_ancestors(self, chebi_id: int, max_depth: int = 20) -> List[int]:
        """
        Get all ancestor ChEBI IDs by recursively traversing parent relationships.
        
        Parameters:
        -----------
        chebi_id : int
            Starting ChEBI ID (without 'CHEBI:' prefix)
        max_depth : int, optional
            Maximum recursion depth to prevent infinite loops (default: 20)
        
        Returns:
        --------
        list of int
            List of all ancestor ChEBI IDs (includes the starting ID)
        
        Examples:
        ---------
        >>> client = ChEBIClient()
        >>> ancestors = client.get_all_ancestors(4167)  # D-glucopyranose
        >>> print(16646 in ancestors)  # Check if carbohydrate ancestor
        True
        """
        # Check cache first
        if chebi_id in self.ancestors_cache:
            logger.debug(f"Cache hit for ChEBI ancestors: {chebi_id}")
            return self.ancestors_cache[chebi_id]
        
        logger.debug(f"Computing all ancestors for ChEBI {chebi_id}")
        ancestors = set([chebi_id])  # Include self
        visited = set()
        queue = [chebi_id]
        depth = 0
        
        try:
            while queue and depth < max_depth:
                current_batch = queue[:]
                queue = []
                
                for current_id in current_batch:
                    if current_id in visited:
                        continue
                    visited.add(current_id)
                    
                    parents = self.get_parents(current_id)
                    for parent_rel in parents:
                        parent_id = parent_rel.get('final_id')
                        if parent_id and parent_id not in ancestors:
                            ancestors.add(parent_id)
                            queue.append(parent_id)
                
                depth += 1
            
            ancestor_list = sorted(list(ancestors))
            logger.debug(f"Found {len(ancestor_list)} ancestors for ChEBI {chebi_id}")
            self.ancestors_cache[chebi_id] = ancestor_list
            self._cache_dirty = True
            
            return ancestor_list
            
        except Exception as e:
            logger.exception(f"Error computing ancestors for ChEBI {chebi_id}: {str(e)}")
            # Return partial results if available
            ancestor_list = sorted(list(ancestors))
            return ancestor_list
    
    def get_main_groups(self, chebi_id: int) -> List[str]:
        """
        Get main groups (children with >1 children) of a ChEBI ID.
        
        Main groups are defined as direct children that themselves have
        more than one child, indicating significant branching in the ontology.
        
        Parameters:
        -----------
        chebi_id : int
            ChEBI ID to get main groups for
        
        Returns:
        --------
        list of str
            List of main group names (children that have >1 children themselves)
        
        Examples:
        ---------
        >>> client = ChEBIClient()
        >>> main_groups = client.get_main_groups(16646)  # carbohydrate
        >>> print(main_groups)
        """
        children = self.get_children(chebi_id)
        main_groups = []
        
        for child in children:
            child_id = child.get('init_id')
            child_name = child.get('init_name', '')
            
            if child_id:
                # Check if this child has more than 1 children
                grandchildren = self.get_children(child_id)
                if len(grandchildren) > 1:
                    main_groups.append(child_name)
        
        return main_groups
    
    def clear_cache(self) -> None:
        """Clear the ChEBI children cache."""
        self.cache.clear()
    
    def _load_cache(self) -> None:
        """
        Load cache from disk if it exists.
        
        The cache is stored as JSON with ChEBI IDs as string keys
        (JSON doesn't support integer keys).
        """
        if not self.cache_file.exists():
            return
        
        try:
            with open(self.cache_file, 'r', encoding='utf-8') as f:
                cache_data = json.load(f)
            
            # Convert string keys back to integers and load into LRU cache
            # Only load non-empty results (don't load cached failures)
            for key_str, value in cache_data.items():
                if value:  # Skip empty lists
                    self.cache.set(int(key_str), value)
            
            print(f"Loaded ChEBI children cache with {len(self.cache)} entries")
        except Exception as e:
            print(f"Warning: Failed to load ChEBI children cache: {str(e)}")
    
    def _load_parents_cache(self) -> None:
        """Load parents cache from disk if it exists."""
        if not self.parents_cache_file.exists():
            return
        
        try:
            with open(self.parents_cache_file, 'r', encoding='utf-8') as f:
                cache_data = json.load(f)
            
            # Convert string keys back to integers and load into LRU cache
            # Only load non-empty results (don't load cached failures)
            for key_str, value in cache_data.items():
                if value:  # Skip empty lists
                    self.parents_cache.set(int(key_str), value)
            
            print(f"Loaded ChEBI parents cache with {len(self.parents_cache)} entries")
        except Exception as e:
            print(f"Warning: Failed to load ChEBI parents cache: {str(e)}")
    
    def _load_ancestors_cache(self) -> None:
        """Load ancestors cache from disk if it exists."""
        if not self.ancestors_cache_file.exists():
            return
        
        try:
            with open(self.ancestors_cache_file, 'r', encoding='utf-8') as f:
                cache_data = json.load(f)
            
            # Convert string keys back to integers and load into LRU cache
            # Only load non-empty results (don't load cached failures)
            for key_str, value in cache_data.items():
                if value:  # Skip empty lists
                    self.ancestors_cache.set(int(key_str), value)
            
            print(f"Loaded ChEBI ancestors cache with {len(self.ancestors_cache)} entries")
        except Exception as e:
            print(f"Warning: Failed to load ChEBI ancestors cache: {str(e)}")
    
    def _save_cache(self) -> None:
        """
        Save cache to disk in JSON format.
        
        Creates the cache directory if it doesn't exist.
        Converts integer keys to strings for JSON compatibility.
        """
        if not self.use_cache or not self._cache_dirty:
            return
        
        try:
            # Create cache directory if it doesn't exist
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            
            # Convert integer keys to strings for JSON
            # Only save non-empty results (don't persist failures)
            cache_data = {str(key): value for key, value in self.cache.items() if value}
            
            # Write to temporary file first, then rename for atomicity
            temp_file = self.cache_file.with_suffix('.tmp')
            with open(temp_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, indent=2, ensure_ascii=False)
            
            # Atomic rename
            temp_file.replace(self.cache_file)
            
            self._cache_dirty = False
            print(f"Saved ChEBI children cache with {len(self.cache)} entries")
        except Exception as e:
            print(f"Warning: Failed to save ChEBI children cache: {str(e)}")
    
    def _save_parents_cache(self) -> None:
        """Save parents cache to disk in JSON format."""
        if not self.use_cache or not self._cache_dirty:
            return
        
        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            # Only save non-empty results (don't persist failures)
            cache_data = {str(key): value for key, value in self.parents_cache.items() if value}
            
            temp_file = self.parents_cache_file.with_suffix('.tmp')
            with open(temp_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, indent=2, ensure_ascii=False)
            
            temp_file.replace(self.parents_cache_file)
            print(f"Saved ChEBI parents cache with {len(self.parents_cache)} entries")
        except Exception as e:
            print(f"Warning: Failed to save ChEBI parents cache: {str(e)}")
    
    def _save_ancestors_cache(self) -> None:
        """Save ancestors cache to disk in JSON format."""
        if not self.use_cache or not self._cache_dirty:
            return
        
        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            # Only save non-empty results (don't persist failures)
            cache_data = {str(key): value for key, value in self.ancestors_cache.items() if value}
            
            temp_file = self.ancestors_cache_file.with_suffix('.tmp')
            with open(temp_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, indent=2, ensure_ascii=False)
            
            temp_file.replace(self.ancestors_cache_file)
            print(f"Saved ChEBI ancestors cache with {len(self.ancestors_cache)} entries")
        except Exception as e:
            print(f"Warning: Failed to save ChEBI ancestors cache: {str(e)}")
    
    def save_cache(self) -> None:
        """
        Manually save all caches to disk.
        
        This is called automatically in batches and on cleanup,
        but can be called manually to ensure data is persisted.
        """
        self._save_cache()
        self._save_parents_cache()
        self._save_ancestors_cache()
    
    def __del__(self):
        """
        Destructor to save all caches when object is garbage collected.
        
        This ensures that any pending cache updates are saved even
        if save_cache() wasn't called explicitly.
        """
        if self.use_cache and self._cache_dirty:
            self._save_cache()
            self._save_parents_cache()
            self._save_ancestors_cache()


# =============================================================================
# Module-level convenience functions for backward compatibility
# =============================================================================

# Global client instance
_default_client = ChEBIClient()


def get_chebi_children(chebi_id: int) -> List[Dict[str, Any]]:
    """
    Get direct children of a ChEBI ID from ChEBI backend API.
    
    Convenience function that uses the default ChEBI client.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID (without 'CHEBI:' prefix)
    
    Returns:
    --------
    list of dict
        List of direct children with 'is a' relationship.
    """
    return _default_client.get_children(chebi_id)


def get_main_groups(chebi_id: int) -> List[str]:
    """
    Get main groups (children with >1 children) of a ChEBI ID.
    
    Convenience function that uses the default ChEBI client.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID to get main groups for
    
    Returns:
    --------
    list of str
        List of main group names
    """
    return _default_client.get_main_groups(chebi_id)


def get_chebi_parents(chebi_id: int) -> List[Dict[str, Any]]:
    """
    Get direct parents of a ChEBI ID from ChEBI backend API.
    
    Convenience function that uses the default ChEBI client.
    
    Parameters:
    -----------
    chebi_id : int
        ChEBI ID (without 'CHEBI:' prefix)
    
    Returns:
    --------
    list of dict
        List of direct parents with 'is a' relationship.
    """
    return _default_client.get_parents(chebi_id)


def get_all_ancestors(chebi_id: int) -> List[int]:
    """
    Get all ancestor ChEBI IDs by recursively traversing parent relationships.
    
    Convenience function that uses the default ChEBI client.
    
    Parameters:
    -----------
    chebi_id : int
        Starting ChEBI ID (without 'CHEBI:' prefix)
    
    Returns:
    --------
    list of int
        List of all ancestor ChEBI IDs (includes the starting ID)
    """
    return _default_client.get_all_ancestors(chebi_id)
