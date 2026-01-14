"""
Configuration file for carbohydrate classification system.

Contains ChEBI ontology IDs and other constants used in the classification logic.
"""

# ChEBI Ontology IDs
CHEBI_ROOT_ID = 78616  # "carbohydrates and carbohydrate derivatives"
CHEBI_CARBOHYDRATE_ID = 16646  # "carbohydrate"
CHEBI_CARBOHYDRATE_DERIVATIVE_ID = 63299  # "carbohydrate derivative"

# API Configuration
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEBI_API_BASE_URL = "https://www.ebi.ac.uk/chebi/backend/api/public"

# Request timeouts (seconds)
API_TIMEOUT = 30

# Rate limiting (seconds between requests)
RATE_LIMIT_DELAY = 0.2

# Result limits
MAX_SYNONYMS = 20  # Maximum number of synonyms to return
MAX_ONTOLOGY_DISPLAY = 10  # Maximum ontology terms to display by default
