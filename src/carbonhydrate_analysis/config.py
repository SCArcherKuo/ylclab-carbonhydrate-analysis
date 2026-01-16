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

# Retry configuration
MAX_RETRIES = 3  # Maximum number of retries for server errors
RETRY_DELAY = 2.0  # Initial delay between retries (seconds)
RETRY_BACKOFF = 2.0  # Backoff multiplier for exponential backoff

# Server error detection and adaptive sleep
SERVER_ERROR_WINDOW = 60  # Time window in seconds to track server errors
SERVER_ERROR_THRESHOLD = 10  # Number of server errors in window to trigger sleep
SERVER_ERROR_SLEEP_TIME = 180  # Sleep time in seconds (3 minutes) when threshold exceeded

# Failed identifiers storage
FAILED_IDENTIFIERS_DIR = "data/cache"  # Directory to store failed identifiers

# Result limits
MAX_SYNONYMS = 20  # Maximum number of synonyms to return
MAX_ONTOLOGY_DISPLAY = 10  # Maximum ontology terms to display by default
