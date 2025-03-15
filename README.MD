# PubMed Pharma Papers

A Python tool to fetch research papers from PubMed and identify those with authors affiliated with pharmaceutical or biotech companies.

## Features

- Search PubMed using its full query syntax
- Identify authors affiliated with pharmaceutical/biotech companies
- Export results to CSV or display in console
- Command-line interface with various options
- Modular design that can be used as a library

## Installation

This project uses [Poetry](https://python-poetry.org/) for dependency management.

```bash
# Clone the repository
git clone https://github.com/yourusername/pubmed-pharma-papers.git
cd pubmed-pharma-papers

# Install dependencies
poetry install

# Activate the virtual environment
poetry shell
```

## Usage

### Command-line Interface

After installation, you can use the `get-papers-list` command:

```bash
# Basic usage
get-papers-list "cancer therapy"

# Save results to a file
get-papers-list "diabetes" -f results.csv

# Enable debug output
get-papers-list "alzheimer's" -d

# Increase maximum results
get-papers-list "obesity" -m 200

# Get help
get-papers-list --help
```

### Query Syntax

The tool supports PubMed's full query syntax. Some examples:

- Simple search: `"diabetes"`
- Author search: `"Smith J[Author]"`
- Date range: `"2020[pdat]:2023[pdat]"`
- Combined search: `"cancer AND (therapy OR treatment) AND 2022[pdat]"`

For more information about PubMed's search syntax, visit the [PubMed Help page](https://pubmed.ncbi.nlm.nih.gov/help/).

### Using as a Library

You can also use this package as a library in your Python code:

```python
from pubmed_pharma_papers import PubMedFetcher

# Initialize the fetcher
fetcher = PubMedFetcher(debug=True)

# Search PubMed
pubmed_ids = fetcher.search_pubmed("diabetes", max_results=50)

# Fetch paper details
papers = fetcher.fetch_paper_details(pubmed_ids)

# Process the papers
for paper in papers:
    print(f"Paper: {paper.title}")
    print(f"Non-academic authors: {[a.name for a in paper.non_academic_authors]}")
    print(f"Companies: {paper.company_affiliations}")
```

## Code Organization

The project is organized as follows:

```
pubmed_pharma_papers/
├── __init__.py         # Package exports
├── core.py             # Core functionality (models and fetcher)
└── cli.py              # Command-line interface
```

- `core.py` - Contains the main functionality:
  - `Author` class - Represents a paper author with affiliation info
  - `Paper` class - Represents a research paper with its metadata
  - `PubMedFetcher` class - Handles PubMed API interactions and processing

- `cli.py` - Implements the command-line interface:
  - Argument parsing
  - CSV export functionality
  - Console output formatting

## How It Works

1. The tool searches PubMed using the Entrez API from Biopython
2. It fetches detailed information for each paper in the search results
3. For each paper, it analyzes author affiliations to identify non-academic authors
4. It uses keywords and heuristics to identify pharmaceutical/biotech companies
5. Results are compiled and exported in the requested format

## Development

### Testing

```bash
# Run tests
poetry run pytest
```

### Type Checking

```bash
# Run mypy
poetry run mypy pubmed_pharma_papers
```

### Formatting

```bash
# Format code with black
poetry run black .

# Sort imports
poetry run isort .
```

## Tools and Libraries Used

- [Poetry](https://python-poetry.org/) - Dependency management and packaging
- [Biopython](https://biopython.org/) - Entrez API access for PubMed queries
- [Rich](https://rich.readthedocs.io/) - Terminal formatting and improved logging
- [Pytest](https://docs.pytest.org/) - Testing
- [Mypy](https://mypy.readthedocs.io/) - Static type checking
- [Black](https://black.readthedocs.io/) - Code formatting
- [iSort](https://pycqa.github.io/isort/) - Import sorting

## Publishing to Test PyPI

To publish the package to Test PyPI:

```bash
# Build the package
poetry build

# Configure Test PyPI repository
poetry config repositories.testpypi https://test.pypi.org/legacy/

# Publish to Test PyPI
poetry publish -r testpypi
```

## License

MIT