# pubmed_pharma_papers/__init__.py

"""
PubMed Pharma Papers - A tool to fetch research papers from PubMed with authors affiliated with pharmaceutical companies
"""

__version__ = "0.1.0"

from pubmed_pharma_papers.core import Author, Paper, PubMedFetcher

__all__ = ["Author", "Paper", "PubMedFetcher"]