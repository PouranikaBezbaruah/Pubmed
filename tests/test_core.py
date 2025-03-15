# tests/test_core.py

"""
Tests for the core functionality of pubmed_pharma_papers
"""
import pytest

from pubmed_pharma_papers.core import Author, Paper, PubMedFetcher


class TestAuthor:
    """Tests for the Author class"""
    
    def test_author_creation(self) -> None:
        """Test that an Author can be created with the correct attributes"""
        author = Author(
            name="Jane Doe",
            affiliation="Pharma Inc., New York, USA",
            email="jane@pharma.com",
            is_corresponding=True,
            is_non_academic=True,
            company="Pharma Inc."
        )
        
        assert author.name == "Jane Doe"
        assert author.affiliation == "Pharma Inc., New York, USA"
        assert author.email == "jane@pharma.com"
        assert author.is_corresponding is True
        assert author.is_non_academic is True
        assert author.company == "Pharma Inc."


class TestPaper:
    """Tests for the Paper class"""
    
    def test_paper_creation(self) -> None:
        """Test that a Paper can be created with the correct attributes"""
        academic_author = Author(
            name="John Smith",
            affiliation="University of Science, USA",
            is_non_academic=False
        )
        
        pharma_author = Author(
            name="Jane Doe",
            affiliation="Pharma Inc., USA",
            email="jane@pharma.com",
            is_corresponding=True,
            is_non_academic=True,
            company="Pharma Inc."
        )
        
        paper = Paper(
            pubmed_id="12345678",
            title="Test Paper",
            publication_date="2023 Jan",
            authors=[academic_author, pharma_author]
        )
        
        assert paper.pubmed_id == "12345678"
        assert paper.title == "Test Paper"
        assert paper.publication_date == "2023 Jan"
        assert len(paper.authors) == 2
        
    def test_non_academic_authors(self) -> None:
        """Test that non_academic_authors returns only non-academic authors"""
        academic_author = Author(
            name="John Smith",
            affiliation="University of Science, USA",
            is_non_academic=False
        )
        
        pharma_author = Author(
            name="Jane Doe",
            affiliation="Pharma Inc., USA",
            is_non_academic=True,
            company="Pharma Inc."
        )
        
        paper = Paper(
            pubmed_id="12345678",
            title="Test Paper",
            publication_date="2023 Jan",
            authors=[academic_author, pharma_author]
        )
        
        non_academic = paper.non_academic_authors
        assert len(non_academic) == 1
        assert non_academic[0].name == "Jane Doe"
        
    def test_company_affiliations(self) -> None:
        """Test that company_affiliations returns unique company names"""
        author1 = Author(
            name="Jane Doe",
            is_non_academic=True,
            company="Pharma Inc."
        )
        
        author2 = Author(
            name="Bob Johnson",
            is_non_academic=True,
            company="Biotech Ltd."
        )
        
        author3 = Author(
            name="Alice Brown",
            is_non_academic=True,
            company="Pharma Inc."  # Duplicate
        )
        
        paper = Paper(
            pubmed_id="12345678",
            title="Test Paper",
            publication_date="2023 Jan",
            authors=[author1, author2, author3]
        )
        
        companies = paper.company_affiliations
        assert len(companies) == 2
        assert "Pharma Inc." in companies
        assert "Biotech Ltd." in companies
        
    def test_corresponding_author_email(self) -> None:
        """Test that corresponding_author_email returns the correct email"""
        author1 = Author(
            name="Jane Doe",
            email="jane@example.com",
            is_corresponding=False
        )
        
        author2 = Author(
            name="Bob Johnson",
            email="bob@example.com",
            is_corresponding=True
        )
        
        paper = Paper(
            pubmed_id="12345678",
            title="Test Paper",
            publication_date="2023 Jan",
            authors=[author1, author2]
        )
        
        assert paper.corresponding_author_email == "bob@example.com"


class TestPubMedFetcher:
    """Tests for the PubMedFetcher class"""
    
    def test_check_non_academic_affiliation(self) -> None:
        """Test that _check_non_academic_affiliation correctly identifies affiliations"""
        fetcher = PubMedFetcher()
        
        # Academic affiliations
        is_non_academic, company = fetcher._check_non_academic_affiliation("Department of Medicine, University of Science, USA")
        assert is_non_academic is False
        assert company is None
        
        is_non_academic, company = fetcher._check_non_academic_affiliation("Hospital General, Medical Center, Spain")
        assert is_non_academic is False
        assert company is None
        
        # Pharmaceutical/biotech affiliations
        is_non_academic, company = fetcher._check_non_academic_affiliation("Research Division, Pharma Inc., New York, USA")
        assert is_non_academic is True
        assert "Pharma Inc." in company
        
        is_non_academic, company = fetcher._check_non_academic_affiliation("Biotech Corporation, San Francisco, CA, USA")
        assert is_non_academic is True
        assert "Biotech Corporation" in company
        
        # None affiliation
        is_non_academic, company = fetcher._check_non_academic_affiliation(None)
        assert is_non_academic is False
        assert company is None