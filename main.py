"""
PubMed Paper Fetcher - A tool to fetch and filter research papers based on author affiliations
"""
import argparse
import csv
import logging
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, Union

import requests
from Bio import Entrez
from rich.console import Console
from rich.logging import RichHandler

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)
logger = logging.getLogger("pubmed_fetcher")

# Set your email for Entrez API
Entrez.email = "your.email@example.com"  # Replace with your email

@dataclass
class Author:
    """Class representing a paper author"""
    name: str
    affiliation: Optional[str] = None
    email: Optional[str] = None
    is_corresponding: bool = False
    is_non_academic: bool = False
    company: Optional[str] = None


@dataclass
class Paper:
    """Class representing a research paper"""
    pubmed_id: str
    title: str
    publication_date: str
    authors: List[Author]

    @property
    def non_academic_authors(self) -> List[Author]:
        """Returns authors affiliated with non-academic institutions"""
        return [author for author in self.authors if author.is_non_academic]

    @property
    def company_affiliations(self) -> Set[str]:
        """Returns a set of unique company affiliations"""
        return {author.company for author in self.non_academic_authors if author.company}

    @property
    def corresponding_author_email(self) -> Optional[str]:
        """Returns the email of the corresponding author if available"""
        for author in self.authors:
            if author.is_corresponding and author.email:
                return author.email
        return None


class PubMedFetcher:
    """Class for fetching and processing papers from PubMed"""

    # List of keywords to identify academic institutions
    ACADEMIC_KEYWORDS = [
        "university", "college", "school", "institute", "academia", "faculty",
        "department", "laboratory", "univ.", "univ", "hospital", "medical center",
        "medical centre", "clinic", "foundation", "institution"
    ]

    # List of keywords to identify pharmaceutical/biotech companies
    COMPANY_KEYWORDS = [
        "pharma", "biotech", "therapeutics", "bioscience", "inc", "corp",
        "co.,", "co.", "ltd", "limited", "llc", "gmbh", "company", "plc",
        "pharmaceutical", "pharmaceuticals", "laboratories", "labs", "ag", "sa"
    ]

    def __init__(self, debug: bool = False):
        """Initialize the PubMed fetcher"""
        self.debug = debug
        if debug:
            logger.setLevel(logging.DEBUG)

    def search_pubmed(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search PubMed for papers matching the query
        
        Args:
            query: PubMed search query
            max_results: Maximum number of results to return
            
        Returns:
            List of PubMed IDs
        """
        logger.debug(f"Searching PubMed with query: {query}")
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            logger.debug(f"Found {len(id_list)} papers")
            return id_list
        except Exception as e:
            logger.error(f"Error searching PubMed: {e}")
            return []

    def fetch_paper_details(self, pubmed_ids: List[str]) -> List[Paper]:
        """
        Fetch detailed information for a list of PubMed IDs
        
        Args:
            pubmed_ids: List of PubMed IDs
            
        Returns:
            List of Paper objects
        """
        if not pubmed_ids:
            return []
            
        logger.debug(f"Fetching details for {len(pubmed_ids)} papers")
        papers = []
        
        try:
            handle = Entrez.efetch(db="pubmed", id=",".join(pubmed_ids), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            for record in records["PubmedArticle"]:
                paper = self._parse_paper_record(record)
                if paper and any(author.is_non_academic for author in paper.authors):
                    papers.append(paper)
            
            logger.debug(f"Processed {len(papers)} papers with non-academic authors")
            return papers
        except Exception as e:
            logger.error(f"Error fetching paper details: {e}")
            return []

    def _parse_paper_record(self, record: Dict) -> Optional[Paper]:
        """
        Parse a PubMed record into a Paper object
        
        Args:
            record: PubMed record from Entrez
            
        Returns:
            Paper object or None if parsing fails
        """
        try:
            article = record["MedlineCitation"]["Article"]
            
            # Extract PubMed ID
            pubmed_id = record["MedlineCitation"]["PMID"]
            
            # Extract title
            title = article["ArticleTitle"]
            
            # Extract publication date
            try:
                date_elements = article["Journal"]["JournalIssue"]["PubDate"]
                pub_date_parts = []
                for key in ["Year", "Month", "Day"]:
                    if key in date_elements:
                        pub_date_parts.append(date_elements[key])
                publication_date = " ".join(pub_date_parts)
            except KeyError:
                publication_date = "Unknown"
            
            # Extract authors
            authors = []
            if "AuthorList" in article:
                for author_data in article["AuthorList"]:
                    if "LastName" not in author_data and "CollectiveName" not in author_data:
                        continue
                    
                    # Get author name
                    if "CollectiveName" in author_data:
                        author_name = author_data["CollectiveName"]
                    else:
                        name_parts = []
                        if "LastName" in author_data:
                            name_parts.append(author_data["LastName"])
                        if "ForeName" in author_data:
                            name_parts.append(author_data["ForeName"])
                        elif "Initials" in author_data:
                            name_parts.append(author_data["Initials"])
                        author_name = " ".join(name_parts)
                    
                    # Get author affiliation and email
                    affiliation = None
                    email = None
                    is_corresponding = False
                    
                    if "AffiliationInfo" in author_data:
                        for affiliation_info in author_data["AffiliationInfo"]:
                            if "Affiliation" in affiliation_info:
                                affiliation = affiliation_info["Affiliation"]
                                # Extract email using regex
                                email_match = re.search(r'[\w\.-]+@[\w\.-]+\.\w+', affiliation)
                                if email_match:
                                    email = email_match.group(0)
                                    is_corresponding = True
                                break
                    
                    # Determine if author is from a non-academic institution
                    is_non_academic, company = self._check_non_academic_affiliation(affiliation)
                    
                    authors.append(Author(
                        name=author_name,
                        affiliation=affiliation,
                        email=email,
                        is_corresponding=is_corresponding,
                        is_non_academic=is_non_academic,
                        company=company
                    ))
            
            return Paper(
                pubmed_id=pubmed_id,
                title=title,
                publication_date=publication_date,
                authors=authors
            )
        except Exception as e:
            logger.error(f"Error parsing paper record: {e}")
            if self.debug:
                logger.exception("Detailed error:")
            return None

    def _check_non_academic_affiliation(self, affiliation: Optional[str]) -> Tuple[bool, Optional[str]]:
        """
        Check if an affiliation is non-academic and extract company name
        
        Args:
            affiliation: Author affiliation string
            
        Returns:
            Tuple of (is_non_academic, company_name)
        """
        if not affiliation:
            return False, None
        
        affiliation_lower = affiliation.lower()
        
        # Check for academic keywords
        for keyword in self.ACADEMIC_KEYWORDS:
            if keyword.lower() in affiliation_lower:
                return False, None
        
        # Check for company keywords
        for keyword in self.COMPANY_KEYWORDS:
            if keyword.lower() in affiliation_lower:
                # Extract company name - this is a simple heuristic
                # Find sentences or comma-separated parts containing the keyword
                parts = re.split(r'[,.;]', affiliation)
                for part in parts:
                    if keyword.lower() in part.lower():
                        company = part.strip()
                        if len(company) > 100:  # Limit length
                            company = company[:100] + "..."
                        return True, company
                
                # If specific company name not found, return the whole affiliation
                if len(affiliation) > 100:  # Limit length
                    return True, affiliation[:100] + "..."
                return True, affiliation
        
        # No academic or company keywords found
        return False, None


def export_to_csv(papers: List[Paper], filename: Optional[str] = None) -> None:
    """
    Export papers to CSV
    
    Args:
        papers: List of Paper objects
        filename: Output filename (if None, print to console)
    """
    if not papers:
        logger.warning("No papers to export")
        return
    
    # Prepare CSV data
    headers = [
        "PubmedID", "Title", "Publication Date", "Non-academic Author(s)",
        "Company Affiliation(s)", "Corresponding Author Email"
    ]
    
    rows = []
    for paper in papers:
        non_academic_authors = "; ".join([author.name for author in paper.non_academic_authors])
        company_affiliations = "; ".join(paper.company_affiliations)
        
        rows.append([
            paper.pubmed_id,
            paper.title,
            paper.publication_date,
            non_academic_authors,
            company_affiliations,
            paper.corresponding_author_email or ""
        ])
    
    # Output to file or console
    if filename:
        try:
            with open(filename, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(rows)
            logger.info(f"Results saved to {filename}")
        except Exception as e:
            logger.error(f"Error writing to CSV file: {e}")
    else:
        # Print to console in a formatted way
        console = Console()
        console.print("\n[bold]Search Results:[/bold]")
        
        for paper in papers:
            console.print(f"\n[bold]Paper ID:[/bold] {paper.pubmed_id}")
            console.print(f"[bold]Title:[/bold] {paper.title}")
            console.print(f"[bold]Publication Date:[/bold] {paper.publication_date}")
            console.print("[bold]Non-academic Authors:[/bold]")
            for author in paper.non_academic_authors:
                console.print(f"  - {author.name}")
                if author.company:
                    console.print(f"    Company: {author.company}")
            if paper.corresponding_author_email:
                console.print(f"[bold]Corresponding Email:[/bold] {paper.corresponding_author_email}")
            console.print("-" * 80)


def main() -> None:
    """Main entry point for the command-line program"""
    parser = argparse.ArgumentParser(
        description="Fetch research papers from PubMed and identify those with authors from pharmaceutical/biotech companies"
    )
    parser.add_argument("query", help="PubMed search query")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")
    parser.add_argument("-f", "--file", help="Output filename for CSV results")
    parser.add_argument("-m", "--max", type=int, default=100, help="Maximum number of results to fetch (default: 100)")
    
    args = parser.parse_args()
    
    # Initialize the fetcher
    fetcher = PubMedFetcher(debug=args.debug)
    
    # Search PubMed
    logger.info(f"Searching PubMed for: {args.query}")
    pubmed_ids = fetcher.search_pubmed(args.query, max_results=args.max)
    
    if not pubmed_ids:
        logger.warning("No papers found matching the query")
        return
    
    # Fetch paper details
    logger.info(f"Fetching details for {len(pubmed_ids)} papers...")
    papers = fetcher.fetch_paper_details(pubmed_ids)
    
    if not papers:
        logger.warning("No papers with non-academic authors found")
        return
    
    logger.info(f"Found {len(papers)} papers with non-academic authors")
    
    # Export results
    export_to_csv(papers, args.file)


if __name__ == "__main__":
    main()