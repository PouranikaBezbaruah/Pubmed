# pubmed_pharma_papers/cli.py

"""
Command-line interface for the PubMed pharma papers tool
"""
import argparse
import csv
import logging
import sys
from typing import List, Optional

from rich.console import Console
from rich.logging import RichHandler

from pubmed_pharma_papers.core import Paper, PubMedFetcher

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)
logger = logging.getLogger("pubmed_pharma_papers")


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