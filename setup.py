from setuptools import setup, find_packages

setup(
    name="pubmed-pharma-papers",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "requests",
        "rich",
    ],
    entry_points={
        "console_scripts": [
            "get-papers-list=pubmed_pharma_papers.cli:main",
        ],
    },
)