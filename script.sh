#!/bin/bash
# Script to publish the package to Test PyPI

echo "Building the package..."
poetry build

echo "Configuring Test PyPI repository..."
poetry config repositories.testpypi https://test.pypi.org/legacy/

echo "Publishing to Test PyPI..."
poetry publish -r testpypi

echo "Package published successfully!"
echo "You can now install it with:"
echo "pip install --index-url https://test.pypi.org/simple/ pubmed-pharma-papers"