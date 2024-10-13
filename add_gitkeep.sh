#!/bin/bash

# Find all empty directories and add a .gitkeep file
find . -type d -empty -exec touch {}/.gitkeep \;

echo "Added .gitkeep files to all empty directories."
