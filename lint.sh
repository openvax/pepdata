#!/bin/bash
set -o errexit

ruff check pepdata/ tests/

echo 'Passes ruff check'
