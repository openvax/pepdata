#!/bin/bash
set -o errexit

pytest --cov=pepdata/ --cov-report=term-missing tests "$@"


