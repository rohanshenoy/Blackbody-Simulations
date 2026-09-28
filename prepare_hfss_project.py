#!/usr/bin/env python3
"""Prepare a cleaned HFSS project. Thin wrapper around bbsim.prepare_project."""
import sys

from bbsim.prepare_project import main

if __name__ == "__main__":
    sys.exit(main())
