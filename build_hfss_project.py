#!/usr/bin/env python3
"""Build an HFSS project from a geometry spec. Thin wrapper around bbsim.build_project."""
import sys

from bbsim.build_project import main

if __name__ == "__main__":
    sys.exit(main())
