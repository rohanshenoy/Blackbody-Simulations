#!/usr/bin/env python3
"""The fibre bore's ports-only mode run (OR-053). Thin wrapper around bbsim.run_modes."""
import sys

from bbsim.run_modes import main

if __name__ == "__main__":
    sys.exit(main())
