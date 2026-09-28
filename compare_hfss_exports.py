#!/usr/bin/env python3
"""Compare HFSS export datasets. Thin wrapper around bbsim.compare_exports."""
import sys

from bbsim.compare_exports import main

if __name__ == "__main__":
    sys.exit(main())
