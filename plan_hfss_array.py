#!/usr/bin/env python3
"""Freeze a job-array run plan, with no AEDT. Thin wrapper around bbsim.array_plan."""
import sys

from bbsim.array_plan import main

if __name__ == "__main__":
    sys.exit(main())
