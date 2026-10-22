#!/usr/bin/env python3
"""Run one frequency headlessly. Thin wrapper around bbsim.run_frequency."""
import sys

from bbsim.run_frequency import main

if __name__ == "__main__":
    sys.exit(main())
