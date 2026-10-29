#!/usr/bin/env python3
"""Freeze a job-array run plan, with no AEDT, and print its guarded submission block; --submit submits its jobs
through the licence guard, and --log-submission RUN KIND ID records a job submitted by hand (no admission). Thin
wrapper around bbsim.array_plan."""
import sys

from bbsim.array_plan import main

if __name__ == "__main__":
    sys.exit(main())
