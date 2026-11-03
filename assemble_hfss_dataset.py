#!/usr/bin/env python3
"""Select one accepted run per task of a job-array plan through the acceptance gate, with no AEDT, and print the report
with its guarded rerun lines; --report reports any plan (a probe's too), --check the same under the assembly's
refusals. --verify TREE [--report-path PATH] re-checks a written tree against BBRsim's loader rules with the standard
library only, reading no plan (exit 0 pass, 1 a failed check, 2 refused). Thin wrapper around bbsim.assemble."""
import sys

from bbsim.assemble import main

if __name__ == "__main__":
    sys.exit(main())
