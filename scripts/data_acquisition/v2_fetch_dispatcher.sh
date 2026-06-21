#!/bin/bash
# v2_fetch_dispatcher.sh - disabled
#
# This legacy dispatcher submitted generic `fetch_srr` jobs in parallel with the
# v2 supervisor, which could lead to the same SRR being fetched concurrently.
# Keep this script as a hard stop so any old cron/screen/manual invocation fails
# loudly instead of launching more jobs.

echo "ERROR: v2_fetch_dispatcher.sh is disabled." >&2
echo "Use scripts/v2_supervisor.sh only." >&2
exit 1
