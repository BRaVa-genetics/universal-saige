#!/bin/bash
# Renamed to 00a_sigma_check.sh (step 0a): it runs straight after step 0. This forwards to it.
echo "05_sigma_check.sh is now 00a_sigma_check.sh (step 0a); forwarding." >&2
exec bash "$(dirname "$0")/00a_sigma_check.sh" "$@"
