#!/usr/bin/env bash
set -euo pipefail

# Run the strain smoke test from the build directory.
# This script expects to be executed from the build directory.

echo "Running strain smoke test..."
./strain
 echo "Strain smoke test passed."
