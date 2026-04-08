#!/usr/bin/env bash
set -euo pipefail

# Run the strain smoke test from the build directory.
# This script expects to be executed from the build directory.

filter_root_warnings() {
  grep -v -E 'std_darwin\.modulemap|__type_traits/add_lvalue_reference\.h|cling::IncrementalParser::CheckABICompatibility|Possible C\+\+ standard library mismatch|IncrementalExecutor::executeFunction|missing the definition of cling::runtime::gCling'
}

echo "Running strain smoke test..."
./strain 2> >(filter_root_warnings >&2)
echo "Strain smoke test passed."
