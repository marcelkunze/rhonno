#!/usr/bin/env bash
set -euo pipefail

# Run the utrain smoke test from the build directory.
# This script expects to be executed from the build directory.

filter_root_warnings() {
  grep -v -E 'std_darwin\.modulemap|__type_traits/add_lvalue_reference\.h|cling::IncrementalParser::CheckABICompatibility|Possible C\+\+ standard library mismatch|IncrementalExecutor::executeFunction|missing the definition of cling::runtime::gCling'
}

echo "Creating minimal ppe.dat for utrain..."
cat > ppe.dat <<'EOF'
1.0 2.0
3.0 4.0
5.0 6.0
EOF

echo "Running utrain smoke test..."
./utrain ppe.dat 2> >(filter_root_warnings >&2)

echo "Cleaning generated files..."
rm -f ppe.dat

echo "Utrain smoke test passed."
