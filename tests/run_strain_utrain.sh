#!/usr/bin/env bash
set -euo pipefail

# Run strain and utrain smoke tests from the build directory.
# This script expects to be executed from the build directory.

readonly ACTION=${1:-all}

filter_root_warnings() {
  grep -v -E 'std_darwin\.modulemap|__type_traits/add_lvalue_reference\.h|cling::IncrementalParser::CheckABICompatibility|Possible C\+\+ standard library mismatch|IncrementalExecutor::executeFunction|missing the definition of cling::runtime::gCling'
}

run_strain() {
    echo "Running strain smoke test..."
    ./strain 2> >(filter_root_warnings >&2)
}

run_utrain() {
    echo "Creating minimal ppe.dat for utrain..."
    cat > ppe.dat <<'EOF'
1.0 2.0
3.0 4.0
5.0 6.0
EOF
    echo "Running utrain smoke test..."
    ./utrain ppe.dat 2> >(filter_root_warnings >&2)
}

cleanup() {
    echo "Cleaning generated files..."
    rm -f ppe.dat parity.trn parity.tst mlp.net xmlp.net sgcs.net sgng.net gcs.net utrain.root
}

case "${ACTION}" in
    strain)
        run_strain
        ;;
    utrain)
        run_utrain
        ;;
    all)
        run_strain
        run_utrain
        ;;
    *)
        echo "Usage: $0 [strain|utrain|all]"
        exit 1
        ;;
esac

cleanup

echo "Smoke test completed successfully."
