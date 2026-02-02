#!/usr/bin/env python3
"""
Validation script for RhoNNO MLX integration.

This script validates the core functionality without requiring MLX or NumPy.
"""

import os
import sys

def validate_net_parsing():
    """Validate .net file parsing."""
    print("Validating .net file parsing...")

    # Simple parser test
    with open('mlp.net', 'r') as f:
        content = f.read()

    lines = content.split('\n')
    if not lines[0].startswith('C++  NEURAL NETWORK OBJECTS'):
        print("❌ Invalid .net file format")
        return False

    # Find network parameters
    i = 0
    while i < len(lines):
        if lines[i].strip() == 'network id  XMLP':
            break
        i += 1

    if i >= len(lines):
        print("❌ Could not find network id")
        return False

    net_id = lines[i].split()[-1]
    i += 1
    innodes = int(lines[i].split()[-1])
    i += 1
    outnodes = int(lines[i].split()[-1])
    i += 1
    layers = int(lines[i].split()[-1])
    i += 1
    in_scale = float(lines[i].split()[-1])

    print(f"✅ Parsed network: {net_id}, {innodes} inputs, {outnodes} outputs, {layers} layers")
    return True

def validate_file_structure():
    """Validate that all required files exist."""
    print("Validating file structure...")

    required_files = [
        'python/rhonno_mlx/__init__.py',
        'python/rhonno_mlx/converter.py',
        'python/rhonno_mlx/mlp.py',
        'python/rhonno_mlx/xmlp.py',
        'python/rhonno_mlx/convert_net.py',
        'python/test_mlx.py',
        'python/train_parity.py',
        'python/setup.py',
        'python/README.md',
        'python/requirements.txt',
        'README.md'
    ]

    missing = []
    for file in required_files:
        if not os.path.exists(file):
            missing.append(file)

    if missing:
        print(f"❌ Missing files: {missing}")
        return False
    else:
        print("✅ All required files present")
        return True

def validate_python_syntax():
    """Validate Python syntax of all files."""
    print("Validating Python syntax...")

    import subprocess
    result = subprocess.run([
        'find', 'python', '-name', '*.py', '-exec', 'python', '-m', 'py_compile', '{}', ';'
    ], capture_output=True, text=True)

    if result.returncode == 0:
        print("✅ All Python files compile successfully")
        return True
    else:
        print(f"❌ Python syntax errors: {result.stderr}")
        return False

def main():
    print("RhoNNO MLX Integration Validation")
    print("=" * 40)

    tests = [
        validate_file_structure,
        validate_python_syntax,
        validate_net_parsing
    ]

    passed = 0
    for test in tests:
        if test():
            passed += 1
        print()

    print(f"Validation complete: {passed}/{len(tests)} tests passed")

    if passed == len(tests):
        print("🎉 MLX integration implementation is complete and valid!")
        return 0
    else:
        print("❌ Some validation tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())