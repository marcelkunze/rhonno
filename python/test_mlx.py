#!/usr/bin/env python3
"""
MLX-based test program for RhoNNO neural networks.

This script demonstrates loading and using MLX-based neural networks
that are compatible with RhoNNO .net files.
"""

import numpy as np
import mlx.core as mx
from rhonno_mlx import TMLP, NetConverter
import argparse


def test_mlp_creation():
    """Test creating a new MLP network."""
    print("Testing MLP creation...")

    # Create a simple MLP
    mlp = TMLP(input_nodes=5, hidden_nodes=10, output_nodes=1)

    # Test with random input
    test_input = mx.random.normal((1, 5))
    output = mlp.predict(test_input)

    print(f"Input shape: {test_input.shape}")
    print(f"Output shape: {output.shape}")
    print(f"Output value: {output}")
    print("✓ MLP creation test passed\n")


def test_net_loading():
    """Test loading a .net file."""
    print("Testing .net file loading...")

    try:
        # Load the existing mlp.net file
        converter = NetConverter()
        model = converter.load_net('../mlp.net')

        # Create MLP from loaded model
        mlp = TMLP(net_file='../mlp.net')

        # Test prediction
        test_input = mx.random.normal((1, 5))
        output = mlp.predict(test_input)

        print(f"Loaded network with {mlp.input_nodes} inputs, {mlp.hidden_nodes} hidden, {mlp.output_nodes} outputs")
        print(f"Test output shape: {output.shape}")
        print("✓ Net loading test passed\n")

    except FileNotFoundError:
        print("mlp.net not found, skipping loading test\n")
    except Exception as e:
        print(f"Error loading net file: {e}\n")


def test_parity_data():
    """Test with parity data if available."""
    print("Testing with parity data...")

    try:
        # Load parity training data
        train_data = np.loadtxt('../parity.trn')
        test_data = np.loadtxt('../parity.tst')

        print(f"Training data shape: {train_data.shape}")
        print(f"Test data shape: {test_data.shape}")

        # For parity data, assume last column is target
        X_train = train_data[:, :-1]
        y_train = train_data[:, -1:]
        X_test = test_data[:, :-1]
        y_test = test_data[:, -1:]

        # Create and test MLP
        input_nodes = X_train.shape[1]
        output_nodes = y_train.shape[1]

        mlp = TMLP(input_nodes=input_nodes, hidden_nodes=10, output_nodes=output_nodes)

        # Test prediction on a few samples
        test_inputs = mx.array(X_test[:5])
        predictions = mlp.predict(test_inputs)

        print(f"Sample predictions: {predictions}")
        print("✓ Parity data test passed\n")

    except FileNotFoundError:
        print("Parity data files not found, skipping test\n")
    except Exception as e:
        print(f"Error with parity data: {e}\n")


def test_conversion_tool():
    """Test the conversion tool."""
    print("Testing conversion tool...")

    try:
        import subprocess
        import os

        # Convert mlp.net to pickle
        result = subprocess.run([
            'python', '-m', 'rhonno_mlx.convert_net',
            '../mlp.net', 'test_model.pkl'
        ], capture_output=True, text=True)

        if result.returncode == 0:
            print("✓ Conversion tool test passed")

            # Load and verify
            with open('test_model.pkl', 'rb') as f:
                import pickle
                model = pickle.load(f)
            print(f"Converted model keys: {list(model.keys())}")

            # Clean up
            os.remove('test_model.pkl')
        else:
            print(f"Conversion failed: {result.stderr}")

        print()

    except Exception as e:
        print(f"Error testing conversion tool: {e}\n")


def main():
    parser = argparse.ArgumentParser(description="Test RhoNNO MLX integration")
    parser.add_argument("--all", action="store_true", help="Run all tests")
    parser.add_argument("--create", action="store_true", help="Test MLP creation")
    parser.add_argument("--load", action="store_true", help="Test net file loading")
    parser.add_argument("--parity", action="store_true", help="Test with parity data")
    parser.add_argument("--convert", action="store_true", help="Test conversion tool")

    args = parser.parse_args()

    if args.all or not any([args.create, args.load, args.parity, args.convert]):
        test_mlp_creation()
        test_net_loading()
        test_parity_data()
        test_conversion_tool()
    else:
        if args.create:
            test_mlp_creation()
        if args.load:
            test_net_loading()
        if args.parity:
            test_parity_data()
        if args.convert:
            test_conversion_tool()


if __name__ == "__main__":
    main()