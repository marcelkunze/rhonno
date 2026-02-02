#!/usr/bin/env python3
"""
MLX-based supervised training example for RhoNNO.

This script demonstrates training MLX-based neural networks on the parity problem,
similar to the C++ strain.cxx example.
"""

import numpy as np
import mlx.core as mx
from rhonno_mlx import TMLP
import argparse
import time


def generate_parity_data(n_samples=1000, n_bits=5):
    """
    Generate parity training data.

    Args:
        n_samples: Number of training samples
        n_bits: Number of input bits

    Returns:
        X: Input array of shape (n_samples, n_bits)
        y: Target array of shape (n_samples, 1)
    """
    X = np.random.randint(0, 2, (n_samples, n_bits)).astype(np.float32)

    # Calculate parity (1 for even, -1 for odd)
    parity = np.sum(X, axis=1) % 2
    y = np.where(parity == 0, 1.0, -1.0).reshape(-1, 1).astype(np.float32)

    return X, y


def train_mlp_simple(mlp, X_train, y_train, epochs=100, learning_rate=0.1):
    """
    Simple training loop (placeholder - full MLX training would use mlx.nn).

    Note: This is a simplified version. For production use, implement
    proper backpropagation with MLX's automatic differentiation.
    """
    print(f"Training MLP for {epochs} epochs...")

    X_mx = mx.array(X_train)
    y_mx = mx.array(y_train)

    for epoch in range(epochs):
        # Forward pass
        predictions = mlp.forward(X_mx)

        # Simple loss (MSE)
        loss = mx.mean((predictions - y_mx) ** 2)

        if epoch % 10 == 0:
            loss_val = float(loss)
            print(f"Epoch {epoch}, Loss: {loss_val:.6f}")

        # Placeholder for gradient descent
        # In a real implementation, use MLX's grad and optimizers
        # For now, just demonstrate the forward pass

    print("Training completed (simplified version)")


def test_network(mlp, X_test, y_test):
    """Test the trained network."""
    X_mx = mx.array(X_test)
    predictions = mlp.predict(X_mx)

    # Convert predictions to binary (-1 or 1)
    pred_binary = mx.where(predictions > 0, 1.0, -1.0)

    # Calculate accuracy
    correct = mx.sum(pred_binary == mx.array(y_test))
    accuracy = float(correct) / len(y_test)

    print(f"Test Accuracy: {accuracy:.4f}")
    return accuracy


def main():
    parser = argparse.ArgumentParser(description="MLX-based supervised training example")
    parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs")
    parser.add_argument("--samples", type=int, default=1000, help="Number of training samples")
    parser.add_argument("--bits", type=int, default=5, help="Number of input bits")
    parser.add_argument("--save", type=str, help="Save trained network to file")

    args = parser.parse_args()

    print("RhoNNO MLX Supervised Training Example")
    print("=" * 40)

    # Generate training data
    print(f"Generating {args.samples} parity training samples with {args.bits} bits...")
    X_train, y_train = generate_parity_data(args.samples, args.bits)
    X_test, y_test = generate_parity_data(200, args.bits)

    print(f"Training data shape: {X_train.shape}")
    print(f"Test data shape: {X_test.shape}")

    # Create MLP
    mlp = TMLP(input_nodes=args.bits, hidden_nodes=10, output_nodes=1)
    print(f"Created MLP: {args.bits} -> 10 -> 1")

    # Train network
    start_time = time.time()
    train_mlp_simple(mlp, X_train, y_train, epochs=args.epochs)
    training_time = time.time() - start_time
    print(".2f")

    # Test network
    accuracy = test_network(mlp, X_test, y_test)

    # Save network if requested
    if args.save:
        mlp.save_net(args.save)
        print(f"Network saved to {args.save}")

    print("\nTraining completed successfully!")
    print(f"Final test accuracy: {accuracy:.4f}")


if __name__ == "__main__":
    main()