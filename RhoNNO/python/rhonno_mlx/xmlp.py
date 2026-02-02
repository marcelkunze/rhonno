"""
MLX-based implementation of TXMLP (Extended Multi-Layer Perceptron).

Compatible with RhoNNO TXMLP class but using MLX for computation.
"""

import mlx.core as mx
import numpy as np
from typing import List, Optional


class TXMLP:
    """MLX-based Extended Multi-Layer Perceptron."""

    def __init__(self, layers=None, input_range=1.0, net_file=None,
                 input_nodes=None, hidden_nodes=None, output_nodes=None,
                 learn_steps=None, transfer_funcs=None):
        """
        Initialize TXMLP.

        Args:
            layers: Number of layers (including input and output)
            input_range: Input scaling range
            net_file: Path to .net file to load
            input_nodes: Number of input nodes (if not loading from file)
            hidden_nodes: List of hidden layer sizes (if not loading from file)
            output_nodes: Number of output nodes (if not loading from file)
            learn_steps: Learning steps for each layer
            transfer_funcs: Transfer functions for each layer
        """
        if net_file:
            from .converter import NetConverter
            converter = NetConverter()
            self.model = converter.load_net(net_file)
            # For now, assume loaded model has the structure
            self.layers = 3  # input + hidden + output
            self.input_range = self.model['input_scale']
        else:
            if layers is None:
                # Assume 2 hidden layers for MLP-like network
                layers = 3
            if hidden_nodes is None:
                hidden_nodes = [10]  # Default single hidden layer
            if isinstance(hidden_nodes, int):
                hidden_nodes = [hidden_nodes]

            self.layers = layers
            self.input_range = input_range
            self.input_nodes = input_nodes
            self.hidden_nodes = hidden_nodes
            self.output_nodes = output_nodes

            if learn_steps is None:
                self.learn_steps = [0.1] * (layers - 1)
            else:
                self.learn_steps = learn_steps

            if transfer_funcs is None:
                self.transfer_funcs = ['fermi'] * (layers - 1)
            else:
                self.transfer_funcs = transfer_funcs

            self._initialize_weights()

    def _initialize_weights(self):
        """Initialize weights and biases for multi-layer network."""
        layer_sizes = [self.input_nodes] + self.hidden_nodes + [self.output_nodes]

        self.weights = []
        self.biases = []

        for i in range(len(layer_sizes) - 1):
            in_size = layer_sizes[i]
            out_size = layer_sizes[i + 1]

            # Xavier initialization
            limit = np.sqrt(6 / (in_size + out_size))
            w = mx.random.uniform(-limit, limit, (in_size, out_size))
            b = mx.zeros((out_size,))

            self.weights.append(w)
            self.biases.append(b)

        self.model = {
            'weights': self.weights,
            'biases': self.biases,
            'input_scale': self.input_range,
            'transfer_funcs': self.transfer_funcs
        }

    def _apply_activation(self, x, func_name):
        """Apply activation function."""
        if func_name == 'fermi' or func_name == 'sigmoid':
            return 1 / (1 + mx.exp(-x))
        elif func_name == 'linear':
            return x
        else:
            # Default to sigmoid
            return 1 / (1 + mx.exp(-x))

    def forward(self, inputs):
        """
        Forward pass through the network.

        Args:
            inputs: Input array of shape (batch_size, input_nodes)

        Returns:
            Output array of shape (batch_size, output_nodes)
        """
        # Scale inputs
        x = inputs * self.input_range

        # Forward through layers
        for i, (w, b) in enumerate(zip(self.weights, self.biases)):
            x = mx.matmul(x, w) + b
            x = self._apply_activation(x, self.transfer_funcs[i])

        return x

    def predict(self, inputs):
        """
        Make predictions.

        Args:
            inputs: Input array

        Returns:
            Predictions
        """
        if isinstance(inputs, np.ndarray):
            inputs = mx.array(inputs)
        elif not isinstance(inputs, mx.array):
            inputs = mx.array(inputs)

        return self.forward(inputs)

    def train_step(self, inputs, targets, learning_rate=0.01):
        """
        Single training step (simplified - no full backprop implemented yet).

        Args:
            inputs: Input batch
            targets: Target outputs
            learning_rate: Learning rate

        Returns:
            Loss value
        """
        # This is a placeholder - full training implementation would require
        # automatic differentiation with MLX
        predictions = self.forward(inputs)
        loss = mx.mean((predictions - targets) ** 2)
        return loss

    def save_net(self, filepath):
        """
        Save network to .net file format.

        Note: This is a simplified version for 2-layer networks.
        """
        # For now, save as XMLP format
        with open(filepath, 'w') as f:
            f.write("C++  NEURAL NETWORK OBJECTS   VERSION 2.0ROOT\n")
            f.write("Filetype text\n\n")
            f.write("network id  XMLP\n")
            f.write(f"innodes     {self.input_nodes}\n")
            f.write(f"outnodes    {self.output_nodes}\n")
            f.write(f"layers    {self.layers}\n")
            f.write(f"in_scale  {self.input_range:.6e}\n\n")

            # Save each layer as a perceptron
            for i, (w, b) in enumerate(zip(self.weights, self.biases)):
                f.write(f"Perceptron ID {i}\n")
                in_nodes = w.shape[0]
                out_nodes = w.shape[1]
                f.write(f"innodes     {in_nodes}\n")
                f.write(f"outnodes    {out_nodes}\n")
                f.write(f"learn_step  {self.learn_steps[i]:.6e}\n")
                transfer_id = 1 if self.transfer_funcs[i] == 'fermi' else 2
                f.write(f"transfer_id {transfer_id}\n\n")

                w_np = np.array(w)
                b_np = np.array(b)

                for j in range(out_nodes):
                    f.write(f"unit number      {j}\n")
                    f.write(f"threshold        {b_np[j]:.6e}\n")
                    f.write("weights\n")
                    for k in range(in_nodes):
                        f.write(f"{w_np[k, j]:.6e}\n")
                    f.write("\n")