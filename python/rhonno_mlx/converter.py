"""
Net file converter for RhoNNO MLX integration.

This module provides tools to convert RhoNNO .net files to MLX models.
"""

import re
import numpy as np

try:
    import mlx.core as mx
    MLX_AVAILABLE = True
except ImportError:
    MLX_AVAILABLE = False
    print("Warning: MLX not available, using NumPy arrays only")


class NetConverter:
    """Converter for RhoNNO .net files to MLX models."""

    def __init__(self):
        self.transfer_functions = {
            0: "user",
            1: "fermi",  # sigmoid
            2: "linear",
            3: "linear_bend",
            4: "sigmoid"
        }

    def parse_net_file(self, filepath):
        """Parse a .net file and return network parameters."""
        with open(filepath, 'r') as f:
            content = f.read()

        # Parse header
        lines = content.split('\n')
        if not lines[0].startswith('C++  NEURAL NETWORK OBJECTS'):
            raise ValueError("Invalid .net file format")

        # Skip to network parameters
        i = 0
        while i < len(lines):
            if lines[i].strip() == 'network id  XMLP':
                break
            i += 1

        if i >= len(lines):
            raise ValueError("Could not find network id")

        # Parse network parameters
        params = {}
        params['net_id'] = lines[i].split()[-1]

        i += 1
        params['innodes'] = int(lines[i].split()[-1])
        i += 1
        params['outnodes'] = int(lines[i].split()[-1])
        i += 1
        params['layers'] = int(lines[i].split()[-1])
        i += 1
        params['in_scale'] = float(lines[i].split()[-1])

        # Parse perceptrons
        perceptrons = []
        while i < len(lines):
            if lines[i].strip().startswith('Perceptron ID'):
                perc_id = int(lines[i].split()[-1])
                i += 1
                innodes = int(lines[i].split()[-1])
                i += 1
                outnodes = int(lines[i].split()[-1])
                i += 1
                learn_step = float(lines[i].split()[-1])
                i += 1
                transfer_id = int(lines[i].split()[-1])

                # Parse units
                units = []
                for unit_idx in range(outnodes):
                    i += 1  # unit number
                    i += 1  # threshold
                    threshold = float(lines[i].split()[-1])
                    i += 1  # weights

                    weights = []
                    for _ in range(innodes):
                        i += 1
                        weights.append(float(lines[i].split()[-1]))

                    units.append({
                        'threshold': threshold,
                        'weights': weights
                    })

                perceptrons.append({
                    'id': perc_id,
                    'innodes': innodes,
                    'outnodes': outnodes,
                    'learn_step': learn_step,
                    'transfer_id': transfer_id,
                    'units': units
                })

            i += 1
            if i >= len(lines):
                break

        params['perceptrons'] = perceptrons
        return params

    def create_mlx_model(self, params):
        """Create an MLX model from parsed parameters."""
        if params['net_id'] == 'XMLP':
            return self._create_xmlp_model(params)
        else:
            raise ValueError(f"Unsupported network type: {params['net_id']}")

    def _create_xmlp_model(self, params):
        """Create TXMLP model from parameters."""
        # For simplicity, assume 2-layer network (input + hidden + output)
        if len(params['perceptrons']) != 2:
            raise ValueError("Only 2-layer networks supported for now")

        hidden_perc = params['perceptrons'][0]
        output_perc = params['perceptrons'][1]

        # Extract weights and biases
        hidden_weights = []
        hidden_biases = []
        for unit in hidden_perc['units']:
            hidden_weights.append(unit['weights'])
            hidden_biases.append(unit['threshold'])

        output_weights = []
        output_biases = []
        for unit in output_perc['units']:
            output_weights.append(unit['weights'])
            output_biases.append(unit['threshold'])

        # Convert to numpy arrays
        hidden_weights = np.array(hidden_weights).T  # (input_size, hidden_size)
        hidden_biases = np.array(hidden_biases)      # (hidden_size,)
        output_weights = np.array(output_weights).T  # (hidden_size, output_size)
        output_biases = np.array(output_biases)      # (output_size,)

        if MLX_AVAILABLE:
            # Create MLX arrays
            hidden_weights_mx = mx.array(hidden_weights)
            hidden_biases_mx = mx.array(hidden_biases)
            output_weights_mx = mx.array(output_weights)
            output_biases_mx = mx.array(output_biases)
        else:
            # Use NumPy arrays
            hidden_weights_mx = hidden_weights
            hidden_biases_mx = hidden_biases
            output_weights_mx = output_weights
            output_biases_mx = output_biases

        return {
            'hidden_weights': hidden_weights_mx,
            'hidden_biases': hidden_biases_mx,
            'output_weights': output_weights_mx,
            'output_biases': output_biases_mx,
            'input_scale': params['in_scale'],
            'transfer_func': self.transfer_functions.get(hidden_perc['transfer_id'], 'fermi')
        }

    def load_net(self, filepath):
        """Load a .net file and return an MLX model."""
        params = self.parse_net_file(filepath)
        return self.create_mlx_model(params)