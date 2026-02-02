# RhoNNO MLX Integration

This package provides MLX-based implementations of neural networks compatible with the RhoNNO C++ library.

## Features

- MLX-based neural network classes (TMLP, TXMLP)
- Conversion tools for .net files to MLX models
- Compatible API with RhoNNO C++ classes

## Installation

```bash
pip install -e .
```

## Requirements

- Python 3.8+
- MLX 0.5.0+
- NumPy 1.21.0+

## Usage

```python
from rhonno_mlx import TMLP, NetConverter

# Convert a .net file to MLX model
converter = NetConverter()
mlp = converter.load_net('mlp.net')

# Or create a new model
mlp = TMLP(input_nodes=5, hidden_nodes=10, output_nodes=1)

# Train and predict
# ... training code ...
predictions = mlp.predict(inputs)
```