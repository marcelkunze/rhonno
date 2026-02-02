# RhoNNO - Neural Network Objects

RhoNNO is a C++ library for neural network implementations, part of the ROOT framework ecosystem.

## Features

- Multi-layer perceptrons (MLP)
- Self-organizing maps and neural gas networks
- Supervised and unsupervised learning algorithms
- Integration with ROOT data analysis framework

## MLX Integration

This project now includes MLX-based Python implementations for improved performance on Apple Silicon and other MLX-supported hardware.

### Python Package

The `python/` directory contains a Python package `rhonno_mlx` that provides:

- `TMLP`: MLX-based Multi-Layer Perceptron
- `TXMLP`: MLX-based Extended Multi-Layer Perceptron
- `NetConverter`: Tools to convert .net files to MLX models

### Installation

```bash
cd python
pip install -e .
```

### Usage

```python
from rhonno_mlx import TMLP, NetConverter

# Load existing .net file
mlp = TMLP(net_file='mlp.net')

# Or create new network
mlp = TMLP(input_nodes=5, hidden_nodes=10, output_nodes=1)

# Make predictions
import mlx.core as mx
inputs = mx.random.normal((batch_size, 5))
outputs = mlp.predict(inputs)
```

### Conversion Tools

Convert existing .net files to MLX format:

```bash
python -m rhonno_mlx.convert_net input.net output.pkl
```

### Testing

Run the test suite:

```bash
python python/test_mlx.py --all
```

Train on parity data:

```bash
python python/train_parity.py --epochs 100
```

## Building the C++ Library

### Prerequisites

- ROOT framework
- CMake 3.0+
- C++17 compatible compiler

### Build Instructions

```bash
mkdir build
cd build
cmake ..
make
```

## Documentation

For detailed API documentation, see the header files in the main directory.

## License

See LICENSE.txt for licensing information.

## Authors

- Johannes Steffens, Bochum University
- M.Kunze, Bochum University