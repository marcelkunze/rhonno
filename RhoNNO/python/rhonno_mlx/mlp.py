"""
MLX-based implementation of TMLP (Multi-Layer Perceptron).

Compatible with RhoNNO TMLP class but using MLX for computation.
"""

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    print("Warning: NumPy not available")

try:
    import mlx.core as mx
    MLX_AVAILABLE = True
except ImportError:
    MLX_AVAILABLE = False
    print("Warning: MLX not available, using NumPy only")


class TMLP:
    """MLX-based Multi-Layer Perceptron with one hidden layer."""

    def __init__(self, input_nodes=None, hidden_nodes=None, output_nodes=None,
                 net_file=None, input_range=1.0, learn_step=0.1):
        """
        Initialize TMLP.

        Args:
            input_nodes: Number of input nodes
            hidden_nodes: Number of hidden nodes
            output_nodes: Number of output nodes
            net_file: Path to .net file to load
            input_range: Input scaling range
            learn_step: Learning step size
        """
        if net_file:
            from .converter import NetConverter
            converter = NetConverter()
            self.model = converter.load_net(net_file)
            # Extract dimensions from loaded model
            self.input_nodes = self.model['hidden_weights'].shape[0]
            self.hidden_nodes = self.model['hidden_weights'].shape[1]
            self.output_nodes = self.model['output_weights'].shape[1]
        else:
            if None in [input_nodes, hidden_nodes, output_nodes]:
                raise ValueError("Must specify input_nodes, hidden_nodes, output_nodes or net_file")

            self.input_nodes = input_nodes
            self.hidden_nodes = hidden_nodes
            self.output_nodes = output_nodes
            self.input_range = input_range
            self.learn_step = learn_step

            # Initialize weights and biases
            self._initialize_weights()

    def _initialize_weights(self):
        """Initialize weights and biases randomly."""
        # Xavier/Glorot initialization
        hidden_limit = np.sqrt(6 / (self.input_nodes + self.hidden_nodes))
        output_limit = np.sqrt(6 / (self.hidden_nodes + self.output_nodes))

        if MLX_AVAILABLE:
            self.model = {
                'hidden_weights': mx.random.uniform(-hidden_limit, hidden_limit,
                                                  (self.input_nodes, self.hidden_nodes)),
                'hidden_biases': mx.zeros((self.hidden_nodes,)),
                'output_weights': mx.random.uniform(-output_limit, output_limit,
                                                  (self.hidden_nodes, self.output_nodes)),
                'output_biases': mx.zeros((self.output_nodes,)),
                'input_scale': self.input_range,
                'transfer_func': 'fermi'
            }
        else:
            self.model = {
                'hidden_weights': np.random.uniform(-hidden_limit, hidden_limit,
                                                  (self.input_nodes, self.hidden_nodes)),
                'hidden_biases': np.zeros((self.hidden_nodes,)),
                'output_weights': np.random.uniform(-output_limit, output_limit,
                                                  (self.hidden_nodes, self.output_nodes)),
                'output_biases': np.zeros((self.output_nodes,)),
                'input_scale': self.input_range,
                'transfer_func': 'fermi'
            }

    def _sigmoid(self, x):
        """Sigmoid activation function."""
        if MLX_AVAILABLE and isinstance(x, mx.array):
            return 1 / (1 + mx.exp(-x))
        else:
            return 1 / (1 + np.exp(-x))

    def _fermi(self, x):
        """Fermi function (sigmoid)."""
        return self._sigmoid(x)

    def forward(self, inputs):
        """
        Forward pass through the network.

        Args:
            inputs: Input array of shape (batch_size, input_nodes)

        Returns:
            Output array of shape (batch_size, output_nodes)
        """
        # Scale inputs
        x = inputs * self.model['input_scale']

        if MLX_AVAILABLE and isinstance(x, mx.array):
            # Hidden layer
            hidden = mx.matmul(x, self.model['hidden_weights']) + self.model['hidden_biases']
            hidden_activated = self._fermi(hidden)

            # Output layer
            output = mx.matmul(hidden_activated, self.model['output_weights']) + self.model['output_biases']
            output_activated = self._sigmoid(output)  # Output also uses sigmoid
        else:
            # NumPy fallback
            # Hidden layer
            hidden = np.dot(x, self.model['hidden_weights']) + self.model['hidden_biases']
            hidden_activated = self._fermi(hidden)

            # Output layer
            output = np.dot(hidden_activated, self.model['output_weights']) + self.model['output_biases']
            output_activated = self._sigmoid(output)  # Output also uses sigmoid

        return output_activated

    def predict(self, inputs):
        """
        Make predictions.

        Args:
            inputs: Input array

        Returns:
            Predictions
        """
        if MLX_AVAILABLE:
            if isinstance(inputs, np.ndarray):
                inputs = mx.array(inputs)
            elif not isinstance(inputs, mx.array):
                inputs = mx.array(inputs)
        else:
            if not isinstance(inputs, np.ndarray):
                inputs = np.array(inputs)

        return self.forward(inputs)

    def save_net(self, filepath):
        """
        Save network to .net file format.

        Args:
            filepath: Path to save the network
        """
        # Convert arrays to numpy for saving
        if MLX_AVAILABLE and isinstance(self.model['hidden_weights'], mx.array):
            hidden_weights = np.array(self.model['hidden_weights'])
            hidden_biases = np.array(self.model['hidden_biases'])
            output_weights = np.array(self.model['output_weights'])
            output_biases = np.array(self.model['output_biases'])
        else:
            hidden_weights = self.model['hidden_weights']
            hidden_biases = self.model['hidden_biases']
            output_weights = self.model['output_weights']
            output_biases = self.model['output_biases']

        with open(filepath, 'w') as f:
            f.write("C++  NEURAL NETWORK OBJECTS   VERSION 2.0ROOT\n")
            f.write("Filetype text\n\n")
            f.write("network id  XMLP\n")
            f.write(f"innodes     {self.input_nodes}\n")
            f.write(f"outnodes    {self.output_nodes}\n")
            f.write("layers    2\n")
            f.write(f"in_scale  {self.model['input_scale']:.6e}\n\n")

            # Hidden layer perceptron
            f.write("Perceptron ID 0\n")
            f.write(f"innodes     {self.input_nodes}\n")
            f.write(f"outnodes    {self.hidden_nodes}\n")
            f.write(f"learn_step  {self.learn_step:.6e}\n")
            f.write("transfer_id 1\n\n")  # Fermi

            for i in range(self.hidden_nodes):
                f.write(f"unit number      {i}\n")
                f.write(f"threshold        {hidden_biases[i]:.6e}\n")
                f.write("weights\n")
                for j in range(self.input_nodes):
                    f.write(f"{hidden_weights[j, i]:.6e}\n")
                f.write("\n")

            # Output layer perceptron
            f.write("Perceptron ID 1\n")
            f.write(f"innodes     {self.hidden_nodes}\n")
            f.write(f"outnodes    {self.output_nodes}\n")
            f.write(f"learn_step  {self.learn_step:.6e}\n")
            f.write("transfer_id 1\n\n")  # Fermi

            for i in range(self.output_nodes):
                f.write(f"unit number      {i}\n")
                f.write(f"threshold        {output_biases[i]:.6e}\n")
                f.write("weights\n")
                for j in range(self.hidden_nodes):
                    f.write(f"{output_weights[j, i]:.6e}\n")
                f.write("\n")