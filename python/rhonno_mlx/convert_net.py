#!/usr/bin/env python3
"""
Command-line tool for converting RhoNNO .net files to MLX models.

Usage:
    python -m rhonno_mlx.convert_net input.net output.pkl
"""

import sys
import pickle
from rhonno_mlx import NetConverter


def main():
    if len(sys.argv) != 3:
        print("Usage: python -m rhonno_mlx.convert_net <input.net> <output.pkl>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        converter = NetConverter()
        model = converter.load_net(input_file)

        # Save as pickle
        with open(output_file, 'wb') as f:
            pickle.dump(model, f)

        print(f"Successfully converted {input_file} to {output_file}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()