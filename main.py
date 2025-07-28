from rings import find_rings
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt


def main():
    # Read test structure
    ats = read('test.extxyz')

    # Find the rings
    rings = find_rings(ats, radii_factor=1.3, bonds=[('P', 'O'), ('Ca', 'O')], repeat=(1, 1, 1))
    print(f'Found {len(rings)} rings')

    # Create a list of structures containing only the rings
    rs = []
    for r in rings:
        ss = ats[r]
        rs.append(ss)

    # Save it to have a look at it
    write('rings.extxyz', rs)


    # Calculate ring lengths
    ring_lengths = [len(r) // 2 for r in rings]

    # Create the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(ring_lengths, bins=range(1, max(ring_lengths) + 2), align='left', 
             color='skyblue', edgecolor='black', rwidth=0.8)

    # Customize the plot
    plt.title('Distribution of Ring Sizes in Phosphate Network', fontsize=12, pad=15)
    plt.xlabel('Number of Atoms in Ring', fontsize=10)
    plt.ylabel('Frequency', fontsize=10)

    # Set x-ticks to show integer values
    plt.xticks(range(min(ring_lengths), max(ring_lengths) + 1))

    # Add grid for better readability
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()

