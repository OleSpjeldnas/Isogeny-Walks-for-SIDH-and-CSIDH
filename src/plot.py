import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Replace 'data.txt' with the path to your actual data file
file_path = '/Users/ole/Downloads/Isogeny-Walk-zk-SNARKs-for-SIDH-and-CSIDH/results_new_2.txt'

# Read the data into a pandas DataFrame
a = 512
df = pd.read_csv(file_path, header=None, names=['i', 'prover_time_seconds', 'verifier_time_ms', 'proof_size_kB'], sep=',')
df['input_size'] = a * (2 ** df['i'])
# Plotting


# Prover Time vs. i
plt.figure(figsize=(10, 6))
plt.plot(df['input_size'], df['prover_time_seconds'], marker='o', linestyle='-', color='blue')
plt.title('Prover Time')
plt.xlabel('Witness Size')
plt.ylabel('Prover Time (s)')
plt.xticks(np.power(2, np.arange(9, 14)), [f'2^{i}' for i in range(9, 14)])  # Set x-ticks to 2^i for i from 9 to 14
plt.grid(False)
plt.show()

# Verifier Time vs. i
plt.figure(figsize=(10, 6))
plt.plot(df['input_size'], df['verifier_time_ms'], marker='o', linestyle='-', color='green')
plt.title('Verifier Time')
plt.xlabel('Witness Size')
plt.ylabel('Verifier Time (ms)')
plt.xticks(np.power(2, np.arange(9, 14)), [f'2^{i}' for i in range(9, 14)])  # Set x-ticks to 2^i for i from 9 to 14
plt.grid(False)
plt.show()

# Proof Size vs. i
plt.figure(figsize=(10, 6))
plt.plot(df['input_size'], df['proof_size_kB'], marker='o', linestyle='-', color='red')
plt.title('Proof Size')
plt.xlabel('Witness Size')
plt.ylabel('Proof Size (kB)')
plt.xticks(np.power(2, np.arange(9, 14)), [f'2^{i}' for i in range(9, 14)])  # Set x-ticks to 2^i for i from 9 to 14
plt.grid(False)
plt.show()