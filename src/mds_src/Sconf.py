def read_numbers(filename):
    with open(filename, 'r') as file:
        numbers = []
        for line in file:
            for item in line.split(' '):
                numbers.append(float(item.strip()))
    return numbers

def calculate_median(numbers,n_numbers):
    numbers.sort()
    n = n_numbers
    mid = n // 2
    if n % 2 == 0:
        return (numbers[mid - 1] + numbers[mid]) / 2.0
    else:
        return numbers[mid]

# Example usage
n_numbers = 4
filename = 'Sconf.txt'  # Ensure this file exists and contains numbers
numbers = read_numbers(filename)
filename = 'Sconf_min.txt'  # Ensure this file exists and contains numbers
min_numbers = read_numbers(filename)
n_pairs = len(numbers) // n_numbers
import numpy as np
n_elms = np.sqrt(n_pairs)
i_elm = 1
j_elm = 1
for i in range(n_pairs):
    median = calculate_median(numbers[i*n_numbers:(i+1)*n_numbers],n_numbers)
    if median < 3:
        phase = 'solid'
        print(f"The pair between element {i_elm}-{j_elm} appears to be {phase}. I suggest that you take the minimum: {min_numbers[i]}")
    else:
        phase = 'liquid'
        print(f"The pair between element {i_elm}-{j_elm} appears to be {phase}. I suggest that you take the median:  {median}")
    j_elm = j_elm + 1
    if j_elm > n_elms:
        j_elm = 1
        i_elm = i_elm + 1
