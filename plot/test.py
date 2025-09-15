import matplotlib.pyplot as plt

def plot_generator(data, plot_type='line', **kwargs):
    # Create the figure and axis
    fig, ax = plt.subplots()

    # Select plot type
    if plot_type == 'line':
        ax.plot(data['x'], data['y'], **kwargs)
    elif plot_type == 'scatter':
        ax.scatter(data['x'], data['y'], **kwargs)
    elif plot_type == 'bar':
        ax.bar(data['x'], data['y'], **kwargs)
    elif plot_type == 'hist':
        ax.hist(data['x'], **kwargs)
    else:
        raise ValueError(f"Unsupported plot type: {plot_type}")

    # Return the figure and axis for further modifications
    return fig, ax
# Example usage
import numpy as np

# Sample data
data = {
    'x': np.random.rand(50),
    'y': np.random.rand(50)
}

# Custom scatter parameters
scatter_kwargs = {
    's': 100,            # size of markers
    'c': 'blue',         # color
    'marker': '^'        # shape (triangle marker)
}

fig, ax = plot_generator(data, plot_type='scatter', **scatter_kwargs)

# Optional: add title, labels, etc.
ax.set_title("Customized Scatter Plot")
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")

plt.show()

