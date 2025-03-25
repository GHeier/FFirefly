import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot_points_and_fit(file, func):
    # Read the data into a Pandas DataFrame
    data = pd.read_csv(file, delim_whitespace=True)
    
    # Access columns by their names
    x = data['Temperature']
    y = data['Max_phi']
    
    # Perform the curve fit
    params, covariance = curve_fit(func, x, y)
    
    # Generate fitted y values
    y_fit = func(x, *params)
    
    # Plot the data and the fit
    plt.scatter(x, y, marker='+', label='Data')  # Use marker='+' for the style
    #plt.plot(x, y_fit, label="Fit", color="red")  # Fitted curve
    plt.xlabel("Temperature")
    plt.ylabel("Max_phi")
    plt.legend()
    plt.title("Curve Fitting")
    plt.show()

# Example usage
if __name__ == "__main__":
    # Define your model function
    def func(x, a, b, c):
        return a * x**2 + b * x + c  # Replace with your actual model
    
    # Call the function with your data file and model
    plot_points_and_fit('Phi_v_T.dat', func)
