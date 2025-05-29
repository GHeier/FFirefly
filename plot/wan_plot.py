import numpy as np
import matplotlib.pyplot as plt
import tbmodels

def wan_plot(filename, mu):

    model = tbmodels.Model.from_wannier_files(hr_file=filename)

#print(model.hamilton(k=[0., 0., 0.]))
    k_path = np.concatenate([
        np.linspace([0, 0, 0], [0.5, 0, 0], 50),   # Γ to X
        np.linspace([0.5, 0, 0], [0.5, 0.5, 0], 50), # X to M
        np.linspace([0.5, 0.5, 0], [0, 0, 0], 50)  # M to Γ
    ])

    bands = model.eigenval(k_path)
    bands = np.array(bands)

    num_kpoints = k_path.shape[0]
    num_bands = bands.shape[1]

    for i in range(num_bands):
        plt.plot(range(num_kpoints), bands[:, i], color='black')

    plt.ylabel("Energy (eV)")
    plt.xlabel("k-point index")
    plt.title("Band Structure")
    plt.grid(True)
    plt.show()
