INPUTS = {
    "hamiltonian/band_structure.hpp": {
        "epsilon": ["float", "int", "Vec"],
        "Bands": ["ptr"],
    },
}


from write import write_export, write_import

write_export.write_export(INPUTS)
write_import.write_import(INPUTS)
