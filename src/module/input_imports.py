INPUTS = {
    "hamiltonian/band_structure.hpp": [
        ("epsilon", ["float", "int", "Vec"]),
    ],
    "objects/CMField/bands.hpp": [
        ("Bands", ["ptr"]),
        (
            "Bands_operator",
            ["float", "ptr", "int", "Vec"],
        ),
    ],
    "objects/CMField/vertex.hpp": [
        ("Vertex", ["ptr"]),
        (
            "Vertex_operator",
            ["complex<float>", "ptr", "Vec", "float=0.0", "string='up'", "string='up'"],
        ),
    ],
    # "objects/surfaces.hpp": [
    #    ("Surface", ["ptr", "func", "float"]),
    #    ("Surface_var_faces", ["vector<Vec>"]),
    # ],
    "objects/CMField/fields.hpp": [
        ("Field_R", ["ptr"]),
        ("Field_R", ["ptr", "string"]),
        ("Field_R_operator", ["float", "ptr", "float"]),
        ("Field_R_operator", ["float", "ptr", "int", "float"]),
        ("Field_R_operator", ["float", "ptr", "Vec", "float=0.0"]),
        ("Field_R_operator", ["float", "ptr", "int", "Vec", "float=0.0"]),
        # ---------------------------------------------------------------------
        ("Field_C", ["ptr"]),
        ("Field_C", ["ptr", "string"]),
        ("Field_C_operator", ["complex<float>", "ptr", "float"]),
        ("Field_C_operator", ["complex<float>", "ptr", "int", "float"]),
        ("Field_C_operator", ["complex<float>", "ptr", "Vec", "float=0.0"]),
        ("Field_C_operator", ["complex<float>", "ptr", "int", "Vec", "float=0.0"]),
    ],
}


from write import write_export, write_import

write_export.write_export_funcs(INPUTS)
print("Successfully wrote export files")
write_import.write_import(INPUTS)
print("Successfully wrote import files")
