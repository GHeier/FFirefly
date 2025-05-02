INPUTS = {
    "hamiltonian/band_structure.hpp": [
        ("epsilon", ["float", "int", "Vec"]),
    ],
    "objects/CMField/bands.hpp": [
        ("Bands", ["ptr"]),
    ],
    "objects/CMField/vertex.hpp": [
        ("Vertex", ["ptr"]),
    ],
    "objects/CMField/fields.hpp": [
        ("Field_R", ["ptr"]),
        ("Field_R", ["ptr", "string"]),
        ("Field_R_operator", ["float", "ptr", "float"]),
        ("Field_R_operator", ["float", "ptr", "int", "float"]),
        ("Field_R_operator", ["float", "ptr", "Vec", "float"]),
        ("Field_R_operator", ["float", "ptr", "int", "Vec", "float"]),
        # ---------------------------------------------------------------------
        ("Field_C", ["ptr"]),
        ("Field_C", ["ptr", "string"]),
        ("Field_C_operator", ["complex<float>", "ptr", "float"]),
        ("Field_C_operator", ["complex<float>", "ptr", "int", "float"]),
        ("Field_C_operator", ["complex<float>", "ptr", "Vec", "float"]),
        ("Field_C_operator", ["complex<float>", "ptr", "int", "Vec", "float"]),
    ],
}


from write import write_export, write_import

write_export.write_export(INPUTS)
write_import.write_import(INPUTS)
