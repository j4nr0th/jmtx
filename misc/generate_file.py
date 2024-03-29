import os
from sys import argv
from generate_directory import translate_file


if __name__ == "__main__":
    if len(argv) != 4:
        print(f"3 arguments are required {len(argv)} were given (argv = {argv})")
        exit(1)

    input_file = argv[1]
    output_file = argv[2]
    target_type = argv[3]

    if not os.path.isfile(input_file):
        print(f"\"{input_file}\" is not a file!")
        exit(1)

    translated_initials = {"_Complex float": "c", "_Complex double": "z", "double" : "d"}
    translated_initials_CAP = {"_Complex float": "C", "_Complex double": "Z", "double" : "D"}
    translated_folder = {"_Complex float": "cfloat", "_Complex double": "cdouble", "double" : "double"}
    if target_type not in translated_initials.keys():
        print(f"\"{target_type}\" was not one of the allowed types ({translated_initials.keys()})")
        exit(1)


    if os.path.exists(output_file):
        print(f"\"{output_file}\" already exists!")
        exit(1)

    f_in = open(input_file, "r")
    f_out = open(output_file, "w")
    translate_file(f_in, f_out, "float", target_type, translated_folder[target_type],
                   translated_initials[target_type], translated_initials_CAP[target_type])

