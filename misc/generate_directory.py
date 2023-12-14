import time
from sys import argv
import os
from typing import TextIO


def translate_file(f_in: TextIO, f_out: TextIO, target: str, replacement: str, folder: str, tran_app: str, tran_app_cap: str) -> None:
    if not f_out.name.endswith(".txt"):
        f_out.write(f"// Automatically generated from {f_in.name} on {time.ctime(time.time())}\n")

    for line in f_in:
        new_line = ((line.replace(target, replacement).replace("jmtx_", f"jmtx{tran_app}_").
                    replace("jmtxs_", f"jmtx{tran_app}s_")).
                    replace("JMTX_", f"JMTX{tran_app_cap}_").
                    replace(f"jmtx{tran_app}_result", "jmtx_result").
                    replace(f"JMTX{tran_app_cap}_RESULT", "JMTX_RESULT").
                    replace(f"JMTX{tran_app_cap}_DEFAULT", "JMTX_DEFAULT").
                    replace(f"JMTX{tran_app_cap}_NODISCARD", "JMTX_NODISCARD").
                    replace(f"jmtx{tran_app}_allocator_callbacks", "jmtx_allocator_callbacks").
                    replace(f"jmtx{tran_app}_allocator_callbacks", "jmtx_allocator_callbacks").
                    replace(f"jmtx{tran_app}_matrix_base", "jmtx_matrix_base").
                    replace(f"jmtx{tran_app}_internal_find_last_leq_value", f"jmtx_internal_find_last_leq_value").
                    replace(f"\"{replacement}", f"\"{folder}").
                    replace(f"\\{replacement}", f"\\{folder}").
                    replace(f"/{replacement}", f"/{folder}").
                    replace(f"_{replacement} ", f"_{folder} ").
                    replace(f"<{replacement}.h>", f"<{target}.h>").
                    replace(f"jmtx{tran_app}_matrix_type_to_str", f"jmtx_matrix_type_to_str").
                    replace(f"sqrtf", "sqrt").
                    replace(f"crealf", "creal").
                    replace(f"cimagf", "cimag").
                    replace(f"cabsf", "cabs").
                    replace("fabsf", "fabs"))
        f_out.write(new_line)

def process_directory(dirname: str, out_dir: str, tgt: str, rpl: str, folder: str, tran_app: str, tran_app_cap: str) -> None:
    assert os.path.isdir(dirname)
    for entry in os.listdir(dirname):
        full_name = os.path.join(dirname, entry)
        out_name = os.path.join(out_dir, entry)
        if os.path.isdir(full_name):
            if not os.path.isdir(out_name):
                os.mkdir(out_name)
            process_directory(full_name, out_name, tgt, rpl, folder, tran_app, tran_app_cap)
        else:
            f_in = open(full_name, "r")
            f_out = open(out_name, "w")
            translate_file(f_in, f_out, tgt, rpl, folder, tran_app, tran_app_cap)
            f_in.close()
            f_out.close()
            print(f"Translated \"{full_name}\" to \"{out_name}\"")

if __name__ == "__main__":
    if len(argv) != 4:
        print(f"3 arguments are required {len(argv)} were given (argv = {argv})")
        exit(1)

    input_directory = argv[1]
    output_directory = argv[2]
    target_type = argv[3]

    if not os.path.isdir(input_directory):
        print(f"\"{input_directory}\" is not a directory!")
        exit(1)

    translated_initials = {"_Complex float": "c", "_Complex double": "z", "double" : "d"}
    translated_initials_CAP = {"_Complex float": "C", "_Complex double": "Z", "double" : "D"}
    translated_folder = {"_Complex float": "cfloat", "_Complex double": "cdouble", "double" : "double"}
    if target_type not in translated_initials.keys():
        print(f"\"{target_type}\" was not one of the allowed types ({translated_initials.keys()})")
        exit(1)
    out_path = os.path.join(input_directory, "..", output_directory)

    if os.path.isfile(out_path):
        print(f"\"{out_path}\" is a file!")
        exit(1)

    if not os.path.isdir(out_path):
        os.mkdir(out_path)

    process_directory(input_directory, out_path, "float", target_type, output_directory,
                      translated_initials[target_type], translated_initials_CAP[target_type])

