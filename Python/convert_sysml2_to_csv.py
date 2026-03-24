"""
Convert a SysML v2 model to a mass properties parts list in csv format
Input: SysML v2 model in contains a nested set of part mass property attributes
Output: csv file containing corresponding mass properties parts list

The part properties are:
    attribute POIconv : String = "-";
    attribute Ipoint : Boolean;
    attribute mass : Real;
    attribute sigma_mass : Real;
    attribute Cx : Real;
    attribute sigma_Cx : Real;
    attribute Cy : Real;
    attribute sigma_Cy : Real;
    attribute Cz : Real;
    attribute sigma_Cz : Real;
    attribute Ixx : Real;
    attribute sigma_Ixx : Real;
    attribute Iyy : Real;
    attribute sigma_Iyy : Real;
    attribute Izz : Real;
    attribute sigma_Izz : Real;
    attribute Ixy : Real;
    attribute sigma_Ixy : Real;
    attribute Ixz : Real;
    attribute sigma_Ixz : Real;
    attribute Iyz : Real;
    attribute sigma_Iyz : Real;

Attribute mass is the mass of a leaf part, or the composite mass of a non-leaf part.
Cx, Cy, Cz, Ixx, Iyy, Izz, Ixy, Ixz, Iyz are the components of the moments of inertia tensor.
The sigma attributes are associated uncertainties for all mass properties.
POIconv is xxx convention. Its value may be "-" or "+".
Ipoint is xxx (True) or yyy (False).
The SysML model should be a parts tree, with a single root part.
In the initial model only mass property values are given for the leaf parts.

Author: Hans Peter de Koning (DEKonsult)
"""

import syside
import sys
import os
import csv
from datetime import datetime, timezone
from typing import Optional

# Create logger for debug, info, warning, error, critical messages
import logging
LOGGER = logging.getLogger()
LOGGER.setLevel(logging.DEBUG)

class SysMLProcessor:
    def __init__(self):
        pass

    def convert_sysml_to_csv(self, sysml_model_path: str, csv_file_path: Optional[str]):
        LOGGER.debug("cwd=", os.getcwd())
        exists_sysml_model = os.path.exists(sysml_model_path)
        if not exists_sysml_model:
            LOGGER.error(f"Cannot find model: {sysml_model_path}")
            sys.exit(2)

        compiler = syside.Compiler()

        model, diag = syside.load_model([sysml_model_path])
        n_packages = len(tuple(model.elements(syside.Package)))
        n_parts = len(tuple(model.elements(syside.PartUsage)))
        LOGGER.info(f"model {sysml_model_path} has {n_packages} package(s) and {n_parts} parts")

        if csv_file_path:
            csv_path = csv_file_path
        else:
            root, ext = os.path.splitext(sysml_model_path)
            csv_path = f"{root}.csv"
        with open(csv_path, mode="w", encoding="utf8", newline="") as csv_file:

            # header = ["shortName", "name", "parentShortName", "POIconv", "Ipoint",
            #         "mass", "Cx", "Cx_sigma", "Cx", "Cx_sigma", "Cz", "Cz_sigma",
            #         "Ixx", "Ixx_sigma", "Ixy", "Ixy_sigma", "Ixz", "Ixz_sigma", "Iyy", "Iyy_sigma", "Iyz", "Iyz_sigma", "Izz", "Izz_sigma"]

            csv_writer = None
            is_first = True
            row_count = 0
            for elem in sorted(model.elements(syside.PartUsage), key=lambda x: x.short_name):

                row_dict: dict[str, any] = dict()

                row_dict["shortName"] = elem.short_name
                row_dict["name"] = elem.name

                owning_ns = elem.owning_namespace
                if row_count < 10:
                    LOGGER.debug("row_count={row_count} owning_ns={owning_ns}")
                if isinstance(owning_ns, syside.PartUsage):
                    row_dict["parentShortName"] = owning_ns.short_name
                else:
                    row_dict["parentShortName"] = "NA"

                for owned_elem in elem.owned_elements.collect():
                    if isinstance(owned_elem, syside.AttributeUsage):
                        feature_expression = owned_elem.feature_value_expression
                        if feature_expression:
                            val, diag = compiler.evaluate_feature(feature_expression, elem)
                            row_dict[owned_elem.name] = val
                        else:
                            row_dict[owned_elem.name] = "NA"

                if is_first:
                    LOGGER.debug("fieldnames={row_dict.keys()}")
                    csv_writer = csv.DictWriter(f=csv_file, fieldnames=row_dict.keys())
                    csv_writer.writeheader()
                    is_first = False
                    row_count += 1

                if row_count < 10:
                    LOGGER.debug("row#{row_count:5d} {row_dict}")
                csv_writer.writerow(row_dict)
                row_count += 1


if __name__ == "__main__":
    n_args = len(sys.argv)
    if  2 <= n_args <= 3:
        # Create message logger on console
        consoleHandler = logging.StreamHandler(sys.stdout)
        consoleHandler.setLevel(logging.INFO)
        consoleHandler.setFormatter(logging.Formatter("%(levelname)-8s: %(message)s"))
        LOGGER.addHandler(consoleHandler)

        start_time = datetime.now(timezone.utc)
        start_time_iso = start_time.isoformat(timespec="seconds").replace("+00:00", "Z")
        LOGGER.info(f"Run started at {start_time_iso}")

        sysml_model_path = sys.argv[1]
        csv_file_path = sys.argv[2] if n_args > 2 else None
        sysml_processor = SysMLProcessor()
        sysml_processor.convert_sysml_to_csv(sysml_model_path, csv_file_path)

        elapsed_duration = datetime.now(timezone.utc) - start_time
        LOGGER.info(f"Run completed in {elapsed_duration.total_seconds()} s")
    else:
        print(f"Usage: python convert_sysml2_to_csv.py SYSML_MODEL_PATH [CSV_FILE_PATH]")
        print()

