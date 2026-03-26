"""
Construct a SysML v2 model from mass properties parts list
Input: csv file containing a mass properties parts list
Output: corresponding SysML v2 model in which each part is a part with mass property attributes

Author: Hans Peter de Koning (DEKonsult)
"""
import os
from collections import namedtuple
from typing import Optional, Union
from datetime import datetime, timezone
import csv
import sys

# Create logger for debug, info, warning, error, critical messages
import logging
LOGGER = logging.getLogger()
LOGGER.setLevel(logging.DEBUG)

INDENT = 4*" "

class Part:
    def __init__(self, short_name: str, name: str):
        self.short_name: str = short_name
        self.name: str = name
        self.children: list[Part] = list()
        self.attributes: list[Attribute] = list()

    def generate_sysml2(self, buffer: list[str], level: int = 1):
        buffer.append(f"{level*INDENT}part <'{self.short_name}'> '{self.name}' {{\n")
        for attribute in self.attributes:
            buffer.append(f"{(level+1)*INDENT}{attribute.get_sysml2_declaration()}\n")

        for child in self.children:
            child.generate_sysml2(buffer, level + 1)

        # Add closing curly bracket
        buffer.append(f"{level*INDENT}}}\n")

    def __repr__(self):
        return f"Part <{self.short_name}>"

class Attribute:
    def __init__(self, kind: str, name: str, value: Optional[str] = None):
        self.kind: str = kind
        self.name: str = name
        self.value: Union[bool, float, str, None] = self.get_python_value(kind, value)

    @staticmethod
    def get_python_value(kind: str, value: Optional[str]) -> Union[bool, float, str, None]:
        python_value = None
        if value == "NA":
            pass
        else:
            if kind == "Boolean":
                if value.upper() == "TRUE":
                    python_value = True
                elif value.upper() == "FALSE":
                    python_value = False
                else:
                    LOGGER.error(f"{value} is not a valid Boolean value")
            elif kind == "Real":
                try:
                    python_value = float(value)
                except ValueError:
                    LOGGER.error(f"{value} is not a valid Real value")
            elif kind == "String":
                python_value = value
            else:
                LOGGER.error(f"unsupported kind: {kind}")
        return python_value

    def get_sysml2_declaration(self) -> str:
        attribute_name = self.name if self.name.isascii() else f"'{self.name}'"
        value = ""
        if self.value is not None:
            if self.kind == "Boolean":
                value = " = true" if self.value else " = false"
            elif self.kind == "Real":
                value = f" = {self.value}"
            elif self.kind == "String":
                value = f" = \"{self.value}\""
        return f"attribute {attribute_name} : {self.kind}{value};"

class SysMLModel:
    def __init__(self):
        self.root_part: Optional[Part] = None
        self.part_dict : dict[str, Part] = dict()

    def read_parts(self, part_details_file_path: str):
        (basename, ext) = os.path.splitext(part_details_file_path)
        delimiter = ""
        if ext == ".tsv":
            delimiter = "\t"
        elif ext == ".csv":
            delimiter = ","
        else:
            LOGGER.error(f"Wrong file type: {ext}. Must be .csv or .tsv format.")

        with (open(part_details_file_path, encoding="UTF8") as part_details_file):
            csv_reader = csv.DictReader(part_details_file, delimiter=delimiter)
            header = csv_reader.fieldnames

            LOGGER.debug(f"part_details_file header={header}")

            has_parent = "parent" in header
            has_pid = "pid" in header

            # Copy header list into attribute_columns and remove non-attribute columns
            attribute_columns: list[str] = header[:]
            for attr_name in ("id", "key", "name", "parent", "pid"):
                try:
                    attribute_columns.remove(attr_name)
                except Exception as e:
                    pass
            LOGGER.debug(f"attribute_columns={attribute_columns}")

            # Reorder attributes if not already in mass, sigma_mass order
            if "mass" in attribute_columns:
                i = attribute_columns.index("mass")
                if i + 1 < len(attribute_columns) and attribute_columns[i+1] != "sigma_mass":
                    reindexed = []
                    half_len = len(attribute_columns)//2

                    for i in range(half_len):
                        reindexed.append(i)
                        reindexed.append(i+half_len)
                    attribute_columns = [attribute_columns[i] for i in reindexed]
                    LOGGER.debug(f"attribute_columns={attribute_columns}")

            for record in csv_reader:
                # LOGGER.debug(f"record={record}")

                part_short_name = record["id"]
                if part_short_name in self.part_dict:
                    LOGGER.error(f"Non-unique part short_name: {part_short_name}. Part is already registered.")
                else:
                    part = Part(part_short_name, record["name"])
                    self.part_dict[part_short_name] = part

                    for attr_name in attribute_columns:
                        attr_kind = self.get_attribute_kind(attr_name)
                        part.attributes.append(Attribute(attr_kind, attr_name, record[attr_name]))

                    if has_parent:
                        parent_short_name = record["parent"]
                    elif has_pid:
                        parent_short_name = record["pid"]
                    else:
                        # csv file has no parent column - construct parent
                        level = part_short_name.count(".")
                        if level == 1:
                            parent_short_name = "NA"
                        else:
                            parent_short_name = part_short_name.rsplit(".", 1)[0]
                    if parent_short_name == "NA":
                        if self.root_part:
                            LOGGER.error(f"{part} has parent=NA, but root part already exists: {self.root_part}. The model must contain only one root part.")
                        else:
                            self.root_part = part
                    else:
                        parent_part = self.part_dict.get(parent_short_name)
                        if parent_part:
                            parent_part.children.append(part)
                        else:
                            LOGGER.error(f"Cannot find parent part with short_name={parent_short_name}")

    @staticmethod
    def get_attribute_kind(attribute_name: str) -> str:
        attribute_kind_dict: dict[str, str] = {
            "POIconv": "String",
            "Ipoint": "Boolean",
        }
        return attribute_kind_dict.get(attribute_name, "Real")

    def report_model(self):
        LOGGER.info(f"Root part is {self.root_part}")
        LOGGER.info(f"Model contains {len(self.part_dict)} parts")
        # LOGGER.info(f"composite structure:\n{self.composite_structure()}")
        LOGGER.info(f"Indented composite structure:\n{self.indented_composite_structure()}")

    def write_sysml2_model(self, sysml2_model_path: str):
        LOGGER.info(f"Generating SysML2 model in {sysml2_model_path}")
        model_name = os.path.basename(os.path.splitext(sysml2_model_path)[0])

        buffer: list[str] = []
        buffer.append(f"package {model_name} {{\n")
        buffer.append(f"{INDENT}private import ScalarValues::Boolean;\n")
        buffer.append(f"{INDENT}private import ScalarValues::Real;\n")
        buffer.append(f"{INDENT}private import ScalarValues::String;\n")
        buffer.append("\n")
        self.root_part.generate_sysml2(buffer)
        buffer.append("}")
        with open(sysml2_model_path, "w", encoding="UTF8") as sysml2_model:
            sysml2_model.write("".join(buffer))

    def composite_structure(self) -> str:
        buffer: list[str] = []
        for part in self.part_dict.values():
            if part.children:
                buffer.append(f"{INDENT}{part} contains {', '.join(str(child) for child in part.children)}")
            else:
                buffer.append(f"{INDENT}{part} is a leaf part")
        return "\n".join(buffer)

    def indented_composite_structure(self) -> str:
        buffer: list[str] = []
        self.add_composite_part(buffer, self.root_part)
        return "\n".join(buffer)

    def add_composite_part(self, buffer: list[str], part: Part, level: int = 1):
        leaf_part_marker = "" if part.children else " (leaf part)"
        buffer.append(f"{level*INDENT}{part}{leaf_part_marker}")
        # Walk composite tree depth-first
        for child in part.children:
            self.add_composite_part(buffer, child, level + 1)

# Main program
if __name__ == "__main__":
    n_args = len(sys.argv)
    if  2 <= n_args <= 3:
        # Create message logger on console
        consoleHandler = logging.StreamHandler(sys.stdout)
        consoleHandler.setLevel(logging.DEBUG)
        consoleHandler.setFormatter(logging.Formatter("%(levelname)-8s: %(message)s"))
        LOGGER.addHandler(consoleHandler)

        start_time = datetime.now(timezone.utc)
        start_time_iso = start_time.isoformat(timespec="seconds").replace("+00:00", "Z")
        LOGGER.info(f"Run started at {start_time_iso}")

        csv_file_path = sys.argv[1]
        if n_args > 2:
            sysml_model_path = sys.argv[2]
        else:
            root, ext = os.path.splitext(csv_file_path)
            sysml_model_path = f"{root}.sysml"

        # Generate a mass properties model
        mass_properties_model = SysMLModel()
        mass_properties_model.read_parts(csv_file_path)
        mass_properties_model.report_model()
        mass_properties_model.write_sysml2_model(sysml_model_path)

        elapsed_duration = datetime.now(timezone.utc) - start_time
        LOGGER.info(f"Run completed in {elapsed_duration.total_seconds()} s")
    else:
        print(f"Usage: python convert_csv_to_sysml2.py CSV_FILE_PATH [SYSML_MODEL_PATH]")
        print()
