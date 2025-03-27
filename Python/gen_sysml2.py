"""
Generate SysML v2 model from mass properties nodes
Input: csv file containing mass properties nodes
Output: corresponding SysML v2 model in which each node is a part with mass property attributes

Author: Hans Peter de Koning (DEKonsult)
"""
import os
from typing import Optional, Union
from datetime import datetime, timezone
import csv
import sys

# Create logger for debug, info, warning, error, critical messages
import logging
LOGGER = logging.getLogger()
LOGGER.setLevel(logging.DEBUG)

INDENT = "    "

class Node:
    def __init__(self, short_name: str, name: str):
        self.short_name: str = short_name
        self.name: str = name
        self.children: list[Node] = list()
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
        return f"Node <{self.short_name}>"

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

class MassPropertiesModel:
    def __init__(self):
        self.root_node: Optional[Node] = None
        self.node_dict : dict[str, Node] = dict()

    def read_nodes(self, node_details_file_path: str):
        with (open(node_details_file_path, encoding="UTF8") as node_details_file):
            csv_reader = csv.reader(node_details_file, delimiter=',')
            header = next(csv_reader)

            LOGGER.debug(f"node_details_file header={header}")

            # Copy header list into attribute_columns and remove non-attribute columns
            attribute_columns: list[str] = header[:]
            attribute_columns.remove("id")
            attribute_columns.remove("name")
            attribute_columns.remove("parent")

            for raw_record in csv_reader:
                record = dict(zip(header, raw_record))
                # LOGGER.debug(f"record={record}")

                node_short_name = record["id"]
                if node_short_name in self.node_dict:
                    LOGGER.error(f"Non-unique node short_name: {node_short_name}. Node is already registered.")
                else:
                    node = Node(node_short_name, record["name"])
                    self.node_dict[node_short_name] = node

                    for attr_name in attribute_columns:
                        attr_kind = self.get_attribute_kind(attr_name)
                        node.attributes.append(Attribute(attr_kind, attr_name, record[attr_name]))

                    parent_short_name = record["parent"]
                    if parent_short_name == "NA":
                        if self.root_node:
                            LOGGER.error(f"{node} has parent=NA, but root node already exists: {self.root_node}. The model must contain only one root node.")
                        else:
                            self.root_node = node
                    else:
                        parent_node = self.node_dict.get(parent_short_name)
                        if parent_node:
                            parent_node.children.append(node)
                        else:
                            LOGGER.error(f"Cannot find parent node with short_name={parent_short_name}")

    @staticmethod
    def get_attribute_kind(attribute_name: str) -> str:
        attribute_kind_dict: dict[str, str] = {
            "POIconv": "String",
            "Ipoint": "Boolean",
        }
        return attribute_kind_dict.get(attribute_name, "Real")

    def report_model(self):
        LOGGER.info(f"Root node is {self.root_node}")
        LOGGER.info(f"Model contains {len(self.node_dict)} nodes")
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
        self.root_node.generate_sysml2(buffer)
        buffer.append("}")
        with open(sysml2_model_path, "w", encoding="UTF8") as sysml2_model:
            sysml2_model.write("".join(buffer))

    def composite_structure(self) -> str:
        buffer: list[str] = []
        for node in self.node_dict.values():
            if node.children:
                buffer.append(f"{INDENT}{node} contains {', '.join(str(child) for child in node.children)}")
            else:
                buffer.append(f"{INDENT}{node} is a leaf node")
        return "\n".join(buffer)

    def indented_composite_structure(self) -> str:
        buffer: list[str] = []
        self.add_composite_node(buffer, self.root_node)
        return "\n".join(buffer)

    def add_composite_node(self, buffer: list[str], node: Node, level: int = 1):
        leaf_node_marker = "" if node.children else " (leaf node)"
        buffer.append(f"{level*INDENT}{node}{leaf_node_marker}")
        # Walk composite tree depth-first
        for child in node.children:
            self.add_composite_node(buffer, child, level + 1)

# Main program
if __name__ == "__main__":
    # Create message logger on console
    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setLevel(logging.DEBUG)
    consoleHandler.setFormatter(logging.Formatter("%(levelname)-8s: %(message)s"))
    LOGGER.addHandler(consoleHandler)

    start_time = datetime.now(timezone.utc)
    start_time_iso = start_time.isoformat(timespec="seconds").replace("+00:00", "Z")
    LOGGER.info(f"Run started at {start_time_iso}")

    # nodes_file_path = "sm_table.csv"
    nodes_file_path = "sm_table_rollup.csv"
    mass_properties_model = MassPropertiesModel()
    mass_properties_model.read_nodes(nodes_file_path)
    mass_properties_model.report_model()
    # mass_properties_model.write_sysml2_model("MassPropertiesSmallModel.sysml")
    mass_properties_model.write_sysml2_model("MassPropertiesSmallModelRollup.sysml")

    duration = datetime.now(timezone.utc) - start_time
    LOGGER.info(f"Finished in {duration.total_seconds()} s")
