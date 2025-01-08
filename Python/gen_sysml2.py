from typing import Union
import csv

INDENT = "    "


class Attribute:
    def __init__(self, kind: str, name: str, value: str):
        self.kind = kind
        self.name = name
        self.value = value

class Node:
    def __init__(self, short_name):
        self.short_name: str = short_name
        self.name: str = ""
        self.children: list[Node] = list()
        self.attributes: list[Attribute] = list()

    def add_attribute(self, kind: str, name: str, value: Union[str, None]):
        self.attributes.append(Attribute(kind, name, value))

    def generate_sysml2(self, buffer: list[str], level: int = 1):
        buffer.append(f"{level*INDENT}part <'{self.short_name}'> '{self.name}'")
        buffer.append(f" {{\n")
        for attribute in self.attributes:
            attribute_name = attribute.name
            if not attribute_name.isascii():
                attribute_name = f"'{attribute_name}'"
            value = "" if attribute.value == "NA" else f" = {attribute.value}"
            buffer.append(f"{(level+1)*INDENT}attribute {attribute_name} : {attribute.kind}{value};\n")
        if self.children:
            for child in self.children:
                child.generate_sysml2(buffer, level + 1)
        buffer.append(f"{level*INDENT}}}\n")

    def __repr__(self):
        return f'<Node: {self.short_name}>'

class MassPropertiesModel:
    def __init__(self):
        self.root_node: Union[Node, None] = None
        self.node_dict : dict[str, Node] = dict()

    def read_edge_list(self, edge_list_file_path: str):
        root_node = None
        with open(edge_list_file_path, encoding="UTF8") as edge_list_file:
            csv_reader = csv.reader(edge_list_file, delimiter=',')
            header = next(csv_reader)
            print("edge_list_file header=", header)

            for edge in csv_reader:
                parent_label = edge[0]
                child_label = edge[1]
                parent_node = self.node_dict.get(parent_label)
                if not parent_node:
                    parent_node = Node(parent_label)
                    self.node_dict[parent_label] = parent_node
                if not self.root_node:
                    self.root_node = parent_node
                child_node = self.node_dict.get(child_label)
                if not child_node:
                    child_node = Node(child_label)
                    self.node_dict[child_label] = child_node
                parent_node.children.append(child_node)

            print("DEBUG: node graph:")
            for key, item in self.node_dict.items():
                print(f"key={key} children={item.children}")

    def read_node_details(self, node_details_file_path: str):
        root_node = None
        with open(node_details_file_path, encoding="UTF8") as node_details_file:
            csv_reader = csv.reader(node_details_file, delimiter='\t')
            header = next(csv_reader)
            print("node_details_file header=", header)

            for record in csv_reader:
                map = dict(zip(header, record))
                # print(f"map={map}")
                node = self.node_dict.get(map["id"])
                if node:
                    node.name = map["name"]
                    node.add_attribute("String", 'POIconv', f"\"{map['POIconv']}\"")
                    value = "NA" if map['Ipoint'] == "NA" else map['Ipoint'].lower()
                    node.add_attribute("Boolean", 'Ipoint', value)
                    for attr_name in ['mass', 'σ_mass', 'Cx', 'σ_Cx', 'Cy', 'σ_Cy', 'Cz', 'σ_Cz', 'Ixx', 'σ_Ixx', 'Iyy', 'σ_Iyy', 'Izz', 'σ_Izz', 'Ixy', 'σ_Ixy', 'Ixz', 'σ_Ixz', 'Iyz', 'σ_Iyz']:
                        node.add_attribute("Real", attr_name, map[attr_name])

    def write_sysml2_model(self, sysml2_model_path: str):
        print(f"=== mass properties model in SysML2 ===")
        print(f"root node={self.root_node}")
        print(f"model contains {len(self.node_dict.keys())} nodes")
        buffer: list[str] = list()
        buffer.append("package MassPropertiesModel {\n")
        buffer.append("    private import ScalarValues::Boolean;\n")
        buffer.append("    private import ScalarValues::Real;\n")
        buffer.append("    private import ScalarValues::String;\n")
        buffer.append("\n")
        self.root_node.generate_sysml2(buffer)
        buffer.append("\n")
        buffer.append("}")
        with open(sysml2_model_path, "w", encoding="UTF8") as sysml2_model:
            sysml2_model.write("".join(buffer))

if __name__ == "__main__":
    mass_properties_model = MassPropertiesModel()
    mass_properties_model.read_edge_list("../el-ok.csv")
    mass_properties_model.read_node_details("../mp-input.tsv")

    mass_properties_model.write_sysml2_model("MassPropertiesModel.sysml")
