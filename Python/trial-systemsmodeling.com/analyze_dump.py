import json
import sys
from typing import Callable
from pprint import pprint

INDENT = "=="


def show_element(elem: dict, level: int = 0) -> None:
    elem_type = elem.get("@type")
    name = elem.get("declaredName", "NONE")
    shortName = elem.get("declaredShortName", "NONE")
    memberName = elem.get("memberName", "NONE")
    print(f"{level * INDENT}type={elem_type} id={elem['@id']} name={name} shortName={shortName} memberName={memberName}")
    pprint(elem)


def get_element(id_ref: dict[str, str], data_dict: dict[str, dict]):
    id = id_ref.get("@id")
    return data_dict[id]


def visit_relationship_tree(elem: dict, function: Callable, data_dict: dict[str, dict], level: int = 0) -> None:
    function(elem, level)
    ownedRelatedElement = elem.get("ownedRelatedElement")
    if ownedRelatedElement:
        for orl_ref in ownedRelatedElement:
            orl = get_element(orl_ref, data_dict)
            visit_relationship_tree(orl, function, data_dict, level + 1)
    for rel_ref in elem['ownedRelationship']:
        rel = get_element(rel_ref, data_dict)
        visit_relationship_tree(rel, function, data_dict, level + 1)


if __name__ == '__main__':
    with open("HPdK_MassPropertiesModelSmall_dump.json", 'r', encoding="utf-8") as json_file:
        data_dict = dict([(x["@id"], x) for x in json.load(json_file)])

        partUsages = [x for x in data_dict.values() if x["@type"] == "PartUsage"]
        # pprint(partUsages)

        pu20 = partUsages[20]
        visit_relationship_tree(pu20, show_element, data_dict)
