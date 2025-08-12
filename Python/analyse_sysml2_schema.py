from typing import Optional
import json


if __name__ == "__main__":
    with open("openapi.json", "r", encoding="UTF8") as json_schema_file:
        json_schema = json.load(json_schema_file)

    schemas = json_schema["components"]["schemas"]
    for schema_name, schema in schemas.items():
        print(f"Schema {schema_name}:")
        print(f"  title: {schema.get('title')}")
        print(f"  description: {schema.get('description')}")
        print(f"  type: {schema.get('type')}")
        print(f"  required_properties: {schema.get('required')}")
        print(f"  properties:")
        for property_name, property in schema.get("properties", {}).items():
            type_name = property.get("type")
            items_type_str = ""
            if type_name == "array":
                items = property.get("items")
                items_type = items.get("type")
                if items_type is None:
                    items_type = items.get("$ref")
                items_type_str = f" of {items_type}"
            print(f"    {property_name}:")
            print(f"      type: {type_name}{items_type_str}")
