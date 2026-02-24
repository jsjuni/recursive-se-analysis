import syside
import os

def main():
    print("CWD=", os.getcwd())
    print("exists=", os.path.exists("MassPropertiesModel.sysml"))

    model, diag = syside.load_model(["MassPropertiesModel.sysml"])
    count = 0
    max_count = 10
    for elem in sorted(model.elements(syside.PartUsage), key=lambda x: x.short_name):
        print(f"sn={elem.short_name:20s} {elem.name} type={type(elem)}")
        count += 1
        if count > max_count:
            print(dir(elem))
            break

    sysml_package = list(model.elements(syside.Package))[0]
    print(sysml_package)
    print(f"sn={sysml_package.short_name:20s} {sysml_package.name} type={type(sysml_package)}")


    # json_string: str = syside.json.dumps(model, options=syside.SerializationOptions.minimal())
    # print("INFO: JSON serialization=")
    # print(json_string)

if __name__ == "__main__":
    main()

