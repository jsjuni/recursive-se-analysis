import syside
import os

def main():
    model, diag = syside.load_model(["NegativeNumbersIssue.sysml"])
    count = 0
    for elem in model.elements(syside.LiteralRational, include_subtypes=True):
        # print(f"sn={elem.short_name:20s} {elem.name} type={type(elem)}")
        name = elem.name if elem.name else "UNNAMED"
        print(f"name={name:30s} type={type(elem).__name__}")
        print(elem)
        count += 1

if __name__ == "__main__":
    main()
