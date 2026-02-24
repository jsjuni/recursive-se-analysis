package NegativeNumbersIssue {
    private import ScalarValues::*;
    part def Part1 {
        doc /* In the concrete textual notation negative numbers 
             * are represented as a unary minus function 
             * applied to an unsignend literal number */
        attribute att1 : Real = 1.2;
        attribute att2 : Real = -1.2;
    }
}
