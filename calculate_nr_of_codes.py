"""
Calculate the number of possible genetic codes.

    n = the number of codons             (e.g. 64)
    k = the number of encoded outputs,   (e.g. for 20 aa and a Stop the value of k is 21)

So n is the number of "items" (codons) being assigned, and k is the number of
distinct "labels"

A genetic code is an assignment of every codon to one amino acid or Stop, using each
at least once. The count is the number of surjective (onto) functions from the n codons
onto the k outputs:

    codes(n, k) = k! * S(n, k)

where S(n, k) is the Stirling number of the second kind (partitions of n items
into k indistinguishable groups) and the k! factor labels those groups.
Everything is computed with exact (arbitrary-precision) integers and then shown
in scientific notation.

Recurrence used for S (same as the original Java program):
    S(n, k) = k * S(n-1, k) + S(n-1, k-1),  with S(0, 0) = 1.

Input : in_parameters.txt   (semicolon-separated: n ; k  ->  codons ; different outputs)
        The input file can list different several pairs of n and k, each pair per line.
        If this file does not exist, it is created automatically with the
        default pair (64 codons ; 21 outputs) shown above.
Output: calculated_nr_of_codes.txt   (n ; k ; possible codes  ->  codons ; encoded outputs ; codes)
Dr. Royal Truman 
"""

import math
import os
from decimal import Decimal

# --- Configuration ---------------------------------------------------------
INPUT_FILE = "in_parameters.txt"
OUTPUT_FILE = "calculated_nr_of_codes.txt"
DECIMALS = 10          # decimal places in the displayed scientific notation
SEP = ";"             # CSV field separator

# Default (n, k) pair used to seed the input file when it is missing.
# These are the example values shown in the module docstring:
# 64 codons and 21 outputs (to represent 20 aa + 1 stop).
DEFAULT_K = 21         # number of outputs
DEFAULT_N = 64         # number of codons

# Title row expected as the first line of the input file (it is skipped when read).
INPUT_HEADER = SEP.join(["number of codons, n", "number of encoded outputs"])


def ensure_input_file() -> None:
    """Make sure INPUT_FILE exists; if not, create it with the default values.

    The input file must be present for the program to read (n, k) pairs. If it
    is missing, a new file is written containing the title row plus one default
    pair (DEFAULT_N ; DEFAULT_K), and a message is displayed to inform the user.
    """
    if os.path.exists(INPUT_FILE):
        return

    default_row = SEP.join([str(DEFAULT_N), str(DEFAULT_K)])
    with open(INPUT_FILE, "w", encoding="utf-8") as f:
        f.write(INPUT_HEADER + "\n")
        f.write(default_row + "\n")

    print(
        f"NOTE: input file {INPUT_FILE!r} was not found, so it has been created "
        f"with the default values ({DEFAULT_N} codons ; {DEFAULT_K} outputs)."
    )


def stirling_second_kind(n: int, k: int) -> int:
    """Return S(n, k) exactly, using an O(n*k) row-by-row DP with a 1-D array."""
    if k > n:
        return 0
    prev = [0] * (k + 1)
    prev[0] = 1  # S(0, 0) = 1
    for i in range(1, n + 1):
        curr = [0] * (k + 1)
        for j in range(1, min(i, k) + 1):
            curr[j] = j * prev[j] + prev[j - 1]
        prev = curr
    return prev[k]


def genetic_codes(n: int, k: int) -> int:
    """Exact number of genetic codes = k! * S(n, k)  (n = codons, k = outputs)."""
    return math.factorial(k) * stirling_second_kind(n, k)


def to_scientific(value: int, decimals: int) -> str:
    """Format an arbitrary-precision integer in scientific notation, rounded."""
    return f"{Decimal(value):.{decimals}e}"


def validate(k: int, n: int) -> str | None:
    """Return an error message if the (k, n) pair is invalid, otherwise None."""
    if k < 1:
        return f"number of outputs ({k}) cannot be < 1"
    if n <= 1:
        return f"number of codons ({n}) must be > 1"
    if k > n:
        return (f"more outputs ({k}) than codons ({n}) "
                "is impossible")
    return None


def main() -> None:
    # Ensure the input file exists before attempting to read it; create it with
    # the default values (and notify the user) if it is missing.
    ensure_input_file()

    header = SEP.join(
        ["number of codons", "number of encoded outputs", "possible number of codes"]
    )
    print(header)

    output_rows = [header]

    with open(INPUT_FILE, encoding="utf-8") as f:
        next(f, None)  # skip the title header row
        # Each remaining line holds one n and k pair (codons ; encoded outputs); the
        # input file may contain a whole list of such pairs to be computed.
        for line_no, raw in enumerate(f, start=2):
            line = raw.strip()
            if not line:
                continue  # ignore blank lines

            parts = [p.strip() for p in line.split(SEP)]
            if len(parts) < 2:
                print(f"ERROR (line {line_no}): could not read a pair -> {line!r}")
                continue
            try:
                n = int(parts[0])
                k = int(parts[1])
            except ValueError:
                print(f"ERROR (line {line_no}): values are not integers -> {line!r}")
                continue

            error = validate(k, n)
            if error:
                print(f"ERROR (line {line_no}): {error}; row ignored")
                continue

            codes = genetic_codes(n, k)
            codes_sci = to_scientific(codes, DECIMALS)

            row = SEP.join([str(n), str(k), codes_sci])
            print(row)
            output_rows.append(row)

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write("\n".join(output_rows) + "\n")


if __name__ == "__main__":
    main()
