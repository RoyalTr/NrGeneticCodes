# Genetic Code Counter

A small, dependency-free Python script that computes the **exact number of possible genetic codes** for a given number of codons and outputs (i.e., amino acid or Stop), using arbitrary-precision integer arithmetic. Results are printed to the console and written to a CSV file in scientific notation.

## What it computes

A *genetic code* here is an assignment of every codon to one output, where **every output is used at least once**. Since amino acids and Stop are distinct labels, this is the number of **surjective (onto) functions** from the set of codons onto the set of outputs:

```
codes(n, k) = k! * S(n, k)
```

where

- `n` = number of codons (e.g. 64)
- `k` = number of outputs (e.g. use 21 to represent 20 aa + 1 Stop)
- `S(n, k)` = the Stirling number of the second kind (the number of ways to partition `n` items into `k` non-empty, indistinguishable groups)

The `k!` factor labels those groups with the distinct outputs.

`S(n, k)` is built up from the recurrence

```
S(n, k) = k * S(n-1, k) + S(n-1, k-1),   with S(0, 0) = 1
```

using an `O(n·k)` row-by-row dynamic program with a one-dimensional array. All intermediate values are exact Python integers (no floating-point loss); only the final display is converted to scientific notation.

## Requirements

- Python 3.10 or newer (the code uses the `X | None` type-hint syntax)
- No third-party packages — only the standard library (`math`, `os`, `decimal`)

## Usage

```bash
python calculate_nr_of_codes.py
```

The script reads its parameters from `in_parameters.txt` and writes results to `calculated_nr_of_codes.txt`, while also printing them to the console.

**If `in_parameters.txt` does not exist**, the script creates it automatically with a default pair (64 codons; 21 outputs), prints a note saying it did so, and then runs on that default.

## Input file format

`in_parameters.txt` is a semicolon-separated file. The **first line is a header** and is skipped. Each remaining line is one `codons ; outputs` pair. You can list as many pairs as you like, one per line.

```
number of codons;number of encoded outputs
64;21
64;10
16;10
```

> **Column order:** the first column is the number of **codons** (`n`), the second is the number of **outputs** (`k`) — matching the order of the parameters in the Stirling-number formula.

## Output

Results are written to `calculated_nr_of_codes.txt` (and echoed to the console) with a header and one row per valid input pair:

```
number of codons;number of encoded outputs;possible number of codes
64;21;1.5101095158e+84
64;20;8.1492025012e+82
10;4;8.1852000000e+5
```

The number of decimal places in the scientific notation is controlled by the `DECIMALS` constant (default `10`).

## Validation

Each input pair is checked before it is computed. A row is reported as an error and skipped if:

- the number of coded outputs is not >= 1,
- the number of codons is not >= 1, or
- there are more coded outputs than codons.

Lines that cannot be parsed as two integers are also reported and skipped, so one bad row does not stop the rest of the file from being processed.

## Configuration

A few constants near the top of the script can be adjusted:

| Constant | Default | Meaning |
| --- | --- | --- |
| `INPUT_FILE` | `in_parameters.txt` | File the script reads pairs from |
| `OUTPUT_FILE` | `calculated_nr_of_codes.txt` | File the results are written to |
| `DECIMALS` | `10` | Decimal places in the scientific-notation output |
| `SEP` | `;` | Field separator used in both files |
| `DEFAULT_N` | `64` | Default number of codons used to seed a missing input file |
| `DEFAULT_K` | `21` | Default number of coded outputs used to seed a missing input file |

## Example

The standard biological case — 64 codons mapped onto 21 labels (20 amino acids plus a stop signal) — yields roughly **1.51 × 10⁸⁴** possible codes.

## Background

This is a Python reimplementation of an earlier Java program; it uses the same Stirling-number recurrence but relies on Python's built-in big integers for exact results of any size.
