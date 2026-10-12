#!/usr/bin/env python3

import sys
from pathlib import Path


METRICS = {
    "jaccard": "Jaccard",
    "bray-curtis": "Bray--Curtis",
    "aitchison": "Aitchison",
    "unweighted-unifrac": "Unweighted UniFrac",
    "weighted-unifrac": "Weighted UniFrac",
}


def parse_file(filename):
    """
    Read PERMANOVA results produced by permanova.R.

    Expected format:

        R2      F       p
        0.52075 1.0975  0.008
    """

    lines = Path(filename).read_text().strip().splitlines()

    if len(lines) < 2:
        raise ValueError(
            f"PERMANOVA result file is empty or incomplete: {filename}"
        )

    # The first line should contain the column names:
    # R2    F    p
    header = lines[0].split()

    if header != ["R2", "F", "p"]:
        raise ValueError(
            f"Unexpected header in {filename}: {lines[0]}"
        )

    # The second line contains the values.
    values = lines[1].split()

    if len(values) != 3:
        raise ValueError(
            f"Expected three values (R2, F, p) in {filename}, "
            f"got: {lines[1]}"
        )

    try:
        r2 = float(values[0])
        f_value = float(values[1])
        p_value = float(values[2])
    except ValueError as e:
        raise ValueError(
            f"Cannot parse PERMANOVA values in {filename}: "
            f"{lines[1]}"
        ) from e

    return r2, f_value, p_value


def metric_name(filename):
    name = Path(filename).name
    name = name.removesuffix("-matrix-permanova.txt")

    return METRICS.get(name, name)


def main():
    if len(sys.argv) < 3:
        raise ValueError(
            "At least two arguments are necessary "
            "(last argument must be output file)."
        )

    files = sys.argv[1:-1]
    output_file = sys.argv[-1]

    if not files:
        print(
            "Usage: python3 permanova-to-latex.py "
            "<metric>-matrix-permanova.txt [...]",
            file=sys.stderr,
        )
        sys.exit(1)

    results = []

    for filename in files:
        try:
            r2, f_value, p_value = parse_file(filename)

            results.append({
                "metric": metric_name(filename),
                "r2": r2,
                "f": f_value,
                "p": p_value,
            })

        except ValueError as e:
            print(f"WARNING: {e}", file=sys.stderr)

    if not results:
        print(
            "ERROR: No valid PERMANOVA results found.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Generate LaTeX
    latex = []

    latex.append(r"\begin{table}[ht]")
    latex.append(r"\centering")
    latex.append(
        r"\caption{PERMANOVA results for the effect of tissue type.}"
    )
    latex.append(r"\label{tab:permanova}")
    latex.append(r"\begin{tabular}{lccc}")
    latex.append(r"\hline")
    latex.append(
        r"Distance metric & $R^2$ & Pseudo-$F$ & $p$ \\"
    )
    latex.append(r"\hline")

    for result in results:
        latex.append(
            f"{result['metric']} & "
            f"{result['r2']:.3f} & "
            f"{result['f']:.3f} & "
            f"{result['p']:.3f} \\\\"
        )

    latex.append(r"\hline")
    latex.append(r"\end{tabular}")
    latex.append(r"\end{table}")

    Path(output_file).write_text(
        "\n".join(latex) + "\n"
    )


if __name__ == "__main__":
    main()
