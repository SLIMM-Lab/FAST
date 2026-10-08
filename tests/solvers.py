#!/usr/bin/env python3
"""Exercise existing example models, solver selection, and PETSc failures."""
import argparse
import math
from pathlib import Path
import re
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("executable", type=Path)
parser.add_argument("--petsc", action="store_true")
args = parser.parse_args()
executable = args.executable.resolve()
repo = Path(__file__).resolve().parents[1]
checks = 0


def run(source, settings="", cli=(), files=None, error=None):
    global checks
    with tempfile.TemporaryDirectory(prefix="fast-solvers-") as temporary:
        directory = Path(temporary)
        model = directory / "model"
        data = (repo / "dat" / (source + ".dat")).read_text()
        # Keep the original header and place settings before the first node.
        data = data.replace("%NODE\n", settings + "\n%NODE\n", 1)
        model.with_suffix(".dat").write_text(data)
        for name, content in (files or {}).items():
            (directory / name).write_text(content)
        completed = subprocess.run(
            [str(executable), str(model), *cli, "-silent"],
            cwd=repo, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=60)
        if error:
            assert completed.returncode != 0, completed.stdout
            assert error in completed.stdout, completed.stdout
            output = model.with_suffix(".pos").read_text()
            assert not output.rstrip().endswith("%END"), output
            checks += 1
            return
        assert completed.returncode == 0, completed.stdout
        output = model.with_suffix(".pos").read_text()
        assert output.rstrip().endswith("%END"), completed.stdout
        assert "%RESULT.CASE.STEP" in output, output
        checks += 1
        return output


def node_values(output, marker, node):
    section = output.split(marker + "\n", 1)[1].split("\n%", 1)[0]
    row = next(line.split() for line in section.splitlines()[1:]
               if line.split() and line.split()[0] == str(node))
    return list(map(float, row[1:]))


def close(actual, expected):
    assert len(actual) == len(expected), (actual, expected)
    assert all(math.isfinite(a) and math.isclose(a, e, rel_tol=2e-5, abs_tol=1e-9)
               for a, e in zip(actual, expected)), (actual, expected)


frame = "frame/LinFrame"
displacement = "%RESULT.CASE.STEP.NODAL.DISPLACEMENT"
baseline = run(frame)
expected = node_values(baseline, displacement, 2)
close(expected, [2.64122e-3, -2.10844e-3, 0, 0, 0, -1.30368e-4])

# All original storage/solver choices must continue to work.
for storage in ("skyline", "unsym_skyline", "sparse", "unsym_sparse", "eigen_unsym_sparse"):
    output = run(frame, "%HEADER.ANALYSIS.MATRIX.TYPE\n'" + storage +
                 "'\n%HEADER.ANALYSIS.MAXIMUM.ITERATIONS\n200\n"
                 "%HEADER.ANALYSIS.TOLERANCE\n1e-12\n")
    close(node_values(output, displacement, 2), expected)

continuum = run("continuum/TestQ4")
close(node_values(continuum, displacement, 4)[:2], [5e-6, -1.5e-6])
heat = run("thermal/TestT3")
close(node_values(heat, "%RESULT.CASE.STEP.NODAL.SCALAR.DATA", 4), [-1.78824])

select = "%HEADER.ANALYSIS.LINEAR.SOLVER\n'petsc'\n"
option_marker = "%HEADER.ANALYSIS.PETSC.OPTIONS\n"
if not args.petsc:
    run(frame, select, error="rebuild with USE_PETSC=1")
else:
    # PETSc methods are passed through, without a FAST-side whitelist.
    for options in ("-ksp_type gmres -pc_type ilu", "-ksp_type cg -pc_type jacobi",
                    "-ksp_type bcgs -pc_type jacobi", "-ksp_type preonly -pc_type lu"):
        output = run(frame, select + option_marker + "'" + options + " -ksp_rtol 1e-12'\n")
        close(node_values(output, displacement, 2), expected)
    output = run(frame, select)  # documented GMRES/ILU defaults
    close(node_values(output, displacement, 2), expected)

    # Model-relative file path works while running from a different directory.
    # CLI must override inline settings, which must override the options file.
    settings = (select + "%HEADER.ANALYSIS.PETSC.OPTIONS.FILE\n'solver.opts'\n" +
                option_marker + "'-ksp_type gmres -pc_type ilu -ksp_rtol 1e-12'\n")
    output = run(frame, settings, files={"solver.opts": "-ksp_type invalid_file_type\n"})
    close(node_values(output, displacement, 2), expected)
    output = run(frame, select + option_marker + "'-ksp_type invalid_inline_type'\n",
                 cli=("-ksp_type", "preonly", "-pc_type", "lu"))
    close(node_values(output, displacement, 2), expected)

    # Check thermal, continuum, a skyline-only dynamic controller, and nonlinear
    # reassembly against the native results, including every printed step.
    for source in ("continuum/TestQ4", "thermal/TestT3", "dynamic/MassSpring", "nonlin/ShallowTruss"):
        native = run(source)
        petsc = run(source, select + option_marker + "'-ksp_type preonly -pc_type lu'\n")
        def result_numbers(text):
            tail = text.split("%RESULT.CASE.STEP\n", 1)[1]
            return [float(n) for n in re.findall(r"[+-]?\d+\.\d+(?:e[+-]?\d+)?", tail, re.I)]
        close(result_numbers(petsc), result_numbers(native))

    run(frame, select + option_marker + "'-ksp_type invalid_solver'\n", error="PETSc error")
    run(frame, select + option_marker + "'-pc_type invalid_pc'\n", error="PETSc error")
    run(frame, select + "%HEADER.ANALYSIS.PETSC.OPTIONS.FILE\n'missing.opts'\n", error="PETSc error")
    run(frame, select + option_marker + "'-ksp_type gmres -pc_type none -ksp_max_it 1 -ksp_rtol 1e-14'\n",
        error="did not converge")

run(frame, "%HEADER.ANALYSIS.LINEAR.SOLVER\n'unknown'\n", error="Unknown linear solver")
run(frame, option_marker + "'-ksp_type cg'\n", error="PETSc options require")
run(frame, "%HEADER.ANALYSIS.LINEAR.SOLVER\nnative\n", error="single-quoted string")
run(frame, "%HEADER.ANALYSIS.LINEAR.SOLVER\n'petsc\n", error="single-quoted string")
print(f"PASS: {checks} model runs and error checks ({'PETSc + native' if args.petsc else 'native'})")
