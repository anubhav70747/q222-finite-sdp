"""Build the three self-contained LaTeX documents with pdfLaTeX."""
from pathlib import Path
import shutil, subprocess
root=Path(__file__).resolve().parent
for part in ("note","article","proofs"):
    out=root/"build"/part
    out.mkdir(parents=True,exist_ok=True)
    for _ in range(2):
        run=subprocess.run(["pdflatex","-interaction=nonstopmode","-halt-on-error",
            "-output-directory="+str(out),"main.tex"],cwd=root/part,
            capture_output=True,text=True,encoding="utf-8",errors="replace")
        if run.returncode:
            print(run.stdout[-6000:])
            raise SystemExit(run.returncode)
    log=(out/"main.log").read_text(errors="replace")
    for bad in ("undefined references","undefined on input","Overfull", "duplicate ignored"):
        if bad in log:
            raise RuntimeError(f"{part}: {bad}; inspect build/{part}/main.log")
    shutil.copy2(out/"main.pdf",root/part/"main.pdf")
    print("Built",part+"/main.pdf")
